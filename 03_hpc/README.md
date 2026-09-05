# 03_hpc — HPC Pipeline Documentation

This directory contains the job-array-generation (using YCRC's dSQ module) notebooks, R execution scripts, and
summary/compilation tooling for running brms models on the Yale Bouchet HPC cluster.

---

## Directory Overview

Every file listed below is present in this directory. Job generation runs
locally, the R scripts run on the cluster (one model per SLURM array task), and
the compile scripts pull the results back.

### Job generation (run locally)

| File | Role |
|---|---|
| **`generate_hpc_jobs.py`** | **Primary entry point.** Reads the shipped analysis dataframe, writes `df_foranalysis_master.csv` (**raw — it does NOT normalize**, see § Normalization) plus the auto-generated `categorical_factor_vars.R` / `monotonic_covariates_vars.R` / `normalization_vars.R`, and emits every SLURM `.txt` job array. Also copies `gelman_normalization.R` and `convergence_gate.R` into both bundles and `hpc_mediation.R` + `helper_scripts/` into the mediation bundle, then prints the transfer / dSQ / retrieve commands to copy-paste. |
| `generate_nonsp_diagnostic_jobs.py` | Builds a follow-up dSQ array that re-runs DHARMa diagnostics on already-fitted single-path models via `nonsp_diagnostic_worker.R`, without re-fitting. `--retrieve` pulls the resulting tarball back. |

#### Default configured scope — every model reported in the paper

The CONFIG block at the top of `generate_hpc_jobs.py` is set to produce every
model the manuscript and its supplement report, and nothing else:

| Setting | Value |
|---|---|
| Single-path model types — `BASE_MODELS` (9) | `nice_covariates_spusers`, `nice_covariates`, `nice_covariates_spusers_iqr`, `empirical_covariates_spusers`, `age_control_spusers`, `true_univariate_spusers`, `nice_covariates_beta_spusers`, `nice_covariates_beta_spusers_iqr`, `nice_covariates_spusers_nonan_caps` |
| Mediation model types — `CUSTOM_MED_TYPES` (8) | the same, minus `nice_covariates_beta_spusers_iqr` |
| Single-path predictors — `PART1_PREDICTOR_GROUPS` | `sp_predictors`, `vch_behavior`, `vch_computations`, `sdt_hppd`, `vch_comp_nominal` (15 variables) |
| Single-path DVs — `PART1_DVS` | `caps_vision`, `hppd_binary` |
| Mediators — `_MED_MEDIATORS` (6) | `vchthreshold`, `vchrate`, `vchbeta`, `vchnu`, `vchnunominal`, `vchbetanominal` |
| Mediation pairings — `MEDIATION_ANALYSES` | `spage → M → hppd_binary` and `avgdose → M → caps_vision` |
| Combined job file — `COMBINED_STEM` | `combined_all_analyses.txt` — 231 job lines (135 single-path + 96 mediation) |

**These lists are not free choices — they are the union of what the supplement
scripts read**, and should be re-derived rather than trusted if either side
changes:

| Config setting | Derived from |
|---|---|
| `BASE_MODELS` | `04_visualizations/supplement/sensitivity_analyses.py` → `MODEL_VARIANTS` + `CANONICAL_MODEL_TYPE` / `SECOND_CANONICAL_TYPE` / `THIRD_CANONICAL_TYPE` |
| `CUSTOM_MED_TYPES` | `supplement/sensitivity_analyses_mediation.py` → `SENSITIVITY_MED_TYPES` + `CANONICAL_MED_TYPE` |
| `vch_comp_nominal` in `PART1_PREDICTOR_GROUPS` | `supplement/regression_results_table_nominal_sensitivity.py` |
| `vchnunominal` / `vchbetanominal` in `_MED_MEDIATORS` | `supplement/mediation_results_table_nominal_sensitivity.py` → `MEDIATORS` |

`vch_omega_nominal` is a single-path predictor only: there is no
`vchomeganominal` key in the `MEDIATORS` dict, so it has no mediation
counterpart to generate. The SDT mediators (`dprimeoverall`,
`criterionoverall`, `meanconffas`) are deliberately out of `_MED_MEDIATORS` — no
supplement figure reads them as mediators, though the SDT variables **are**
single-path predictors via `sdt_hppd`.

> **dSQ array indices are 0-based over the job file's lines.** `dSQBatch.py:82`
> uses `for i, l in enumerate(tf): if i == tid`, so task ID *N* runs line *N+1*.
> `combined_all_analyses.txt` opens with two `#` comment lines, so a 233-line file
> yields `--array 2-232` — 231 tasks, comments excluded. An array starting at 2
> is correct, not off-by-one.

That is **40 job lines → 64 fitted models**: 16 mediation jobs plus 24 single-path
lines. Each single-path line fits *both* DVs, because `generate_nonsp_job_array()`
bundles DVs three-per-line (`bundle_by=3`) and this scope has only two — so 24
lines produce 48 models. Keep that in mind when sizing an array or reading `sacct`.

The wider sensitivity sets are **preserved as commented-out lines directly above
each setting**, so restoring any of them is a one-line uncomment: the full-sample
`nice_covariates`, `empirical_covariates_spusers`, `age_control_spusers` and
`true_univariate_spusers` model types; the nominal-parameter HGF predictors
(`vch_comp_nominal`); the SDT mediators (`dprimeoverall`, `criterionoverall`,
`meanconffas`) and nominal HGF mediators; and the `lshs_total` / `baggot_total`
mediation outcomes.

**Nothing mirrors this table.** Both compile scripts read their filters
straight out of this file's CONFIG — see § Compile filters below, which is also
why there is no separate compile config to keep in sync.

### Cluster execution (one model per array task)

| File | Role |
|---|---|
| `nonsp_predictors.R` | Single-path regressions: `dv ~ predictor + covs`. Parses model-type keyword suffixes to apply row-level filters (see § Row-Level Subsetting Keywords). Average predictive comparisons use the **analytic** `E[Y] = g⁻¹(η)` — no DV-step sampling. |
| `hpc_mediation.R` | Two-equation mediation: `spvar → mediator → dv`. Uses `MC_L = 1000` **mediator** draws; the DV step is **analytic** `E[Y] = g⁻¹(η)`. Sources every `.R` in `helper_scripts/`. |
| `gelman_normalization.R` | The Gelman 2SD transform, the complete-case restriction (`drop_incomplete_model_rows()`), formula-term parsing and `resolve_raw_column()`. Sourced by **both** execution scripts. A deliberate port of `modules/data_prep.py :: normalize_analysis_df()` — change one and you must change the other. Hand-written; edit `03_hpc/gelman_normalization.R`, never a copy under a job bundle. |
| `convergence_gate.R` | The Rhat / bulk-ESS / tail-ESS / divergence thresholds and the check that decides whether a model may write its **summary result tables**. Sourced by **both** execution scripts. Thresholds must stay in sync with `06_submission/verification/mcmc_validity_review.py` and the project brms skill. Hand-written; edit `03_hpc/convergence_gate.R`, never a copy under a job bundle. |
| `normalization_vars.R` | **Auto-generated** by `generate_hpc_jobs.py` from `master_config.py`. Carries which variables each normalization rule applies to (`need_non_normalized_vars`, `inplace_normalized_vars`, `normalize_in_place_vars`, `ensure_normalized_copy_vars`, `gelman_skip_vars`). Membership lives here; arithmetic lives in `gelman_normalization.R`. |
| `nonsp_diagnostic_worker.R` | Worker invoked by each task of the diagnostic array built by `generate_nonsp_diagnostic_jobs.py`. |
| `helper_scripts/` | Function definitions sourced by `hpc_mediation.R` at startup — `process_and_save_mediation_model_function.R`, `calculate_mediation_effect_function.R`, `dharma_check.R`, `plot_posterior_w_ci_shading_function.R`. |

### Compiling results (run locally, after jobs finish)

**You do not invoke either compile script by hand.** `generate_hpc_jobs.py`
runs them for you and supplies the filters. Both files must stay in this
directory — the generator resolves them by path, and each is also imported by
the other side of the config coupling described below.

The full sequence, start to finish:

| # | What you do | What happens |
|---|---|---|
| 1 | `python generate_hpc_jobs.py` | Writes `df_foranalysis_master.csv`, the job `.txt` arrays and the R-script bundles; prints the transfer and submit commands. Then prompts to compile — **answer `N`**, nothing has been submitted yet. Declining prints the compile commands for later. |
| 2 | Run the printed `rsync` commands | Job bundles land on the cluster. |
| 3 | Run the printed `dsq` command, then the `sbatch` command **it** prints | dSQ does not submit; it emits an `sbatch` line you run yourself. Walltime is per job: `90:00` single-path, `360:00` mediation. |
| 4 | `squeue -u <user>` | Wait for the array. |
| 5 | `ssh -MNf bouchet`, then re-run `generate_hpc_jobs.py` and answer **`y`** | Pulls and compiles. Re-running the generator is idempotent — it just rewrites the same job files — and is the easiest way to get the filters right without retyping them. Running the two commands from step 1 does exactly the same thing. |

Each script opens one SSH connection, tarballs the relevant CSVs *on the
cluster*, transfers a single archive, and extracts it locally. With 8000+ CSVs
that is far faster than per-file `scp`, and it means there is no separate
retrieve step to run.

| File | Role | Writes to |
|---|---|---|
| `compile_nonsp_results.py` | Single-path regressions: pulls the `summary_dfs` + `diagnostics` CSVs and compiles them into two flat tables, one row per model. Also exports `load_generator_config()`. | `results/sensitivity_analyses_single_paths/existingresults_manuscript.csv`<br>`results/sensitivity_analyses_single_paths/existingresults_manuscript_counterfactual.csv` |
| `compile_mediation_results.py` | Mediation models: pulls the summary CSVs, DHARMa PNGs, posterior-predictive PNGs and `.RData` fits, preserving the remote directory layout. Run as `compile_mediation_results.py subset`. | `results/{dv}/mediation_models/{model_name}/`<br>e.g. `results/hppd_binary/mediation_models/…`, `results/caps_vision/mediation_models/…` |

Those are the paths `04_visualizations/0X_all_figures.py` reads, which is why
the mediation extraction preserves the remote tree rather than flattening it.

> **An SSH ControlMaster must be open before step 5** (`ssh -MNf bouchet`).
> Bouchet requires Duo, so the scripts cannot authenticate on their own.

### Compile filters — one source of truth

**Both compile scripts read their filters out of `generate_hpc_jobs.py`'s CONFIG.**
There is no separate compile config to keep in sync. Adding a model type to
`BASE_MODELS` or `CUSTOM_MED_TYPES` is the only edit required for it to be
generated, submitted **and** pulled.

`compile_nonsp_results.load_generator_config()` does the work. It execs
`generate_hpc_jobs.py` with `HPC_JOBS_CONFIG_ONLY=1`, which:

- skips the dataframe load, so reading the settings needs no data on disk, and
- stops at the `END OF CONFIG DERIVATION` marker — **before the first line that
  writes anything** — so importing the settings can never emit a job file,
  touch the network, or prompt.

It returns the generator's own `MODEL_VARIANTS`, `ALL_PREDS`, `HPPD_CAPS_DVS`,
`VCH_AS_DVS`, `MEDIATION_ANALYSES` and `DVS`, and raises a descriptive error if
any of those is missing — which is what happens if a setting is moved *below*
the marker.

| Script | Subset filters |
|---|---|
| `compile_nonsp_results.py subset` | `ALL_PREDS` (suffixed `_normalized`), `MODEL_VARIANTS`, `HPPD_CAPS_DVS` (+ `VCH_AS_DVS` when `INCLUDE_SP_VCH_JOBS`) |
| `compile_mediation_results.py subset` | `_derive_model_names(MEDIATION_ANALYSES, DVS)` |

> **Anything hoisted above the marker must stay above it.** `DVS`, the predictor
> group resolution (`SP_PREDS` / `VCH_BEH` / `VCH_COMP` / `_resolve_groups` /
> `ALL_PREDS` / `HPPD_CAPS_DVS` / `PART2_PREDS` / `VCH_AS_DVS`) and
> `MEDIATION_ANALYSES` live above the marker for exactly this reason.
> They are pure lookups into `iv_type_dict` and need no data.

### Convergence gate — a model that did not converge writes no results

Both R scripts refuse to write their **summary result tables** when the fitted
model fails any of the four MCMC validity checks. The reason is printed to the
console, so it lands in the SLURM `.out` file and the failure is legible later.

| Check | Threshold |
|---|---|
| Rhat | < 1.01 |
| bulk ESS | ≥ 1000 |
| tail ESS | ≥ 1000 |
| divergent transitions (post-warmup) | 0 |

Rhat and ESS are checked on the **focal coefficients only** — the spvar and
mediator paths the model exists to estimate:

| script | gated parameters |
|---|---|
| `nonsp_predictors.R` | `b_{spvar}`, plus `b_hu_{spvar}` when the hu submodel varies |
| `hpc_mediation.R` | a (`spvar → mediator`), c′ (`spvar → DV`), b (`mediator → DV`), plus the hu counterparts of c′ and b when hu varies |

The divergence check is not parameter-specific and applies to every model,
including covariate-only runs with no focal coefficient.

**Withheld on failure:**
`results/summary_dfs/{dv}.csv`, `results/summary_dfs/{dv}_counterfactual.csv`,
`path_coefficients_summary.csv`, `hu_paths_summary.csv`,
`mc_mediation_summary.csv`, `path_counterfactual_summary.csv`.

**Still written on failure:** the fit `.RData`, the diagnostics CSV, the
convergence / PP-check / DHARMa PNGs, and the files
`process_and_save_mediation_model()` writes. A blocked model leaves a complete
record of *why* it was blocked and no reportable numbers — the absence of a
summary table is the signal.

A named focal parameter that is not in the draws is a **failure**, not a skip:
that means a naming bug or a model that did not fit what it was asked to, and
neither should quietly produce a summary table.

Thresholds live in `convergence_gate.R` and match
`06_submission/verification/mcmc_validity_review.py`
(`RHAT_MAX` / `ESS_MIN` / `DIV_MAX`) and the project brms skill — one set of
numbers in three places; change one, change all three. The two are
complementary: the gate stops an unusable model producing a reportable number at
all, while `mcmc_validity_review.py` still reviews **every** coefficient in
whatever does get written. A model can pass the gate and still be flagged
downstream on a covariate; that is expected.

---

### Generated output (not tracked)

| Path | Contents |
|---|---|
| `output/mediation_analyses/` | Mediation job bundle: `df_foranalysis_master.csv`, the copied `hpc_mediation.R` + `helper_scripts/`, the generated `*_vars.R` files, and one `.txt` per analysis plus `mediation_combined.txt`. Rebuilt from scratch on every `generate_hpc_jobs.py` run. |
| `../data/final/nonsp_predictor_analyses/` | Single-path job bundle and the combined job-array `.txt`. |

### R package requirements

The job lines load the cluster's R module and prepend the personal library
(`$HOME/R/<major>.<minor>`). Most packages — `brms`, `rstan`, `DHARMa`,
`tidyverse`, `bayestestR`, `posterior`, `bayesplot`, `glue`, `zoo` — resolve from
the site-wide `R-bundle-CRAN` module. Two must be installed into the personal
library:

```r
remotes::install_github("Pakillo/DHARMa.helpers")         # DHARMa.helpers
install.packages("distributions3")                        # ONLY for the two *_deprecated.R scripts
```

`distributions3` is **not** attached by the live scripts. Attaching it masks
**`brms::Beta`** with
`distributions3::Beta(alpha, beta)`: the `exclude = c("Gamma", "pdf")` list caught
`stats::Gamma` but not `Beta`, so `get_brms_family("beta")` returned a
`distributions3` object instead of a brms family. Do not re-add the attach line to
the live scripts.

---

## MCMC Validity Checks

`mcmc_validity_review.py` now lives in `06_submission/verification/` — it is a
submission-time audit, not part of producing any manuscript output, and that
directory is gitignored. It scans `existingresults_manuscript.csv` and every
mediation `summary_*.csv` for Rhat > 1.01, ESS < 1000, or divergent transitions
and writes a Markdown report. Its thresholds must stay in sync with
`convergence_gate.R`. See `06_submission/README.md`.

Note the division of labour: the gate in `convergence_gate.R` decides whether a
model writes results **at all**, while the review script re-checks **every**
coefficient in what was written.

## Pipeline Architecture

### Non-SP Predictor Pipeline (`nonsp_predictors.R`)

```
generate_hpc_jobs.py
        │
        │  (1) Build df_foranalysis_master.csv  [single master df, ALL VALUES RAW — no normalization
        │       happens in Python; see § Normalization]
        │  (2) Write categorical_factor_vars.R   [auto-generated from master_config.CATEGORICAL_FACTOR_VARS]
        │      Write normalization_vars.R        [auto-generated from master_config VARIABLE_REGISTRY
        │                                         + BASE_COVARIATE_SETS]
        │      Copy  gelman_normalization.R      [hand-written; the transform itself]
        │      Copy  convergence_gate.R          [hand-written; the Rhat/ESS/divergence gate]
        │  (3) Build covariates_master_dict      [maps model type → covariate list]
        │  (4) Call generate_nonsp_job_array()   [one txt line per (predictor, dv, model_type)]
        │
        └─→  output/*.txt  →  dsq → SLURM → nonsp_predictors.R (per job)
                                                        │
                                                        │  source() the four files above
                                                        │  row-level keyword filters (spusers, iqr, …)
                                                        │  drop rows incomplete on ANY model term
                                                        │  gelman_normalize_df()  ← normalization happens HERE
                                                        │  verify_gelman_normalization()
                                                        │  as.factor() all listed variables
                                                        │
                                                        └─→  HPC_PARENT/{model_type}/{predictor}/{dv}/
                                                               summary_dfs/{dv}.csv
                                                               diagnostics/{dv}_diagnostics.csv
                                                               fits/{dv}_fit.RData
                                                               ...
```

**Key design principles:**
- One master df (`df_foranalysis_master.csv`) is passed to ALL jobs, carrying **raw values**.
  Row-level subsetting (e.g. SP users only) happens **inside R** based on keywords in the model
  type name string — not by building multiple dfs in Python.
- **Normalization is the LAST preparation step, and it happens in R.** Because subsetting is
  R-side, normalizing in Python would centre and scale every subset model on a larger sample
  than the one entering the likelihood. Each script therefore filters rows, drops rows
  incomplete on any model term, and only then Gelman-normalizes — so mean 0 / sd 0.5 holds on
  exactly the rows brms fits. See § Normalization.
- `generate_hpc_jobs.py` prints the exact TRANSFER → HPC → RETRIEVE → COMPILE commands
  needed after each run. When `PROMPT_TO_COMPILE = True` it offers three sequential prompts:
  (1) run `compile_nonsp_results.py`, (2) run `compile_mediation_results.py subset`.
- The `covariates_master_dict` key is the HPC results directory name. Two entries with the
  same key will overwrite each other's HPC results — always use unique keys.
- `_covs_for_variant()` strips R-side subsetting suffixes before resolving the base covariate
  set, so `nice_covariates_spusers_iqr` maps to the same covariates as `nice_covariates`.
- **brms family / distribution strings come entirely from `VARIABLE_REGISTRY` in `master_config.py`** —
  there are no hardcoded family lists in `generate_hpc_jobs.py`. To change how a variable is
  modeled, update its `distribution` field in `master_config.py`. See § "Family Lookup" below.

### Mediation Pipeline (`hpc_mediation.R`)

Two-equation path model: mediator ~ spvar + covs, then dv ~ mediator + spvar + covs.
Outputs are pulled back by `compile_mediation_results.py` and read by the figure
scripts in `04_visualizations/`.

The `MEDIATORS` dict maps short names (e.g. `'vchnu'`) to `(raw_col, norm_col)` 2-tuples.
The `DVS` dict maps short names to column name strings.
Both mediator family and DV family are derived at job-generation time via `_registry_distribution()`
which reads from `VARIABLE_REGISTRY` — no family values are hardcoded.

`generate_hpc_jobs.py` creates `_normalized` columns in the mediation master CSV only for
mediators that are actually referenced in the current `MEDIATION_ANALYSES` config (derived
at runtime — not from a hardcoded list). To add a new mediator, add it to `MEDIATORS` and
ensure its `distribution` is set in `VARIABLE_REGISTRY`.

**`helper_scripts/` convention:** `hpc_mediation.R` sources every `*.R` file in `helper_scripts/`
at startup via `list.files(helpers_dir, pattern="\\.R$", full.names=TRUE)`. R's `list.files()`
defaults to `recursive=FALSE`, so only top-level files in the directory are sourced — **subdirectories
are ignored**. Only function-definition scripts belong at the top level of `helper_scripts/`.
Ad-hoc scripts with top-level executable code (e.g. `load(FIT_PATH)`) must never sit at the
top level of `helper_scripts/`, or they will crash every mediation job at startup.
`generate_hpc_jobs.py` uses `shutil.copytree()`, so anything placed in a subdirectory is
still shipped to the cluster but is never sourced.

`process_and_save_mediation_model_function.R` writes the
`dharma_quantiles_vs_spvar_*` and `dharma_quantiles_vs_mediator_*` columns natively
for every model fit through this pipeline, so no backfill step is needed here.

#### MC Integration Design: `posterior_predict` vs `posterior_epred`

The causal mediation estimand (Imai et al. 2010 g-formula) requires integrating
E[Y(x, M(x'))] over the mediator's conditional distribution:

> ∫ E[Y | X=x, M=m, Z] · p(m | X=x', Z, θ) dm

The key design question is how to evaluate each component. Both the mediator step
and the DV step have a choice: use the conditional mean (posterior_epred-style) or
sample a noisy draw from the full predictive distribution (posterior_predict-style).

**Step-by-step: what happens in each posterior draw (S = 8,000 draws)**

For each posterior draw θ_s and each of the N subjects:

*Step 1 — Compute E[M | X, Z, θ_s]:*
Apply the a-path model's inverse link to the linear predictor to get the conditional
mean µ of the mediator. For example:
- Student-t mediator (e.g. `vch_beta`): µ = b_int + b_predictor × X + b_covs
- ZIBeta mediator (e.g. `vch_bl_yes_0`): µ = plogis(b_int + b_predictor × X + b_covs)

This is the same in all approaches. The question is what you do with µ next.

*Step 2 — Mediator draws (L = 1000 per subject per posterior draw):*

**Plug-in / posterior_epred for mediator (NOT what we do):**
Use µ directly as the mediator value. This is fast but biased whenever the DV is a
nonlinear function of M, because E[f(M)] ≠ f(E[M]) (Jensen's inequality). For example,
with a logistic DV and a wide mediator distribution, plugging in µ systematically
misestimates the average indirect effect.

**Noisy draw / posterior_predict for mediator (what we do — and what BayesGmed does):**
Sample L = 1000 independent mediator values from the full conditional distribution:
- Student-t: M\* ~ Student-t(ν_s, µ, σ_s)
- ZIBeta: flip Bernoulli(zi_s) → if structural zero, M\* = 0; else M\* ~ Beta(µ×φ_s, (1−µ)×φ_s)

These samples include the observation-level scatter (σ_s, φ_s) around µ. Averaging
E[Y] over L draws then approximates the integral over p(M | X, Z, θ_s) without
Jensen's bias. L = 1000 was chosen to be generous; L = 250 is already sufficient
to attenuate MC noise to a negligible level for our dataset sizes (see empirical
comparison below). Increasing L further converges toward the posterior_epred result.

*Step 3 — DV evaluation given M\*:*

This is where our pipeline and BayesGmed diverge.

**posterior_epred for DV — what we do (analytic; the only route as of 8/21/26):**
Compute E[Y | X, M\*, Z, θ_s] deterministically — no observation noise:
- Bernoulli DV: E[Y] = plogis(η_Y)
- Hurdle NB DV: E[Y] = (1 − plogis(η_hu)) × exp(η_mu) / (1 − P₀), where
  P₀ = (φ / (exp(η_mu) + φ))^φ is the NB's own probability of a zero. The count
  component is zero-truncated, so the untruncated NB mean must be renormalised by
  1/(1 − P₀). The `/(1 − P₀)` factor was missing before 8/20/26 (it underestimated
  E[Y] and every contrast derived from it).
  Ref: <https://github.com/paul-buerkner/brms/blob/master/R/posterior_epred.R>

**posterior_predict for DV — BayesGMed-compatible; REMOVED 8/21/26:**
Behind an `add_observation_noise` flag, both scripts could instead draw Y\* from the
full predictive distribution:
- Bernoulli DV: Y\* ~ Bernoulli(plogis(η_Y)) → a 0 or 1
- Normal DV: Y\* ~ Normal(η_Y, σ_Y) → a noisy continuous value
- Hurdle NB: flip the hurdle Bernoulli, then if nonzero draw
  Y\* ~ **ZeroTruncatedNegBin**(exp(η_mu), φ_s) via `distributions3::rztnbinom()`.
  The zero truncation is what makes this a hurdle NB rather than a zero-inflated
  NB, and is what makes the MC route target the same estimand as the analytic
  route above. Before 8/20/26 both scripts used a plain `rnbinom()`, which can
  return 0 after the hurdle has been cleared — inaccurate zero inflation,
  inconsistent with the distribution brms actually fits. Fixed 8/21/26.
  Ref: <https://github.com/paul-buerkner/brms/blob/master/R/posterior_epred.R>

*Step 4 — Average and compute effects:*
Average E[Y] (or Y\*) over the L mediator draws and N subjects → one NIE/NDE value
per posterior draw. The full posterior of NIE/NDE across S draws gives the HDI and p_direction.

**Tradeoffs**

The two DV-step approaches target the **same population-average estimand**. Individual
observation-level Y noise cancels out when averaging over N subjects and L mediator
draws. The practical difference is Monte Carlo efficiency:

- **posterior_epred for DV** computes the conditional mean exactly. No extra noise is
  introduced. The HDI width reflects only parameter uncertainty and mediator-distribution
  integration uncertainty — it does not grow with decreasing L.

- **posterior_predict for DV** adds Bernoulli/NB scatter to each Y draw. This noise
  cancels in expectation (same NIE point estimate) but inflates HDI width. The inflation
  shrinks at rate 1/√L: cranking L → ∞ would cause posterior_predict to converge to the
  same HDI as posterior_epred. At L = 1 (equivalent to BayesGmed's per-draw loop), the
  noise can dominate entirely for small-N or zero-heavy outcomes.

#### DV-step: analytic, and the observation-noise route that was removed

Because the two routes target the same estimand and posterior_predict only ever added
Monte Carlo error, an inflated HDI, and runtime, the `add_observation_noise` switch and
its entire Monte Carlo DV branch were **removed from both scripts on 8/21/26**. The live
scripts evaluate the DV analytically and expose no such flag.

The removed route is preserved verbatim — with both 8/21/26 sampler fixes applied — in
the two deprecated observation-noise scripts, which are not shipped in this
repository. They were reference material:
they are not referenced by `generate_hpc_jobs.py`, are never submitted to the cluster,
and produce no number in the manuscript.

**The mediator is still drawn.** Step 2 above is untouched — L = 1000 mediator draws per
(posterior draw × subject) is the Monte Carlo integration the estimator is built on, and
is what avoids the Jensen bias of plugging in µ. Only the *DV-step* draw was removed.

**Consequence for `distributions3`.** No live script calls `rztnbinom()`, so
`distributions3` is not attached — it masks `brms::Beta` (see § Dependencies).

**BayesGmed implementation**

BayesGmed (Belaid et al. 2023; https://github.com/belayb/BayesGmed) uses
posterior_predict for **both** the mediator and the DV throughout. The implementation
lives in Stan `generated quantities` blocks using `_rng` sampling functions
(`normal_rng`, `bernoulli_logit_rng`, etc.) — see the Stan model files in
`inst/stan/` of that repository (e.g. `NY_NM_single.stan`, `BY_BM_single.stan`).
Their effective L equals N (one loop iteration per subject per MCMC draw, no
separate MC sample size). For the Normal/Binary families they support, observation
noise is modest relative to effect size and the N-subject average sufficiently
attenuates it. BayesGmed does not directly support ZIBeta mediators or hurdle-NB
outcomes, which is why we implement the g-formula estimand ourselves.

**Empirical comparison (ZIBeta mediator → hurdle-NB DV, N = 137)**

| Approach | NIE median | 94% HDI | P(dir) | Notes |
|---|---|---|---|---|
| Plug-in E[M], epred DV (Jensen-biased) | 0.040 | [−0.003, 0.124] | 0.987 | Fast; biased for nonlinear DV |
| Noisy M\*, predict DV, L=1 (BayesGmed-style) | 0.051 | [−0.248, 0.372] | 0.615 | Signal buried in NB/Bernoulli noise |
| Noisy M\*, epred DV, L=250 (prior approach) | 0.046 | [−0.003, 0.124] | 0.986 | Correct integral; low MC noise |
| **Noisy M\*, predict DV, L=1000 (current default)** | **≈0.046** | **≈[−0.003, 0.124]** | **≈0.986** | BayesGMed-compatible; HDI converges to epred result as L→∞ |

The Jensen bias from the full plug-in approach is negligible here (ZIBeta structural
zero-inflation ~4%, smooth Beta over the remainder where the logit link is nearly
linear). The L=1 posterior_predict result is not a different causal answer — it is
the same estimand approximated with a Monte Carlo estimator that is too noisy for
our dataset size. The comparison analysis that produced this table was run during
model development and is not part of the reproduction pipeline.

### Counterfactual Settings

> **Terminology.** Despite the name used throughout this pipeline, these outputs are
> **not causal counterfactuals**. What is computed is an *average predictive comparison*
> in the sense of Gelman & Pardoe (2007), *Average Predictive Comparisons for Models with
> Nonlinearity, Interactions, and Variance Components*, Sociological Methodology
> 37(1):23–51, [doi:10.1111/j.1467-9531.2007.00181.x](https://doi.org/10.1111/j.1467-9531.2007.00181.x):
> for each posterior draw, `E[Y]` is evaluated for every observed subject with the focal
> predictor set to `X0` and then to `X1` — every other covariate held at that subject's own
> observed value — then averaged over subjects and differenced. It is a model-based
> predictive contrast that rests on no identification assumptions, and the study design is
> cross-sectional, so it must not be reported or read as a causal effect. The
> `counterfactual` naming is retained only because it is baked into output filenames,
> CSV column names, and `compile_nonsp_results.py`.

`nonsp_predictors.R` computes `E[Y] = g⁻¹(η)` in closed form, per family, at `CF_X0 = 0.0`
and `CF_X1 = 0.5` (Gelman units, so `X1` is +1 raw SD). There is no sampling and there are
no MC tuning constants — `add_observation_noise` and `CF_MC_L` were removed on 8/21/26
(see § DV-step: analytic, and the observation-noise route that was removed).

For `hurdle_negbinom_huvary` the closed form is

```
E[Y] = (1 − plogis(η_hu)) × exp(η_mu) / (1 − P₀),   P₀ = (φ / (exp(η_mu) + φ))^φ
```

The `1/(1 − P₀)` factor is the zero-truncation renormalisation of the count component;
omitting it (as this code did before 8/20/26) underestimates `E[Y]` and every contrast
derived from it.

**Runtime.** The analytic route is effectively instantaneous — a handful of vectorised
matrix operations on the S × N posterior grid. The previous Monte Carlo route, for
reference, ran ~6 min (`bernoulli`, N = 186) to ~39 min (`hurdle_negbinom_huvary`,
N ≈ 228) per job at `CF_MC_L = 1000`, `S = 36000`. Walltimes in `dsq` submissions were
sized for that route and are now generously over-provisioned.

### Sensitivity Analysis Pipeline (`04_visualizations/supplement/sensitivity_analyses.py`)

A combined four-stage script that mirrors `generate_hpc_jobs.py` for a fixed set of
SP-user-restricted sensitivity model types. Designed to never collide with primary
manuscript results.

#### Model types

| Model type | Base covariate set | Notes |
|---|---|---|
| `empirical_covariates_spusers` | Empirical covariates without `location_europe` | Standard empirical |
| `age_control_spusers` | `age_v2` only | Age-only control |
| `true_univariate_spusers` | None | No covariates |
| `nice_covariates_beta_spusers` | `nice_covariates` + `vch_beta` covariate | Beta-control: tests whether effects survive controlling for decision noise |
| `nice_covariates_beta_spusers_iqr` | `nice_covariates` + `vch_beta` covariate, IQR outlier filter | Beta-control with additional outlier exclusion |

All types apply the `spusers` row filter (`psycheduse_yn == "Yes"`) inside `nonsp_predictors.R`.
The `nice_covariates_beta_*` types are appended to `MODEL_VARIANTS` directly (not derived from
`BASE_MODELS × MODIFIERS`) because they require the `_beta` suffix logic that only applies to
this base model. Their results were generated by the primary manuscript pipeline
(`generate_hpc_jobs.py`) and already exist on HPC — the generate stage is not needed for them.

#### Four stages

```
Stage 1 — generate:
    Reads df_foranalysis_master.csv (built by generate_hpc_jobs.py).
    Writes SLURM .txt job arrays to:
        data/final/sensitivity_analyses/
    Prints TRANSFER + HPC submission commands.

Stage 2 — compile:
    SSHes into Bouchet via ControlMaster (~/.ssh/bouchet.sock).
    Tarballs summary_dfs/ + diagnostics/ CSVs for all four model types.
    Extracts, merges DHARMa diagnostics (all 11 columns incl. q25/q50/q75).
    Writes:
        results/sensitivity_analyses_single_paths/existingresults_sensitivity.csv
        results/sensitivity_analyses_single_paths/existingresults_sensitivity_counterfactual.csv

Stage 3 — heatmap:
    Reads existingresults_sensitivity_counterfactual.csv.
    Cells show ± posterior probability of direction (P(contrast > 0) or P(contrast < 0)).
    Red border = any of 11 DHARMa columns p < 0.05.
    Writes one PNG per DV to:
        results/sensitivity_analysis_heatmaps/sensitivity_analyses/

Stage 4 — pdf:
    Assembles per-DV diagnostic PDFs from:
            results/diagnostics/{dv}_diagnostic_compilation.png
    Writes multi-page PDFs to:
        results/supplement/sensitivity_analyses/
```

#### Usage

```bash
# All enabled stages (controlled by CONFIG toggles):
/usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses.py

# Individual stages:
/usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses.py generate
/usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses.py compile
/usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses.py heatmap
/usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses.py pdf
```

#### Key differences from `generate_hpc_jobs.py`

- Job arrays written to `data/final/sensitivity_analyses/` — never touches
  `data/final/nonsp_predictor_analyses/`
- Compiled CSVs named `existingresults_sensitivity*.csv`, not `existingresults_manuscript*.csv`
- Heatmap reads from counterfactual CSV (response-scale posterior probability), not
  the fixed-effects log-scale coefficient
- DHARMa flagging includes per-quantile p-values (q25/q50/q75) that `compile_nonsp_results.py`
  does not yet include
- PDF output goes to `results/supplement/sensitivity_analyses/`, not
  `04_visualizations/figures/`

#### Adding new model types

`nice_covariates_beta_spusers` and `nice_covariates_beta_spusers_iqr` are appended directly to
`MODEL_VARIANTS` at the bottom of the list comprehension (not derived from `BASE_MODELS × MODIFIERS`)
because the `_beta` suffix needs special handling that only applies to this base covariate set.
Their HPC jobs were generated by the primary `generate_hpc_jobs.py` pipeline, so the generate
stage of `sensitivity_analyses.py` is not needed for them — run compile + heatmap only.

---

### Mediation Sensitivity Pipeline (`04_visualizations/supplement/sensitivity_analyses_mediation.py`)

A combined three-stage script for HPC mediation sensitivity analyses. Produces compound
path heatmaps (A path, B path, C' path, NIE) across covariate specifications, with the
canonical `nice_covariates_spusers` column prepended from the primary results tree.

#### Covariate types

| Type | Covariate set | Notes |
|---|---|---|
| `empirical_covariates_spusers` | Empirical covariates (no `location_europe`), SP users | |
| `nice_covariates` | `age_v2`, `sex_v2`, `raven_total`, `mental_illness2_v2` — full sample | Promoted to SECOND_CANONICAL column |
| `nice_covariates_spusers_iqr` | nice_covariates, SP users, IQR outlier filter | Promoted to THIRD_CANONICAL column |
| `age_control_spusers` | `age_v2` only, SP users | |
| `true_univariate_spusers` | No covariates, SP users | |
| `nice_covariates_beta_spusers` | nice_covariates + `vch_beta`, SP users | Beta-control; `vchbeta` mediator row left blank (circular) |

The `nice_covariates_beta_spusers_iqr` type is omitted — no mediation results exist on HPC for it.

#### Usage

```bash
/usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses_mediation.py generate
/usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses_mediation.py compile
/usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses_mediation.py recompile  # force re-pull from HPC
/usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses_mediation.py heatmap
```

Output heatmaps: `results/supplement/sensitivity_analyses_mediation/{a,b,cprime,nie}_compound_heatmap.png`

---

## Family / Distribution Lookup

**`generate_hpc_jobs.py` derives all brms family strings from `VARIABLE_REGISTRY` in `modules/master_config.py`.**
There are no hardcoded family lists anywhere in the script. Two helpers enforce this:

| Helper | Role |
|---|---|
| `_dv_settings(dv)` | Returns the `settings` string passed to `nonsp_predictors.R` for a given DV column. Reads `VARIABLE_REGISTRY[dv]['distribution']`; has small override sets for pipeline-specific `_hierarchial` suffix variants. |
| `_registry_distribution(col, role)` | Looks up `VARIABLE_REGISTRY[col]['distribution']` for mediation mediators and DVs. Raises a descriptive `ValueError` if the column is missing or has a null distribution — so errors surface at job-gen time, not on HPC. |

**Valid `distribution` strings in `VARIABLE_REGISTRY`** (and their brms equivalents):

| Registry string | brms family | Use case |
|---|---|---|
| `bernoulli` | `bernoulli()` | Binary (0/1) outcomes |
| `logistic` | ⚠️ Not a valid brms family — use `bernoulli` | (legacy; do not use in `master_config.py`) |
| `zero_negbinomial` | `zero_inflated_negbinomial()` | Count with excess zeros |
| `negbinomial` | `negbinomial()` | Overdispersed count, no zero-inflation |
| `hurdle_negbinom_huvary` | `hurdle_negbinomial()` + varying `hu` submodel | Count with excess zeros AND predictor effects on hurdle probability. `bf()` gets two formulae: `dv ~ sp_and_covs` AND `hu ~ sp_and_covs`. Extracts spvar→hu and mediator→hu paths separately in addition to standard mu-submodel paths. |
| `student_t` | `student()` | Heavy-tailed continuous |
| `gamma` | `Gamma(link="log")` | Positive continuous (incl. HGF parameters) |
| `zero_inflated_beta` | `zero_inflated_beta()` | Proportion with excess zeros |
| `ordinal` | `cumulative()` | Ordinal response (Bürkner & Vuorre 2019) |
| `lognormal` | `lognormal()` | Log-normal continuous |

Notes:
- `nonsp_predictors.R` uses `"ordinal"` to trigger `cumulative()`; `hpc_mediation.R` accepts both `"ordinal"` and `"cumulative"` as synonyms.
- `nonsp_predictors.R` accepts both `"bernoulli"` and `"logistic"` as strings (both map to `bernoulli()` family). `master_config.py` only ever emits `"bernoulli"` — the `"logistic"` branch exists only to prevent silent failures from legacy job arrays.
- `hurdle_negbinom_huvary` vs `hurdle_negbinom`: the `_huvary` suffix signals that the `hu` (hurdle/zero-probability) submodel receives its own formula with all predictors and covariates. A plain `hurdle_negbinom` registry string (if ever added) would use a fixed `hu ~ 1` intercept only.

**To change how a variable is modeled:** update its `distribution` field in `master_config.py`. Re-run `generate_hpc_jobs.py` to regenerate job arrays with the updated family.

---

## Result File Schemas

All result CSVs are written by the HPC R scripts during job execution. The per-job output directory structure is:

```
{HPC_PARENT}/{model_type}/{predictor}/{dv}/
    summary_dfs/{dv}.csv                   ← fixed-effects table (nonsp_predictors.R)
    summary_dfs/{dv}_counterfactual.csv    ← counterfactual E[Y] contrast (nonsp_predictors.R)
    diagnostics/{dv}_diagnostics.csv       ← DHARMa diagnostics (nonsp_predictors.R)
    fits/{dv}_fit.RData                    ← brms fit object
{results_dir}/{model_name}/
    *.html                                 ← mediation path diagrams (hpc_mediation.R)
    path_coefficients_summary.csv          ← A/B/C' path coefficients, link scale, 94% HDI (hpc_mediation.R)
    path_counterfactual_summary.csv        ← A/B/C' path effects, response scale, 94% HDI (hpc_mediation.R)
    mc_mediation_summary.csv               ← MC integration causal effects, response scale (hpc_mediation.R)
    hu_paths_summary.csv                   ← hurdle hu-submodel paths (hpc_mediation.R, huvary only)
    convergence_diagnostics.csv            ← divergences / Rhat / ESS (hpc_mediation.R)
```

### `{dv}.csv` — fixed-effects table (`nonsp_predictors.R`)

One row per fixed-effect coefficient. Produced for **every model family**.

| Column | Type | Description |
|---|---|---|
| `var` | string | brms coefficient name with `b_` prefix stripped (underscores also stripped by brms from variable names, e.g. `capsvision_avglifedose`) |
| `Estimate` | float | Posterior mean |
| `Est.Error` | float | Posterior SD |
| `Q2.5` | float | 2.5th percentile of posterior |
| `Q97.5` | float | 97.5th percentile of posterior |
| `Rhat` | float | Gelman–Rubin convergence statistic (target < 1.01) |
| `Bulk_ESS` | float | Bulk effective sample size |
| `Tail_ESS` | float | Tail effective sample size |
| `model` | string | `settings` string (e.g. `"hurdle_negbinom_huvary"`) |
| `covariates` | string | Space-separated covariate string passed to the model |
| `N` | int | Number of observations used in fit |
| `num_divergents` | int | Number of divergent transitions |
| `prob_below_0` | float | Posterior P(coefficient < 0) |
| `prob_above_0` | float | Posterior P(coefficient > 0) |
| `hdi_lower_94` | float | Lower bound of 94% HDI |
| `hdi_upper_94` | float | Upper bound of 94% HDI |

Note: brms strips underscores from response variable names in coefficient name strings. `caps_vision` → `capsvision`, `avg_life_dose` → `avglifedose`. When parsing this file downstream, use `gsub("_", "", col_name)` to match R's naming.

### `{dv}_diagnostics.csv` — DHARMa diagnostics (`nonsp_predictors.R`)

One row per model. Per-quantile p-values are saved alongside the combined
BH-adjusted p, because a single quantile can deviate while the combined test
does not: a quantile line turns red on the residual plot whenever its own
p < 0.05, and only the per-quantile column explains why.

| Column | Description |
|---|---|
| `dharma_uniformity_pval` | `testResiduals()$uniformity` — KS test on residual uniformity |
| `dharma_dispersion_pval` | `testResiduals()$dispersion` |
| `dharma_outlier_pval` | `testResiduals()$outliers` |
| `dharma_heteroscedasticity_pval` | Combined BH-adjusted p from `testQuantiles(sim, predictor = rank(spvar)/n, plot=FALSE)` — spvar is rank-transformed before the call to match what `plotResiduals()` shows internally. When `spvar` is empty, one column per covariate: `dharma_hetero_pval_{cov}`. |
| `dharma_heteroscedasticity_q25_pval` | Individual q25 unadjusted p-value (same `testQuantiles` call). A quantile line turns **red** in the image when this < 0.05, even if the combined p is > 0.05. |
| `dharma_heteroscedasticity_q50_pval` | Individual q50 unadjusted p-value |
| `dharma_heteroscedasticity_q75_pval` | Individual q75 unadjusted p-value |
| `dharma_quantiles_pval` | Combined BH-adjusted p from `plotResiduals(sim)` vs rank-transformed fitted values |
| `dharma_quantiles_q25_pval` | Individual q25 unadjusted p-value from vs-fitted quantile test |
| `dharma_quantiles_q50_pval` | Individual q50 unadjusted p-value |
| `dharma_quantiles_q75_pval` | Individual q75 unadjusted p-value |

When `spvar` is empty (covariate-only analyses), per-quantile columns are emitted for each covariate:
`dharma_hetero_q25_pval_{cov}`, `dharma_hetero_q50_pval_{cov}`, `dharma_hetero_q75_pval_{cov}`.

**`hpc_mediation.R` diagnostic columns** (produced by `process_and_save_mediation_model_function.R`,
one row per response: DV and mediator):

| Column | Description |
|---|---|
| `dharma_ks_pval` | `testUniformity()$p.value` |
| `dharma_outlier_pval` | `testOutliers()$p.value` |
| `dharma_dispersion_pval` | `testDispersion()$p.value` |
| `dharma_zeroinflation_pval` | `testZeroInflation()$p.value` |
| `dharma_heteroscedasticity_pval` | Combined BH-adjusted p from rank-transformed `testQuantiles` vs fitted values |
| `dharma_heteroscedasticity_q25_pval` | Individual q25 unadjusted p-value |
| `dharma_heteroscedasticity_q50_pval` | Individual q50 unadjusted p-value |
| `dharma_heteroscedasticity_q75_pval` | Individual q75 unadjusted p-value |

Heatmap flags: `p < 0.05` → asterisk + red cell. Letter codes: L=linearity, U=uniformity, O=outlier, D=dispersion, H=heteroscedasticity, Q=quantiles.

**Interpretation note:** Always check the per-quantile columns (`_q25/q50/q75`) alongside the combined p.
The combined BH-adjusted p can be non-significant while individual quantiles flag (e.g., q75 p=0.033 → red
line in plot). A model is flagged if **any** of the four heteroscedasticity columns (combined + 3 quantiles)
has p < 0.05.

### `{dv}_counterfactual.csv` — average predictive comparison (`nonsp_predictors.R`)

Computed by `average_predictive_comparison()`, a top-level function in
`nonsp_predictors.R`: set the focal predictor to `X0` and then `X1` for every
observed unit, hold that unit's own covariates at their observed values, take
`posterior_epred()` on both, and average the difference over subjects — one value
per posterior draw.

That is four lines of `posterior_epred()` rather than a hand-built linear
predictor, deliberately: a prefix match on draw-column names **silently scores
monotonic `mo()` covariates as zero for every subject**, because it looks for
`b_mo(<var>)` where brms writes `bsp_mo<var>`. Any covariate set containing
`highest_education_balanced` — `empirical_covariates`
and `empirical_covariates`, **not** `empiric_covariates` — predates the fix
and must be re-run.

> Not a causal counterfactual — see the terminology note under [Counterfactual Settings](#counterfactual-settings).

One row per model. Computes E[Y|X=X1] − E[Y|X=X0] integrating both submodels for hurdle/ZI families.
Produced for all families **except** `ordinal`, `lognormal_hierarchial`, `lognormal_norm`,
`student_t_lognormal_hierarchial`, `student_t_lognormal_norm`. Only written when `spvar` is non-empty.

- **X0 = 0.0** (Gelman-normalized mean, = raw mean of predictor)
- **X1 = 0.5** (Gelman-normalized = +1 raw SD above mean; 0.5 Gelman units × 2×SD = 1 raw SD)
- Note for future analyses: consider X0=−0.5, X1=+0.5 (±1 raw SD symmetric contrast)

| Column | Type | Description |
|---|---|---|
| `dv` | string | Outcome variable column name |
| `settings` | string | brms family / distribution string (e.g. `"hurdle_negbinom_huvary"`) |
| `spvar` | string | Predictor variable (Gelman-normalized column name) |
| `X0` | float | Counterfactual X value for "control" condition (0.0 = normalized mean) |
| `X1` | float | Counterfactual X value for "treated" condition (0.5 = +1 raw SD) |
| `contrast_label` | string | Human-readable label, e.g. `"X1=0.5 vs X0=0.0 (+1 raw SD)"` |
| `scale` | string | Units of the contrast. Actual values: `response_count` (hurdle NB), `probability` (bernoulli/logistic), `response_scale` (student_t), `response_exp_scale` (gamma), `proportion_scale` (beta / zero_inflated_beta). Set from `CF_SCALE_BY_FAMILY` in `nonsp_predictors.R`; an unlisted family is a hard stop, never a blank label. |
| `mean` | float | **Reported point estimate.** Posterior mean of E[Y|X1] − E[Y|X0]. Same summary of the posterior as brms' `Estimate` in the coefficient CSVs |
| `median` | float | Posterior median of the same contrast — retained so the choice of point estimate can be revisited without refitting. Not reported. |
| `hdi_lower_94` | float | Lower bound of 94% HDI of the contrast posterior |
| `hdi_upper_94` | float | Upper bound of 94% HDI of the contrast posterior |
| `prob_above_0` | float | P(contrast > 0) |
| `prob_below_0` | float | P(contrast < 0) |
| `N_obs` | int | Number of observations used |
| `N_draws` | int | Number of posterior draws used |

### `mc_mediation_summary.csv` — MC integration causal mediation (`hpc_mediation.R`)

One row per causal effect. Uses Imai et al. (2010) g-formula (L=1000 mediator draws per observation×draw).
Produced for **all family combinations** whenever `hpc_mediation.R` runs successfully.

- **X0 = 0.0** (Gelman-normalized mean), **X1 = 0.5** (+1 raw SD)
- Effect scale: response scale of the DV (counts for `hurdle_negbinom`; probability for `bernoulli`)

| Column | Type | Description |
|---|---|---|
| `effect` | string | Effect label: `"NIE Indirect"`, `"NDE Direct"`, `"TE Total Effect"`, `"PMed Proportion mediated"` |
| `mean` | float | **Reported point estimate.** Posterior mean of the causal effect. Means are exactly additive (NIE + NDE = TE); medians are not. Applies to the `PMed` row too. |
| `median` | float | Posterior median of the same effect — retained for revisiting the choice without refitting. Not reported. |
| `hdi_low` | float | Lower bound of 94% HDI |
| `hdi_high` | float | Upper bound of 94% HDI |
| `p_above_0` | float | P(effect > 0) |
| `p_below_0` | float | P(effect < 0) |
| `p_direction` | float | max(p_above_0, p_below_0) — probability of direction |
| `n_draws` | int | Number of posterior draws used in MC integration |

PMed (proportion mediated) = NIE / TE. Interpret with caution when TE is near zero or crosses zero; PMed can exceed [0,1] in that case.

### `path_coefficients_summary.csv` — individual path coefficients (`hpc_mediation.R`)

One row per path. Produced for **all family combinations** whenever `hpc_mediation.R` runs successfully.
Reports the A, B (mu), and C' (mu) path coefficients on the **link/latent scale** (NOT response scale).
For `hurdle_negbinom_huvary` models, also includes B (hu) and C' (hu) paths on the logit scale.

Path labels:
- `a`       — A path: spvar → mediator (mediator submodel; log or logit link)
- `b`       — B path (mu): mediator → DV mu submodel (log, logit, or identity link)
- `c_prime` — C' path (mu): spvar → DV mu submodel (direct effect)
- `b_hu`    — B path (hu): mediator → DV hurdle probability (logit; huvary models only)
- `c_hu`    — C' path (hu): spvar → DV hurdle probability (logit; huvary models only)

Naming note: `var_brms_col` uses brms' naming convention (underscores stripped from response names, preserved in predictor names). Reported point estimate = posterior **mean** (consistent with `hu_paths_summary.csv` and `mc_mediation_summary.csv`, and with the `Estimate` column of the single-path coefficient CSVs, which is brms' own posterior mean). The posterior **median** is stored alongside it in every summary CSV so the choice can be revisited without refitting, but is not reported anywhere. The reported summary is named once, in `POINT_ESTIMATE_COL` in `modules/master_config.py`; every Python reader goes through `point_estimate()` from that module rather than naming a column itself.

| Column | Type | Description |
|---|---|---|
| `path` | string | Path label: `"a"`, `"b"`, `"c_prime"`, `"b_hu"`, `"c_hu"` |
| `var_brms_col` | string | Raw brms posterior draw column name (e.g. `"b_vchnunormalized_vch_threshold_normalized"`) |
| `mean` | float | **Reported point estimate.** Posterior mean |
| `median` | float | Posterior median — retained, not reported |
| `hdi_lower_94` | float | Lower bound of 94% HDI |
| `hdi_upper_94` | float | Upper bound of 94% HDI |
| `prob_above_0` | float | P(coefficient > 0) |
| `prob_below_0` | float | P(coefficient < 0) |
| `p_direction` | float | max(prob_above_0, prob_below_0) |
| `n_draws` | int | Number of posterior draws |

**Important:** These are link-scale coefficients, not response-scale effects. For mediation
diagrams, use `path_counterfactual_summary.csv` instead (response-scale A/B/C' paths).

### `convergence_diagnostics.csv` — MCMC convergence record (`hpc_mediation.R`)

One row per fitted model. Produced for **every** `hpc_mediation.R` run, written immediately
after `brm()` returns — before the path, MC, and predictive-comparison blocks — so a failure
in any of those still leaves a convergence record on disk.

| Column | Type | Description |
|---|---|---|
| `model_name` | string | Canonical model name (also the containing directory) |
| `dv` | string | Outcome variable |
| `mediator` | string | Mediator variable (mediator-submodel response) |
| `spvar` | string | SP predictor |
| `N` | int | `nobs(fit)` — rows actually used after keyword filters and NA dropping |
| `iter` | int | Total MCMC iterations per chain |
| `warmup` | int | Warmup iterations per chain |
| `num_divergents` | int | Post-warmup divergent transitions across all chains |
| `max_rhat` | float | Largest Rhat across all parameters |
| `min_bulk_ess` | float | Smallest bulk ESS across all parameters |
| `min_tail_ess` | float | Smallest tail ESS across all parameters |

Divergences are counted from `brms::nuts_params(fit, pars = "divergent__")`, which returns
post-warmup draws only — the same method used by `nonsp_predictors.R`, so the two pipelines
are directly comparable. Rhat and ESS come from `posterior::summarise_draws()`.

Parameters whose Rhat is `NA` (constant or deterministic quantities such as `lprior` when no
explicit prior is set) are excluded from `max_rhat`/`min_*_ess` and named in the job log. If
the Rhat/ESS summary fails outright, those three columns are `NA` and the reason is logged;
`num_divergents` is still recorded.

**Pitfall if you modify this code:** do not drop the explicit `Parameter == "divergent__"`
filter. Other sampler parameters returned by `nuts_params()` (`accept_stat__`, `stepsize__`)
legitimately take the value 1, so summing `Value == 1` over an unfiltered frame silently
inflates the count.

### `path_counterfactual_summary.csv` — counterfactual path effects, response scale (`hpc_mediation.R`)

One row per path. Produced for **all supported family combinations** (skipped gracefully for `cumulative`/ordinal).
Reports individual path effects on the **DV response scale**, computed deterministically by
evaluating E[Y] at fixed contrast points (no additional MC loop — reuses MC block parameters).

Contrast definitions (consistent with `mc_mediation_summary.csv`):
- **X0 = 0.0** (Gelman-normalized mean of spvar), **X1 = 0.5** (+1 raw SD)
- **M_norm = 0.0** (mediator at its raw sample mean), **M_norm = 0.5** (+1 raw SD of mediator = +0.5 Gelman units)
- All effects are averaged over observations (marginalizes over the covariate distribution)

Path definitions:
- `A path`  — E[Mediator | X=X1] − E[Mediator | X=X0]: effect of +1 raw SD spvar on mediator
- `B path`  — E[DV | M_norm=0.5, X=X0] − E[DV | M_norm=0.0, X=X0]: effect of +1 raw SD mediator on DV, spvar held at mean. For `hurdle_negbinom_huvary`, integrates **both** mu and hu submodels through E[Y] = (1−plogis(η_hu)) × exp(η_mu) — no separate B_hu needed.
- `C' path` — E[DV | M_norm=0.0, X=X1] − E[DV | M_norm=0.0, X=X0]: direct effect of spvar on DV, mediator fixed at mean

Scale notes:
- A path: mediator response scale — proportions (0–1) for `beta`/`zero_inflated_beta`; raw positive units for `gamma`; Gelman-normalized units for `student_t`/`gaussian` (1 unit = 2 raw SDs of mediator)
- B and C' paths: DV response scale — directly comparable to NIE/NDE in `mc_mediation_summary.csv`

| Column | Type | Description |
|---|---|---|
| `effect` | string | Path label: `"A path"`, `"B path"`, `"C' path"` |
| `mean` | float | **Reported point estimate.** Posterior mean of the path effect |
| `median` | float | Posterior median of the path effect — retained, not reported |
| `hdi_low` | float | Lower bound of 94% HDI |
| `hdi_high` | float | Upper bound of 94% HDI |
| `p_above_0` | float | P(effect > 0) |
| `p_below_0` | float | P(effect < 0) |
| `p_direction` | float | max(p_above_0, p_below_0) — probability of direction |
| `n_draws` | int | Number of posterior draws |

Note: `hdi_low`/`hdi_high` naming matches `mc_mediation_summary.csv` (both use `effect_summary_hpc()`), NOT `hdi_lower_94`/`hdi_upper_94` used in `path_coefficients_summary.csv`.

**Supported family combinations.** The A / B / C' paths in this file are computed
with `posterior_epred()`, so they carry **no family branching at all** — brms
applies each submodel's own mean function, whatever it is. The constraint that
remains is the NIE/NDE/TE Monte-Carlo block, which draws the mediator and so
needs a family-specific sampler and a family-specific `E[Y]`:

mediator ∈ {`zero_inflated_beta`, `beta`, `student_t`, `gamma`} × DV ∈
{`hurdle_negbinom_huvary`, `hurdle_negbinomial`, `bernoulli`, `student_t`}.
Cumulative/ordinal DVs are skipped with a log message.

There are no `zero_negbinomial`, `negbinomial` or `gaussian` DV branches:
`get_brms_family()` has no `switch` entry for either negbinomial family (the
strings appear only in its error message), and no `VARIABLE_REGISTRY` variable
declares `gaussian`. `student_t` keeps the identity-link branch.

> ⚠️ **51 `VARIABLE_REGISTRY` variables declare a family this script cannot
> construct** — 11 `zero_negbinomial` (`phq9_tot`, `hppd_sx_count`,
> `caps_vision_frequency`, …) and 40 `lognormal` (the `pwPE_*` family). A
> mediation job for any of them `stop()`s at family construction. Open item; see

### `hu_paths_summary.csv` — hurdle hu-submodel paths (`hpc_mediation.R`)

One row per hu-submodel path. Only written when DV family is `"hurdle_negbinom_huvary"`.
Reports the effect of spvar and mediator on the **hurdle probability** P(Y=0), supplementing
the count-submodel (mu) paths reported by `process_and_save_mediation_model()`.

| Column | Type | Description |
|---|---|---|
| `path` | string | Path label: `"c_hu"` (spvar→hu) or `"b_hu"` (mediator→hu) |
| `var_brms_col` | string | Raw brms posterior draw column name (e.g. `"b_hu_capsvision_avglifedose"`) |
| `mean` | float | **Reported point estimate.** Posterior mean |
| `median` | float | Posterior median — retained, not reported |
| `hdi_lower_94` | float | Lower bound of 94% HDI |
| `hdi_upper_94` | float | Upper bound of 94% HDI |
| `prob_above_0` | float | P(coefficient > 0) |
| `prob_below_0` | float | P(coefficient < 0) |
| `p_direction` | float | max(prob_above_0, prob_below_0) |
| `n_draws` | int | Number of posterior draws |

The hu link is **logit**: positive coefficients → higher P(Y=0) (hurdle more likely); negative → more nonzero responses.
brms naming for hu paths: `b_hu_{response_brms}_{predictor_brms}` (underscores stripped from both).

### Notes for downstream scripts (`04_visualizations`, `05_results_narrative`)

1. **brms strips underscores** from all variable names in coefficient column names. Always apply `gsub("_", "", name)` when constructing coefficient name strings. E.g. `caps_vision` → `capsvision`, `vch_bl_yes_0_normalized` → `vchblyes0normalized`.
2. **`{dv}.csv`** contains mu-submodel coefficients only. For `hurdle_negbinom_huvary` models, hu-submodel paths are in `hu_paths_summary.csv` in the mediation results directory.
3. **`mc_mediation_summary.csv`**, **`path_coefficients_summary.csv`**, and **`hu_paths_summary.csv`** live at `{results_dir}/{model_name}/`, not inside `summary_dfs/` — they are mediation-level outputs, not per-model-type outputs.
4. **HDI naming conventions** (project-canonical):
   - brms coefficient tables (`{dv}.csv`, `hu_paths_summary.csv`, `path_coefficients_summary.csv`): `hdi_lower_94` / `hdi_upper_94`
   - mediation effect tables (`mc_mediation_summary.csv`): `hdi_low` / `hdi_high`
   - Forest plot standardized CI columns: `standardized_lci` / `standardized_uci`
   - **How they are computed:** `compute_hdi_94()` — a thin wrapper over
     `bayestestR::hdi(ci = 0.94)` returning `c(lower = ..., upper = ...)`. It is defined in
     `nonsp_predictors.R` and in each `helper_scripts/` file that uses it. It **errors** on any
     non-finite draw (NA/NaN/Inf) rather than dropping it silently, and errors when the sample is
     too short to support a 94% interval. Every one of these errors carries the condition class
     `HDI_ERROR_CLASS` (`"hdi_integrity_error"`). Model steps wrapped in `tryCatch()` —
     the counterfactual block in `nonsp_predictors.R` and the coefficient-table block in
     `process_and_save_mediation_model_function.R` — re-raise on that class instead of logging
     and continuing, so a bad posterior **fails the SLURM task** (non-zero exit, findable via
     `dsqa -s FAILED`) rather than exiting 0 with a missing CSV. Unrelated errors in those same
     blocks are still caught and logged as before. Requires the `bayestestR` package (already a
     transitive dependency of `performance`, and loaded explicitly in `nonsp_predictors.R` and
     `hpc_mediation.R`). It replaced a hand-rolled Kruschke sliding-window implementation;
     the two were verified bit-identical (max abs diff 0) across 79 fixed-effect parameters from
     18 saved fits spanning the bernoulli, hurdle_negbinomial, and zero_inflated_negbinomial
     families, so this change does not move any published number.
5. **`prob_above_0` / `prob_below_0`** are raw posterior proportions (not one-tailed p-values). `p_direction = max(prob_above_0, prob_below_0)` is the probability of direction (pd) used in Bayesian reporting.
6. **Path scale distinction:** `path_coefficients_summary.csv` = link-scale coefficients (log/logit); `path_counterfactual_summary.csv` = response-scale effects (same scale as `mc_mediation_summary.csv`). For mediation diagrams displaying individual path strength, use `path_counterfactual_summary.csv`. For the B path in `hurdle_negbinom_huvary` models, `path_counterfactual_summary.csv` is the **only** file with a single integrated B path that correctly combines mu and hu submodel effects.

---

## Normalization

**Gelman normalization:** `(x - mean(x)) / (2 * sd(x))` — divide by **2×SD**, not 1×SD.

### Where it happens

**The HPC master CSV ships RAW. Normalization happens on the cluster, in R, as the
last preparation step before `brm()`.**

Why: the row-level subsetting keywords are applied inside the R scripts, after the
CSV is read. Normalizing in Python meant every subsetted model was centred and
scaled on a *larger* sample than the one entering the likelihood — mean was not 0
and sd was not 0.5 on the fitted rows, so `X1 = 0.5` did not mean "+1 raw SD" of
the analysis sample. Each R script now runs:

```
read raw CSV
  → row-level keyword filters          (spusers, iqr, …)
  → drop rows incomplete on ANY model term   drop_incomplete_model_rows()
  → gelman_normalize_df()                    ← normalization
  → verify_gelman_normalization()            ← stops on any violated invariant
  → as.factor() / as.ordered()
  → brm()
```

The complete-case step uses every variable the model references — DV, focal
predictor, mediator and `mediator_in_dv` (mediation), and all covariates.
brms would drop those rows anyway,
but only *after* normalization; doing it first makes the normalization sample
and the estimation sample the same set.

| Concern | Lives in | Kind |
|---|---|---|
| Which variables each rule applies to | `master_config.py` → exported to `normalization_vars.R` | auto-generated |
| The arithmetic | `03_hpc/gelman_normalization.R` | hand-written |
| The Python original it is ported from | `modules/data_prep.py :: normalize_analysis_df()` | hand-written |

`gelman_normalization.R` is a deliberate line-for-line port: same four rules, same
order, same "skip columns already at mean 0 / sd 0.5" guard, same NaN handling.
Running it on the full unsubset dataframe reproduces the Python output — 3,939 of
3,989 numeric columns bit-identical, the rest differing only by R-vs-NumPy
summation in the last ulp. **If you change a rule in one, change it in the other.**

| Script | Normalizes |
|---|---|
| `03_hpc/nonsp_predictors.R`, `03_hpc/hpc_mediation.R` | their own analysis sample, at run time |
| `03_hpc/generate_hpc_jobs.py` | **nothing** — it exports rule membership only |
| `04_visualizations/supplement/sensitivity_analyses.py` | sensitivity master df, in Python |
| `04_visualizations/supplement/sensitivity_analyses_mediation.py` | mediation sensitivity df, in Python |

The two sensitivity scripts still call `normalize_analysis_df()`. That is
harmless — Gelman standardization is invariant under a prior positive affine
transform, so R re-normalizing on the analysis sample lands on the same values
either way — but it means the master CSVs they build are pre-normalized while the
one `generate_hpc_jobs.py` builds is not.

Do not add normalization anywhere else — extend `normalize_analysis_df()` and its
R port together.

### Which convention applies to a variable

**Declared per variable in `VARIABLE_REGISTRY` (`modules/master_config.py`).**
`NEED_NON_NORMALIZED` and `INPLACE_NORMALIZED` are *derived* from those flags, the
same way `CATEGORICAL_FACTOR_VARS` is — they are not hand-maintained lists. To change
how a variable is normalized, edit its registry entry.

| Registry flag | Raw column | `{col}_normalized` | Applies when |
|---|---|---|---|
| `need_non_normalized=True` | preserved on its original scale | separate normalized column | the brms family constrains the scale: `gamma` (positive), `beta` / `zero_inflated_beta` (0–1), `hurdle_negbinomial` / `negbinomial` (non-negative integer) |
| `inplace_normalized=True` | Gelman-normalized | alias of the raw column | the family is an unbounded real (`student_t`, `gaussian`) |
| neither | Gelman-normalized in place | created for predictors | ordinary continuous covariate or predictor |

A registry key whose *name* already ends in `_normalized` names a column that is
normalized by construction — leave both flags `False` on those entries.

Current membership:

| List | Variables |
|---|---|
| `NEED_NON_NORMALIZED` | `caps_vision`, `mean_conf_fas`, `vch_bl_yes_0`, `vch_bl_yes_75`, `vch_hit_rate`, `vch_nu`, `vch_nu_avg`, `vch_nu_nominal` |
| `INPLACE_NORMALIZED` | `criterion_overall`, `d_prime_overall`, `vch_beta`, `vch_beta_avg`, `vch_beta_nominal`, `vch_omega`, `vch_omega_avg`, `vch_omega_nominal`, `vch_threshold` |

**Categorical variables are excluded from normalization.** Variables in
`CATEGORICAL_FACTOR_VARS` are written to the master CSV with raw values intact; R calls
`as.factor()` on them via the auto-generated `categorical_factor_vars.R`. Variables in
`MONOTONIC_COVARIATES` are likewise excluded — `mo()` requires raw integers.

### Empty string → NA for character factor covariates

`read.csv()`'s default `na.strings` is `"NA"` only, so a blank field in a
**character** column arrives as `""`, not `NA`. For a factor covariate that is
silently wrong three times over:

1. the missing group survives `drop_incomplete_model_rows()`, which tests for NA;
2. `as.factor()` turns `""` into a real level;
3. `""` sorts first, so it becomes the brms **reference level** — every other
   coefficient would be reported relative to the people who did not answer.

Numeric columns are unaffected (`read.csv` already yields NA for a blank numeric
field), which is why this went unnoticed until the first character-valued factor
covariate was added. Both `nonsp_predictors.R` and `hpc_mediation.R` now convert
`""` → NA for every **character-typed** column in `binary_factor_vars`,
immediately after sourcing `categorical_factor_vars.R` and well before any
filtering, complete-case deletion or normalization. The conversion is logged
(`[blank->NA] <col>: <n> empty string(s) converted to NA`).

The only character-typed entry in `binary_factor_vars` is
`monitor_check_operationalized_final`; every other factor var is numeric with
zero blanks. The conversion guards the whole class of bug rather than that one
column.

### Verification

`normalize_analysis_df()` calls `verify_normalization()` before returning, which raises
`ValueError` if any registry-declared invariant is violated:

- a `NEED_NON_NORMALIZED` variable whose `{col}_normalized` column is missing or is not at mean 0 / sd 0.5
- a `NEED_NON_NORMALIZED` variable whose **raw** column *is* at mean 0 / sd 0.5 — its family cannot accept normalized values, and the raw values are no longer recoverable from the dataframe
- an `INPLACE_NORMALIZED` variable whose raw column or `_normalized` alias is missing or not normalized

`verify_gelman_normalization()` in `gelman_normalization.R` asserts the same invariants
R-side, and is called by both execution scripts immediately after normalizing — so a
violation stops the job rather than fitting a model on a mis-scaled dataframe.
Normalization is applied only to columns not already at mean 0 / sd 0.5 (that pair is the
fixed point of the transform), so calling either function twice is safe.

**`_normalized` columns are not in the master CSV.** They are created by
`gelman_normalize_df()` on the cluster. `generate_hpc_jobs.py` validates predictor and
mediator columns against `_WILL_EXIST_IN_R` — its prediction of what the R normalizer will
produce — rather than against `df_foranalysis.columns`, which would reject every one.

**CRITICAL — always use `_normalized` as pred_col in job gen scripts, even for in-place predictors:**

`compile_nonsp_results.py` appends `_normalized` to every raw predictor name from
`compile_nonsp_results.py` when building HPC `find` commands (line: `[f"{p}_normalized" for p in cfg["ALL_PREDS"]]`).
This means the HPC results directory must be named `vch_beta_normalized/`, not `vch_beta/`.

Job gen scripts must therefore always pass the `_normalized` column as pred_col — e.g.,
`pred_col = 'vch_beta_normalized'`, NOT `pred_col = 'vch_beta'` — even though both columns
contain identical values (vch_beta is normalized in-place). If you pass the un-suffixed name,
results land in `vch_beta/` and compile cannot find them. This is a silent failure with no error message.

---

## Row-Level Subsetting Keywords

Both R execution scripts (`nonsp_predictors.R` and `hpc_mediation.R`) parse
keywords embedded in the model type / model name string to apply row-level filters
to the dataframe.

**Filters are applied in this order:**

| Keyword | Filter applied |
|---|---|
| `spusers` | Keep `psycheduse_yn == "Yes"` |
| `nopsychosis` | Keep `psych_spectrum_v2 < 1` |
| `iqr` | Drop rows outside `[Q1 − 1.5×IQR, Q3 + 1.5×IQR]`. In `nonsp_predictors.R` the fence applies to spvar. In `hpc_mediation.R` it applies to **spvar and the mediator**: a row is kept only if both values sit inside their own fence. Both fences are computed on the frame **as it stands when the filter runs** — i.e. after `spusers` / `nopsychosis`, and before either fence is applied, so the mediator's quartiles describe the frame entering the step rather than whatever the spvar fence left behind. The two fences are independent; only their row-dropping composes, so the retained set does not depend on the order they are applied in. |
| `nocurrenthppd` | Keep `persist_vis_current == 0` |

| `nonan_caps` | Keep rows with non-missing `caps_bl_1` |

The keyword set is defined by `R_SIDE_SUFFIXES` in `modules/master_config.py`;
this table and the `grepl()` blocks in both R scripts must stay in step with it.
Every keyword is reachable from a configured model type — an unused keyword is an
unused way to subset the sample, which this pipeline does not carry.

These columns must be present in `df_foranalysis_master.csv` for the VCH-QC keywords
to work. Regenerate the master CSV (re-run the generation notebook) any time
`modules/data_prep.py` changes the list membership.

---

## Updating `r_side_suffixes` / `R_SIDE_SUFFIXES`

**`R_SIDE_SUFFIXES` in `modules/master_config.py` is the single source of truth.**
`generate_hpc_jobs.py` imports `R_SIDE_SUFFIXES` directly from there — do **not**
hardcode keywords anywhere else.

When adding a new subsetting keyword, update **all** of the following:

1. `modules/master_config.py` — add the new suffix to `R_SIDE_SUFFIXES` (with a descriptive comment)
2. `nonsp_predictors.R` — add a `grepl()` block after the `nocurrenthppd` block
3. `hpc_mediation.R` — same (uses `model_name`, `message()`)

A new filter must not depend on the normalized scale — filters run *before*
normalization, so no `_normalized` column exists yet. Use the raw column
(`resolve_raw_column()` maps a `{x}_normalized` argument to `{x}`), and prefer a
scale-free criterion such as a quantile fence. This is exactly why `nooutlier`
was retired. **Removing a suffix** is the inverse: drop it from `R_SIDE_SUFFIXES`
(so `_covs_for_variant()` can no longer resolve such a variant), add it to
`_RETIRED_KEYWORDS` in `generate_hpc_jobs.py`, and `stop()` on it in both R
scripts to catch job files already staged on the cluster.

Edit `hpc_mediation.R` in this directory, not the copy under
`output/mediation_analyses/` — that one is overwritten on every
`generate_hpc_jobs.py` run.

Choose keyword names carefully: no new keyword should be a substring of any existing keyword
(or vice versa). This eliminates the need for ordering constraints or `else if` chains in R.

---

## ⚠️ CATEGORICAL VARIABLES — MANDATORY POLICY ⚠️

> **THIS IS THE MOST COMMON SOURCE OF SILENTLY WRONG BRMS RESULTS IN THIS PIPELINE.**
> **READ THIS SECTION BEFORE ADDING ANY COVARIATE OR PREDICTOR.**

### The problem

brms is an R package. R does not know whether a numeric column is a continuous measure or a
discrete factor unless you tell it explicitly with `as.factor()`. If a multi-level categorical
variable (e.g. `mental_illness2_v2` with levels 0/1/2/3 for none/anxiety/depression/psychosis)
is passed to brms as a plain numeric column, brms fits a **single linear slope across the integer
codes** — treating the categories as if they form a meaningful continuum. **The model will run
without any error or warning. The coefficient estimates will be wrong. You will not know.**

### The rule

**Every variable with discrete, unordered levels must be declared in exactly one place:**

> **`modules/master_config.py` → `CATEGORICAL_FACTOR_VARS`**

`generate_hpc_jobs.py` reads this list and automatically:
1. Excludes those variables from Gelman normalization
2. Writes `categorical_factor_vars.R` (defining `binary_factor_vars`) to both output directories
3. Enforces a low-cardinality safety check (hard stop if any model covariate has ≤3 levels and is undeclared)

Both `nonsp_predictors.R` and `hpc_mediation.R` `source()` the auto-generated file and call
`as.factor()` on every variable in the list before fitting brms. **Do not add factor variables
to the R scripts directly** — they will be overwritten the next time `generate_hpc_jobs.py` runs,
and the Python normalization exclusion will be out of sync.

### Current categorical variables

| Variable | Levels | Notes |
|---|---|---|
| `mental_illness2_v2` | 0 = none, 1+ = various Dx | Raw from Python (not normalized); `as.factor()` in R |
| `sex_v2` | 1/2 (binary) | Raw; `as.factor()` in R |
| `amph_lifetime`, `race_bipoc`, `race_asian`, `inhalants_lifetime`, `coke_lifetime` | 0/1 (binary) | Raw; `as.factor()` in R |
| `psych_spectrum_v2` | 0/1 (binary) | Raw; `as.factor()` in R |

See `modules/master_config.py :: CATEGORICAL_FACTOR_VARS` for the authoritative list.

### Ordinal covariates are a separate policy — `MONOTONIC_COVARIATES`

An **ordered** categorical variable does not belong in `CATEGORICAL_FACTOR_VARS`.
`as.factor()` would spend K−1 free parameters and throw away the ordering; a plain
numeric column would impose equal spacing between adjacent levels. The middle
course is a brms monotonic effect: one slope plus a K−1 simplex that lets the data
decide the spacing while keeping the order.

> **Single source of truth: `MONOTONIC_COVARIATES` in `modules/master_config.py`.**

`generate_hpc_jobs.py` reads it and automatically:

1. wraps the term as `mo(<var>)` when building the formula RHS (`get_covs()`);
2. excludes it from Gelman normalization — `mo()` requires raw integer codes;
3. writes `monotonic_covariates_vars.R` (defining `monotonic_covariate_vars`) into
   **both** job bundles.

Both `nonsp_predictors.R` and `hpc_mediation.R` source that file and call
`as.ordered(as.integer(x))` on every listed variable **after** row filtering and
complete-case deletion, immediately before `brm()`.

| Variable | Levels | Notes |
|---|---|---|
| `highest_education_balanced` | 1–6 balanced education ordinal | Raw; `as.ordered()` in R, `mo()` in the formula |

**Convert after subsetting, never before.** The conversion has to see the final
analysis sample, because which levels are *present* determines the number of
simplex parameters. brms accepts a bare numeric in `mo()` and fits it without
warning, using `x - min(x)` — which keeps an empty interior level, where
`as.ordered(as.integer(x))` applied after subsetting drops it and renumbers. A
row filter that emptied one level would therefore make the two execution scripts
fit different monotonic specifications for the same covariate, with no error and
no difference in the job line.

### Safety net

`generate_hpc_jobs.py` runs a **low-cardinality check** before saving `df_foranalysis_master.csv`:
any numeric column with ≤3 observed levels that is NOT in `CATEGORICAL_COVARIATES` triggers a hard
stop with a descriptive error. **Do not bypass this check** by just pressing `y` — investigate
whether the variable is genuinely continuous or needs to be added to the categorical lists.

---

## HPC Job Submission Pattern

```python
combined_stem = 'my_combined_file'
combined_path = combine_job_files([...], combined_stem)
print(f'\nTRANSFER:\n  scp -r {os.path.abspath(OUTPUT_BASE)} msg74@transfer-bouchet.ycrc.yale.edu:{HPC_PARENT}')
print(f'\nTo run on HPC:\n  ssh msg74@bouchet.ycrc.yale.edu')
print(f'  cd {HPC_BASE}; module load dSQ; dsq --job-file {combined_stem}.txt --mem-per-cpu 4g -t 90:00 --mail-type ALL')
```

HPC remote base: `/nfs/roberts/scratch/pi_arp29/msg74/aim1_baseline_final/nonsp_predictor_analyses`

Check job status: `squeue -u msg74` / `sacct -u msg74`
