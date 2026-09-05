# 02 — HGF Modeling

Hierarchical Gaussian Filter (HGF) modeling of the visual conditioned-hallucination
(VCH) task for *A Computational and Behavioral Exploration of Serotonergic
Psychedelics' Persisting Perceptual Effects*.

This directory contains everything needed to go from raw trial-level task data to
the three HGF parameters used throughout the manuscript — **ν** (prior weighting),
**β** (decision precision), and **ω₂** (belief-updating rate) — plus the full
Bayesian-workflow validation reported in Supplementary Figures S2 and S3.

Models are fit in Julia (Turing.jl / ActionModels.jl / HierarchicalGaussianFiltering.jl);
staging, import, and figures are Python. Every script below opens with a **ROLE IN
THE MANUSCRIPT** block stating which result or panel it produces and its position
in the pipeline — that block is the ground truth if anything here goes stale.

---

## Order of operations

Pick the branch below that matches what you're actually trying to do. Most readers
only need the first one.

### "I just want the manuscript figures"

Run nothing in this directory. Supplementary Fig. **S1** (`beta_sigmoid_creator.py`,
in `04_visualizations/supplement/`) is pure math and reads nothing. **S2** and **S3**
are multi-panel assemblies (`hgf_ppc_oos_assembly.py`, `hgf_6panel_assembly.py`, also
in `04_visualizations/supplement/`) that read panel images and CSVs already vendored
in `julia_outputs/` and `model_comparison/bms/bms_summary.csv` — a fresh clone has
everything they need. Just run the four `04_visualizations/supplement/hgf_*.py`
scripts as listed in Step 4 of the top-level README.

### "I'm re-running the HGF fit itself" (new data, or to re-derive ν/β/ω₂)

This is the one linear dependency chain in the directory. Each stage requires the
previous one to have finished:

1. **`hgf_pipeline.py`** (local) — stages per-subject CSVs, builds the SLURM job
   array, prints the transfer/submit commands.
2. **`model_fitting_singleagent_forarray.jl`** (HPC, one task per participant) —
   fits the HGF, writes chains + posterior medians + state trajectories.
3. **`import_hgf_results_unified.py`** (local) — pulls the HPC output together,
   applies the R-hat gate, derives pwPE, writes the wide and long dataframes that
   every downstream statistical model and figure reads.

See **Script reference** below for exact commands and outputs at each stage.

### "I want to regenerate S2/S3 from scratch" instead of trusting the vendored files

Four validation analyses feed S2/S3. Two of them need nothing but the Stage 1
staging CSVs; the other two need Stage 1–3 to have already produced fitted
parameters, because they use each participant's own posterior medians as
generative truth:

| Track | Needs Stage 1–3 fit first? | Feeds |
|---|---|---|
| **A — Prior-based recovery** (`prior_recovery_vch_mcmc.*`, `prior_recovery_aic_bic.*`) | No — simulates from the prior, not from real per-subject data | S3a, S3b |
| **B — Model selection / BMS** (`bms_vch.*`) | No — but optionally warm-starts its MAP fits from Stage 3 medians if present | S3c |
| **C — Posterior-based recovery** (`param_recovery_vch.*`) | **Yes** — needs `vch_nu`/`vch_beta`/`vch_omega` from Stage 3 as the generative truth | S3d, S3e |
| **D — Posterior predictive checks** (`ppc_classic_vch.*`) | **Yes** — same reason as C | S3f, S2b |

Within each track, run the `.jl` file (HPC or local, per script) before the `.py`
file that compiles and plots its output — see **Script reference**. After all four
tracks have produced their panel images/CSVs, run the four assembly scripts in
`04_visualizations/supplement/` to rebuild S2 and S3 themselves.

### "I want to re-derive the empirical stimulus-intensity mapping"

Rarely needed — only if you're questioning the four hardcoded ν/β/ω-independent
constants that convert nominal contrast (0/25/50/75 %) into the detection
probabilities the HGF actually sees. Run
`out_of_set_data/empirical_condition_analysis.py`. It requires the out-of-set
cohort's raw trial data, which is **not distributed** with this repository (see
*Data availability* below) — the derived CSVs and figures it already produced are
included so the derivation can be checked without re-running it.

---

## Script reference

> **HPC tip, applies to every `.jl` script below run on Bouchet.** Bouchet has
> heterogeneous CPU types and Julia's precompilation cache is CPU-target-specific.
> Export `JULIA_CPU_TARGET=generic` when precompiling *and* when running, or tasks
> landing on a non-matching node will spend their entire wall time recompiling and
> produce no output.

### Stage 0 — stimulus intensity derivation

<details><summary><code>out_of_set_data/empirical_condition_analysis.py</code></summary>

- **Prerequisites:** none to run downstream scripts against its output (already
  computed). To re-run it yourself: the out-of-set cohort's raw trial CSV, not
  distributed here.
- **What it does:** verifies which computation reproduces the pipeline's hardcoded
  `empirical_condition` values, then plots per-condition detection probability for
  the full QC-passing sample and the non-hallucinator subset (n=29) used as the
  production mapping.
- **Outputs:** `empirical_condition_verification.csv`,
  `empirical_condition_full_vs_nonhall.csv`,
  `figures/empirical_condition_by_condition.png`,
  `figures/empirical_condition_by_condition_spusers.png`.
- **Used by:** `figures/empirical_condition_by_condition_spusers.png` is read
  directly by `04_visualizations/supplement/hgf_ppc_oos_assembly.py` as panel **a**
  of Supplementary Fig. S2. The four derived scalars themselves are hand-copied
  into `VCH_CONDITION_TO_EMPIRICAL`/`VCH_STAGING_MAPPINGS` (`hgf_pipeline.py`),
  `VCH_CORRECTED_MAPPING` (`model_fitting_singleagent_forarray.jl`),
  `VCH_STIM_MAPPINGS` (`ppc_classic_vch.jl`), and `EMPIRIC_INTENSITIES`
  (`prior_recovery_vch_mcmc.jl`) — see *Stimulus intensity* below.

</details>

### Stage 1–3 — main HGF fit

<details><summary><code>hgf_pipeline.py</code> — Stage 1, local</summary>

- **Prerequisites:** none (entry point). Julia environment instantiated
  (`julia --project=. -e "using Pkg; Pkg.instantiate()"`).
- **What it does:** parses each participant's compressed JSON task blob out of the
  wide df, writes one staging CSV per participant, builds the dSQ job-array `.txt`
  file, and prints the exact transfer/submit/pull commands.
- **Outputs:** `data_n_cmnds/<project>/<timepoint>/vch_data/<record_id>.csv`
  (gitignored — per-subject data), a dSQ job file.
- **Used by:** `model_fitting_singleagent_forarray.jl` reads the staged CSVs on
  HPC. `load_public_wide_df()`, defined here, is the shared dataframe-loading entry
  point every other script in this directory imports.

</details>

<details><summary><code>model_fitting_singleagent_forarray.jl</code> — Stage 2, HPC</summary>

- **Prerequisites:** Stage 1 has staged this participant's CSV and transferred it
  (plus this script and the Julia environment) to the cluster.
- **What it does:** fits one participant's HGF via MCMC (NUTS, MAP init, 4 chains ×
  1000 iterations, `MCMCSerial`), then replays the fit to extract trial-by-trial
  belief trajectories. One SLURM array task per participant (two bundled per task).
  Recomputes `empirical_condition` from `condition` at runtime rather than trusting
  the staging CSV, so a stale staging file can never silently propagate an old
  stimulus mapping into a new fit.
- **Outputs:** `results/<model_type>/<project>/<timepoint>/vch/` — `<id>.jls` (full
  chains), `<id>.csv` (flat posterior samples), `<id>_medians.csv`,
  `<id>_rhats.csv`, `<id>state_trajectories.csv`. See *Output files* below.
- **Used by:** `import_hgf_results_unified.py`, after these results are pulled back
  locally.

</details>

<details><summary><code>import_hgf_results_unified.py</code> — Stage 3, local</summary>

- **Prerequisites:** Stage 2 output pulled back from HPC into the local
  `results/` tree.
- **What it does:** collects posterior medians and state trajectories per
  participant, enforces the R-hat gate (0.9–1.1; failing parameters → `NaN`),
  derives pwPE and time-aligns it (the raw HGF output leads the response by one
  trial), and merges everything into the wide and long analysis dataframes. Backs
  up the input wide df to `data/final/backups/` before overwriting it.
- **Outputs:** `data/final/df_withstates_public_<date>.csv` (wide, one row per
  participant — the primary input to every statistical model in this manuscript),
  `data/final/vch_master_withstates_<date>.csv` (long, one row per trial).
- **Used by:** every script in `03_hpc/`, `04_visualizations/`, and
  `05_results_narrative/` that reads `vch_nu`, `vch_beta`, `vch_omega`, or any
  `vch_*` state-trajectory column. Also the required upstream for validation
  tracks C and D below (they need this script's fitted medians as generative
  truth).

</details>

<details><summary><code>helper_functions/create_agent.jl</code> — shared, not run directly</summary>

- **What it does:** defines the generative model itself — the HGF submodel plus
  the perceptual/response model (the "Learn → Perceive → Respond" steps described
  under *The task and the model*). `include`d by every `.jl` script in this
  directory: Stage 2, both `bms_vch.jl` and `ppc_classic_vch.jl`. (The two
  `prior_recovery_*.jl` and `param_recovery_vch.jl` scripts reimplement the same
  forward pass inline rather than including this file, so that they don't have to
  fight its I/O assumptions when simulating instead of fitting — check each
  script's own header if the two implementations ever need to be reconciled.)

</details>

### Validation track A — prior-based recovery and model identifiability (S3a, S3b)

<details><summary><code>prior_recovery_vch_mcmc.jl</code> — HPC</summary>

- **Prerequisites:** none from Stage 1–3 — parameters are sampled from each
  candidate model's own prior, not from real participant fits.
- **What it does:** for each of 500 simulation indices × 4 true models, draws
  generative parameters from that model's prior, simulates a 360-trial VCH
  session, fits all four candidate models by MCMC, and computes bridge-sampling
  log model evidence for the confusion matrix.
- **Outputs:**
  `param_recovery/prior_based_mcmc/results/sim{i}_{true_model}.csv` (gitignored).
- **Used by:** `prior_recovery_vch_mcmc.py` (compile + plot) and
  `prior_recovery_aic_bic.jl` (rescoring by AIC/BIC, same chains).

</details>

<details><summary><code>prior_recovery_vch_mcmc.py</code> — local, generate/compile/plot</summary>

- **Prerequisites:** `DO_GENERATE` to build the job array before the HPC run;
  `DO_COMPILE`/`DO_PLOT` after results are pulled back.
- **What it does:** builds the 2000-task dSQ array and prints transfer/submit/pull
  commands; after pulling, stacks per-simulation CSVs into one summary and draws
  the generative-vs-recovered scatter panels.
- **Outputs:**
  `param_recovery/prior_based_mcmc/results/prior_recovery_mcmc_summary.csv`,
  `param_recovery/prior_based_mcmc/figures/prior_recovery_scatter_{2,3}level_{empiric,nominal}.png`
  (vendored into `julia_outputs/` for the public repo).
- **Used by:** the four `prior_recovery_scatter_*.png` files are read by
  `04_visualizations/supplement/hgf_param_recovery_assembly.py`, which writes
  `param_recovery_scatterplots_all_models.png` — panel **a** of S3. Its confusion
  matrix uses bridge-sampling evidence, which the manuscript does **not** report;
  `prior_recovery_aic_bic.py`'s BIC version is the one that appears as S3b.

</details>

<details><summary><code>prior_recovery_aic_bic.jl</code> — HPC, patch step</summary>

- **Prerequisites:** `prior_recovery_vch_mcmc.jl` has already written the chains
  for this `(sim_index, true_model)`. Refits nothing — replays the original RNG
  stream to reconstruct the generative parameters, cross-checks them against the
  saved CSV (fatal error on mismatch beyond 1e-8), then scores the existing chains.
- **What it does:** computes AIC and BIC per candidate model from the saved MCMC
  chains, appending `aic_{model}`, `bic_{model}`, `aic_winner`, `bic_winner` columns
  to the existing per-simulation result CSV in place.
- **Outputs:** updates
  `param_recovery/prior_based_mcmc/results/sim{i}_{true_model}.csv` in place.
- **Used by:** `prior_recovery_aic_bic.py`.

</details>

<details><summary><code>prior_recovery_aic_bic.py</code> — local, plot</summary>

- **Prerequisites:** `prior_recovery_aic_bic.jl` has patched all 2000 result CSVs.
- **What it does:** recompiles the summary CSV from the per-simulation files (so it
  always reflects the latest columns) and draws row-normalized confusion matrices
  — rows are the generating model, columns are the model that won after inversion.
- **Outputs:**
  `param_recovery/prior_based_mcmc/figures/aic_bic/model_identifiability_confusion_bic.png`
  (and an `_aic.png` variant, not published) — vendored into `julia_outputs/`.
- **Used by:** `04_visualizations/supplement/hgf_6panel_assembly.py`, panel **b**
  of S3 — the manuscript reports the BIC version.

</details>

### Validation track B — model selection / RFX-BMS (S3c)

<details><summary><code>bms_vch.jl</code> — local (fast — MAP only, no MCMC)</summary>

- **Prerequisites:** Stage 1 staging CSVs. Optionally warm-starts its MAP fits
  from Stage 3's posterior medians if `results/vch/<model_type>/.../
  <id>_medians.csv` already exist — falls back to random-restart MAP otherwise.
- **What it does:** for every participant × every one of the four candidate
  models, computes a Laplace-approximate log model evidence at the MAP estimate.
  A subject whose Hessian isn't positive definite for any candidate model is
  dropped from BMS (their MCMC estimates from Stage 3 remain valid everywhere
  else).
- **Outputs:** `model_comparison/bms/results/<id>_lme.csv` (gitignored).
- **Used by:** `bms_vch.py`.
- **No seed is set** on the random-restart MAP fallback, so re-running this
  script can shift Ef/PXP/BOR slightly. Add a seed if you need an exact match
  to a prior run.

</details>

<details><summary><code>bms_vch.py</code> — local, run/BMS/plot</summary>

- **Prerequisites:** `DO_RUN=True` to invoke `bms_vch.jl` first (skips subjects
  already done); `DO_BMS`/`DO_PLOT` to run group-level RFX-BMS and draw bar charts.
- **What it does:** variational E/M on a Dirichlet over model frequencies →
  expected frequency (Ef), exceedance probability (XP), Bayesian Omnibus Risk
  (BOR), protected exceedance probability (PXP) (Stephan et al. 2009; Rigoux et
  al. 2014).
- **Outputs:** `model_comparison/bms/bms_summary.csv` (tracked — see *Data and
  paths*) and `model_comparison/figures/bms_ef_pxp.png` (not published as-is).
- **Used by:** `04_visualizations/supplement/hgf_bms_modified.py` re-renders
  `bms_summary.csv` with cosmetic changes only (title removed, axis labels
  capitalised) into `bms_ef_pxp_modified.png` — panel **c** of S3.

</details>

### Validation track C — posterior-based parameter recovery (S3d, S3e)

<details><summary><code>param_recovery_vch.jl</code> — HPC</summary>

- **Prerequisites:** Stage 3 output — reads each participant's own fitted
  posterior medians (`vch_nu`, `vch_beta`, `vch_omega[, vch_omega3]`) from the
  staged run as generative truth, and that participant's real stimulus sequence.
- **What it does:** simulates 10 synthetic response sequences per participant from
  their own fitted medians, then refits the HGF to each under identical priors and
  MCMC settings as Stage 2.
- **Outputs:** `param_recovery/results/{model_type}/` — `{id}_iter{i}.jls/.csv`,
  `{id}_iter{i}_medians.csv`, `{id}_iter{i}_rhats.csv` (gitignored).
- **Used by:** `param_recovery_vch.py`.

</details>

<details><summary><code>param_recovery_vch.py</code> — local, generate/compile/plot</summary>

- **Prerequisites:** `DO_GENERATE` before the HPC run; `DO_COMPILE`/`DO_PLOT`
  after pulling results. A run with any parameter's R-hat outside 0.9–1.1 is
  flagged `converged = False` in the summary CSV and excluded from every figure
  (the flag stays in the CSV even though the exclusion is invisible in the plots).
- **What it does:** compiles recovered medians against generative truth, computes
  per-run beta/nu posterior correlation, draws recovery scatter plots.
- **Outputs:** `param_recovery/figures/corr_gen_vs_rec_2level_empiric.png` (panel
  **d**), `param_recovery/figures/pair_beta_nu_2level_empiric.png` (panel **e**) —
  vendored into `julia_outputs/`.
- **Used by:** `04_visualizations/supplement/hgf_6panel_assembly.py`, directly —
  these two are read straight from `julia_outputs/param_recovery/figures/`, not
  through an intermediate assembly script.

</details>

### Validation track D — posterior predictive checks (S3f, S2b)

<details><summary><code>ppc_classic_vch.jl</code> — local</summary>

- **Prerequisites:** a medians CSV written by `ppc_classic_vch.py`'s `DO_RUN` step,
  which itself requires Stage 3 output (pulls posterior medians from the wide df).
- **What it does:** simulates one synthetic response sequence per participant at
  their posterior median parameters (no MAP, no covariance — deliberately
  conservative), computes per-condition and per-block detection rates for both the
  real and simulated sequences.
- **Outputs:** `param_recovery/ppc_classic/results/<model_type>/<id>_ppc_classic.csv`
  (gitignored).
- **Used by:** `ppc_classic_vch.py`'s `DO_COMPILE`/`DO_PLOT` steps.

</details>

<details><summary><code>ppc_classic_vch.py</code> — local, run/compile/plot</summary>

- **Prerequisites:** Stage 3 output. Run one `model_type` at a time — for
  `2level_nominal`, `PARAM_COLS` points at the *same* wide-df columns as
  `2level_empiric` (because `import_hgf_results_unified.py` writes nominal results
  under the empiric column names), so those columns hold whichever variant was
  imported most recently. Import → PPC → import the other, or the nominal figures
  are silently drawn from empiric estimates.
- **What it does:** writes the medians CSV, invokes `ppc_classic_vch.jl` once per
  model type, then plots observed vs. simulated group-mean detection rates with
  bootstrapped 94% CIs (10,000 resamples) plus per-participant spaghetti lines.
- **Outputs:**
  `param_recovery/ppc_classic/results/<model_type>/ppc_classic_<model_type>_all.csv`,
  `param_recovery/ppc_classic/figures/<model_type>/ppc_classic_<model_type>_{conditions,blocks}.png`
  — vendored into `julia_outputs/` for `2level_empiric` and `2level_nominal`.
- **Used by:** two downstream assemblies —
  `04_visualizations/supplement/hgf_ppc_assembly.py` stacks the `2level_empiric`
  conditions/blocks images into `ppc_2level_stacked.png` (panel **f** of S3);
  `04_visualizations/supplement/hgf_ppc_oos_assembly.py` reads the compiled
  `..._all.csv` for both `2level_empiric` and `2level_nominal` to build panel **b**
  of S2 from scratch.

</details>

---

## The task and the model

Participants viewed visual gratings at one of four contrast levels — 0 %, 25 %,
50 %, or 75 % of their individually QUEST-calibrated detection threshold — and
reported on each trial whether they detected a stimulus. Twelve blocks of 30
trials; 360 trials total.

The HGF models how a participant learns the *probability* of detecting the
stimulus over trials, and how that learned expectation combines with the current
stimulus to produce a percept.

### Architecture

```
[xvol]  ← log-volatility of xprob   (3-level only)
   ↓
[xprob] ← learned log-odds of detection probability
   ↓
[xbin]  ← binary detection response (0/1)
```

Each trial runs three steps (implemented in `helper_functions/create_agent.jl`):

1. **Learn.** The participant's *previous response* is fed to the HGF as its
   observation. This closed loop — the model learns from its own prior percept
   rather than from ground truth — is what distinguishes the CH-task HGF from a
   standard binary HGF.

2. **Perceive.** Blend the top-down prediction with the bottom-up stimulus,
   weighted by ν:

   ```
   xbin_pred = σ(xprob_posterior_mean)
   belief    = xbin_pred + 1/(1+ν) × (stimulus − xbin_pred)
             = ν/(1+ν) × xbin_pred  +  1/(1+ν) × stimulus
   ```

   ν → ∞ gives a purely prior-driven percept; ν = 0 a purely stimulus-driven one;
   ν = 1 weights them equally.

3. **Respond.** Pass the belief through a unit-square sigmoid whose steepness is
   set by β:

   ```
   P(yes) = 0.5 + 0.5 × tanh(β × (belief − 0.5))
   ```

### Parameters and priors

Priors are truncated normals, identical across every model variant so that no
variant is advantaged in model comparison. Centres come from posteriors of an
earlier CH-task HGF, widened to σ = 1.

| Symbol | Julia name | Prior | Constraint | Interpretation |
|---|---|---|---|---|
| β | `action_precision` | N(3.41, 1) | > 0.001 | Consistency of responses with the belief |
| ν | `prior_posterior_weight` | N(0.7265, 1) | > 0 | Weight on the prior vs. the likelihood |
| ω₂ | `xprob_volatility` | N(−5.1683, 1) | < −0.5 | Log tonic volatility of xprob; higher = faster updating |
| ω₃ | `xvol_volatility` | N(−6, 1) | < −0.5 | Log tonic volatility of xvol (3-level only) |

**MCMC:** NUTS, 4 chains × 1000 iterations, MAP initialisation, `MCMCSerial`.
R̂ is required to fall within 0.9–1.1 for every estimate entering any analysis;
`import_hgf_results_unified.py` sets failing parameters to `NaN`.

---

## Stimulus intensity: empirical vs. nominal

Each trial enters the HGF as a stimulus intensity in [0, 1] — the ground-truth
probability that a signal is present. Two conventions are used in the manuscript:

| Convention | Column | Values | Where used |
|---|---|---|---|
| **Empirical** (primary) | `empirical_condition` | {0.000, 0.418, 0.712, 0.899} | All primary results |
| **Nominal** (sensitivity) | `condition` | {0.000, 0.250, 0.500, 0.750} | Supplementary Table S5 |

The 25/50/75 % empirical values are hit rates from the non-hallucinator subset
(n = 29) of an independent, out-of-set COPE normative cohort — participants
reporting neither visual nor auditory hallucinations, so that the detection
likelihood reflects ordinary perception rather than a mixed clinical sample.
Derived by `out_of_set_data/empirical_condition_analysis.py` (Stage 0 above).

> **The 0 % condition is always exactly 0.0 — never the empirical false-alarm rate.**
> At 0 % contrast there is no stimulus, so the ground-truth detection probability
> is 0. Substituting a measured false-alarm rate would misrepresent the task
> structure to the HGF.

The empirical convention was chosen over the nominal one because QUEST-targeted
and observed detection probabilities diverged, the divergence replicated across
independent samples, detection rates were stable across the session, and PPCs
under the nominal convention departed visibly from observed behaviour. See
`out_of_set_data/README_empirical_condition_redetermination.md` for the full
derivation, and Supplementary Fig. S2.

### Where these values live

The mapping is hardcoded in four places, deliberately, so that no script depends
on a value it did not itself declare:

- `hgf_pipeline.py` → `VCH_CONDITION_TO_EMPIRICAL`, `VCH_STAGING_MAPPINGS`
- `model_fitting_singleagent_forarray.jl` → `VCH_CORRECTED_MAPPING`
- `ppc_classic_vch.jl` → `VCH_STIM_MAPPINGS`
- `prior_recovery_vch_mcmc.jl` → `EMPIRIC_INTENSITIES`

**If you change the mapping, change all four and re-stage.** As a safeguard,
`model_fitting_singleagent_forarray.jl` recomputes `empirical_condition` from the
`condition` column at runtime, ignoring whatever the staging CSV holds — so a
stale staging file can never silently propagate an old mapping into a new fit.

---

## Output files

Per-participant results land in `results/<model_type>/<project>/<timepoint>/<modality>/`:

| File | Contents |
|---|---|
| `<id>.jls` | Full MCMCChains object (Julia serialization) |
| `<id>.csv` | All posterior samples, flat CSV |
| `<id>_medians.csv` | Posterior median per parameter |
| `<id>_rhats.csv` | Gelman–Rubin R̂ per parameter |
| `<id>state_trajectories.csv` | Trial-by-trial HGF states at the posterior medians |

### Key state-trajectory columns

| Column | Description |
|---|---|
| `xbin_pred` | Predicted P(yes) *before* observing the response — the model's per-trial behavioural prediction |
| `xbin_pe` | Value prediction error: `response − xbin_pred` |
| `xprob` | Posterior mean log-odds of detection after the update — the core learned belief |
| `xprob_pred`, `xprob_pred_pp` | Prediction and its precision, before the trial's update |
| `belief` | The blended perceptual belief (step 2 above) |

**Precision-weighted prediction error** is the central derived quantity:

```
pwPE(t) = xbin_pe(t) / xprob_pred_pp(t)
```

The raw HGF output leads the response by one trial, so `import_hgf_results_unified.py`
shifts `pwPE` by −1 to align trial *t*'s pwPE with trial *t*'s response. The
alignment is validated on every run: pwPE should be positive when `response = 1`
and negative when `response = 0`.

The shift is applied to the per-trial state frame **before** it is either merged
into the long df or aggregated into the wide df, so both consumers see the same
series. This matters for the per-block columns: aggregating before the shift
would make `vch_pwPE_block_1..12` medians of a series one trial out of step with
the `pwPE` released in the long df. `vch_pwPE_median` and `vch_pwPE_mean` are
unaffected either way — a whole-series median or mean is invariant to a
one-trial shift. Clipping (`|pwPE| > 5` → `NaN`) is applied before the shift, so
those `NaN`s travel with the values they belong to.

### Wide-df column naming

| Pattern | Example | Contents |
|---|---|---|
| `vch_{param}` | `vch_nu` | Posterior median |
| `vch_{state}_median` / `_mean` | `vch_xprob_median` | Summary across all trials |
| `vch_{state}_block_{n}` | `vch_xprob_block_3` | Median within block *n* |
| `vch_{param}_log` | `vch_nu_log` | log1p-transformed parameter |

3-level columns carry a `_3lev` suffix. The `_nominal` variant shares column
names with the primary empiric run — only its filenames differ.

---

## Data and paths

### Layout

Every path is resolved from the scripts' own location, so the repository can be
cloned or moved anywhere without editing any file:

```
hppd_manuscript_public/
├── data/final/
│   ├── df_public_<date>.csv          ← wide df: input to the pipeline
│   ├── vch_master_public.csv         ← long df: one row per trial
│   └── backups/                      ← automatic pre-import copies
└── 02_hgf_modeling/                  ← this directory
    ├── data_n_cmnds/                 ← staged per-subject CSVs + job arrays (gitignored)
    ├── results/vch/<model_type>/<tp>/vch/    ← pulled HPC fits (gitignored)
    ├── model_comparison/
    │   └── bms/bms_summary.csv       ← TRACKED (see below); everything else gitignored
    ├── param_recovery/               ← recovery + PPC outputs (gitignored)
    └── julia_outputs/                ← TRACKED — vendored copies of the figures/CSVs
                                          S2 and S3 are assembled from (see below)
```

### `julia_outputs/` — why it exists

The panel images and CSVs that `04_visualizations/supplement/hgf_param_recovery_
assembly.py`, `hgf_ppc_assembly.py`, and `hgf_ppc_oos_assembly.py` read are
**vendored** into `julia_outputs/` — checked in rather than regenerated on clone
— because regenerating them means re-running the Julia validation pipeline
(Validation tracks A–D above) on a cluster. The directory layout under
`julia_outputs/` mirrors the Julia project's own output tree, so a file's path
there is the path it had when it was produced; the *Used by* line in each
script's entry above states exactly which vendored path that script's downstream
consumer reads. `model_comparison/bms/bms_summary.csv` is tracked the same way
but lives outside `julia_outputs/`, directly at that path, via a `.gitignore`
exception — see the comment above `!model_comparison/bms/bms_summary.csv` in
`02_hgf_modeling/.gitignore`.

`out_of_set_data/behavioral_data_OUT_OF_SET_with_metadata.csv`, vendored under
`julia_outputs/out_of_set_data/`, is trial-level task data keyed by `sudo_rec`, a
pseudonymised record id — no direct identifiers, the same class of data as the
shipped `vch_master_public.csv`. See *Data availability* below for why it isn't
tracked in git despite sitting in that directory.

### The primary, wide df (row per participant)

Every script here obtains the dataframe through one entry point:

```python
from hgf_pipeline import load_public_wide_df
df, path = load_public_wide_df("hppd_manuscript", require=["vch_nu", "vch_beta"])
```

`require` names the columns the caller depends on; anything missing raises
immediately, naming itself, instead of failing later with an opaque error.

### Data availability

Participant-level data beyond the released wide and long dfs is not distributed.

- Per-subject task CSVs (`data_n_cmnds/`) and fit results (`results/`) are
  regenerated by the pipeline from the released wide df.
- `out_of_set_data/behavioral_data_OUT_OF_SET_with_metadata.csv` — trial-level
  data from the independent COPE normative cohort — belongs to a separate study
  and is not redistributed (excluded via `.gitignore`, even though a local copy
  may be sitting in the directory on the machine that produced this repo). Its
  column dictionary (`behavioral_data_OUT_OF_SET_README.txt`) and all derived
  summaries and figures **are** included, so the derivation of the empirical
  stimulus intensities can be checked without the raw data.

---

## References

- Mathys, C. D. et al. (2014). Uncertainty in perception and the Hierarchical
  Gaussian Filter. *Frontiers in Human Neuroscience*, 8, 825.
- Stephan, K. E. et al. (2009). Bayesian model selection for group studies.
  *NeuroImage*, 46, 1004–1017.
- Rigoux, L. et al. (2014). Bayesian model selection for group studies — revisited.
  *NeuroImage*, 84, 971–985.
- Hoffman, M. D. & Gelman, A. (2014). The No-U-Turn Sampler. *JMLR*, 15, 1593–1623.
- Vehtari, A., Gelman, A. & Gabry, J. (2017). Practical Bayesian model evaluation
  using leave-one-out cross-validation and WAIC. *Statistics and Computing*, 27,
  1413–1432.
