# Computationally Characterizing Psychedelics' Persisting Perceptual Abnormalities

Analysis code for *A Computational and Behavioral Exploration of Serotonergic
Psychedelics' Persisting Perceptual Effects*.

**Preprint:** https://doi.org/10.64898/2026.08.10.743999

Everything needed to reproduce every figure, table and reported statistic in the
manuscript is here, starting from the de-identified analysis dataframes in
`data/final/`. 

---

## What you need

| | |
|---|---|
| **Python** | 3.12 — `pandas numpy scipy matplotlib seaborn statsmodels pingouin scikit-learn python-docx pillow svgutils` |
| **R** | 4.4 — `brms rstan DHARMa DHARMa.helpers bayestestR posterior tidyverse glue zoo distributions3` |
| **Julia** | 1.11 — only for re-deriving the HGF parameters (Step 1); `Project.toml` / `Manifest.toml` are pinned |
| **An HPC (if wanting to rerun in timely manner)** | Steps 1 and 2 fit thousands of MCMC models. Both stages generate SLURM job arrays; neither is practical on a laptop. |

---

## Before you run anything: three things to set

**1. Your cluster paths and login.** Every remote path derives from one constant
per file. 

| File | Constants |
|---|---|
| `03_hpc/generate_hpc_jobs.py` | `HPC_PARENT`, `HPC_USER`, `HPC_TRANSFER`, `HPC_LOGIN` |
| `03_hpc/compile_nonsp_results.py` | `HPC_BASE`, `HPC_HOST` |
| `03_hpc/compile_mediation_results.py` | `HPC_HOST` |
| `03_hpc/generate_nonsp_diagnostic_jobs.py` | `HPC_PARENT`, `HPC_LOGIN`, `HPC_TRANSFER`, `HPC_HOST` |
| `02_hgf_modeling/hgf_pipeline.py` | `HPC_BASE`, `HPC_USER`, `HPC_LOGIN_HOST`, `HPC_TRANSFER_HOST`, `JULIA_MODULE_*` |

`HPC_HOST` is an **SSH alias**, not a hostname — the pull scripts reuse an open
`ControlMaster` (`ssh -MNf <alias>`) rather than authenticating per file.

**2. The module lines in the generated job files.** `generate_hpc_jobs.py` writes
`module load foss/2022b && module load R/4.4.1-foss-2022b` into every job line and
sets `R_LIBS_USER=$HOME/R/4.4`. Change those to your cluster's module names for R.
`02_hgf_modeling/hgf_pipeline.py` writes its own module line from
`JULIA_MODULE_VCH` (default `Julia/1.11.4-linux-x86_64`) — change that to your
cluster's Julia module name.

**3. Nothing else.** All in-repo paths are relative to the repository root and
resolve from wherever you launch a script.

### The data

`data/final/df_public_<date>.csv` is the analysis dataframe, **already fully
prepared** — every derived column is present. Read it directly; there is no
preparation step to run. Scripts pick the most recent `df_public_*.csv` and print
which one they chose.

Two smaller files ship beside it:

| File | Used by |
|---|---|
| `vch_master_public.csv` | trial-level VCH data — HGF state trajectories (Figure 6), `results_narrative.py` |
| `df_recruit_public_<date>.csv` | de-identified recruitment export (raw data with unused participants but redacted PII) — `consort_diagram.py` (S9), `results_narrative.py` |

`modules/master_config.py` is the single source of truth for variable metadata:
brms family, normalization rule, plot labels, covariate sets, colour palettes. Any
question about how a variable is handled is answered there first.

---

## Run order

### Step 1 — HGF parameters *(optional — skip unless you are re-deriving them)*

**Nothing in Steps 2–5 needs this step.** The three HGF parameters the manuscript
uses — ν (prior weighting), β (decision precision) and ω₂ (belief-updating rate) —
are already columns of the shipped `df_public_*.csv`, and the panels behind
Supplementary Figures S2 and S3 are vendored in `02_hgf_modeling/julia_outputs/`.
A fresh clone reproduces every figure without running any Julia.

Run this only if you want to repeat the HGF inversion — on these data, or on new
data of your own. Note that we failed to set a seed so results will vary from exact numbers but
this should be negligible (make sure that HGF inversion succeeded for every participant though
-- you might need to run more iterations or for longer for a few fits)

```bash
python 02_hgf_modeling/hgf_pipeline.py          # stage 1: per-subject CSVs + job array
#   → transfer, then run the printed dSQ command on the cluster
#     (model_fitting_singleagent_forarray.jl)
python 02_hgf_modeling/import_hgf_results_unified.py   # stage 3: R-hat QC, merge back
```

> **Check R-hat before using anything you refit.** Every estimate entering an
> analysis must have R-hat within 0.9–1.1.
> `import_hgf_results_unified.py` runs that check and sets failing parameters to
> `NaN` rather than dropping them silently — so a participant whose chains did
> not converge simply disappears from downstream models unless you notice. Read
> the R-hat summary it prints, and **refit any model that failed** (more
> iterations, or a different initialisation) before treating its parameters as
> usable.

Model comparison and recovery (`bms_vch`, `prior_recovery_*`, `param_recovery_*`,
`ppc_classic_*`) are validation runs reported in the supplement. They do not feed
any parameter used in the main analyses.

### Step 2 — Fit the BRMs

One entry point builds every job array the manuscript needs — single-path
regressions and mediation models, across all covariate specifications.

We used the dSQ module on Yale's Bouchet cluster to handle job arrays. If that is
not available on yours, you will need to adapt: the generated `.txt` files are
just lists of shell commands, one model per line, and can be run however your
scheduler prefers — or locally, slowly.

**You never call the two `compile_*.py` scripts yourself.** `generate_hpc_jobs.py`
runs them for you, and passes them the exact DV / predictor / model-type filters
it just generated jobs from. Leave both files where they are — deleting either
one breaks step 4 below.

#### 1. Generate the jobs — and decline the compile prompt

```bash
python 03_hpc/generate_hpc_jobs.py
```

It writes the analysis dataframe, the job arrays, and a copy of the R scripts
into each job bundle, then prints the `rsync` → `dsq` → `sbatch` sequence.

It then asks whether to pull and compile results. **Answer `N`.** Nothing has
been submitted yet, so there is nothing on the cluster to pull. On `N` it prints
the two compile commands so you have them for step 4.

#### 2. Transfer and submit

Run the `rsync` commands it printed, then the `dsq` command.

**dSQ is two steps:** `dsq` does not submit anything — it prints an `sbatch`
command, which you then run yourself. Walltime is **per job**, not per array:
`-t 90:00` for single paths, `-t 360:00` for mediation.

#### 3. Wait

```bash
squeue -u <your-username>      # jobs still running
sacct  -u <your-username>      # what finished, and how
```

Thousands of models. Hours, not minutes.

#### 4. Compile, once the array is done

Open an SSH ControlMaster first — the compile scripts reuse it rather than
authenticating per file:

```bash
ssh -MNf bouchet
```

Then either re-run the generator and answer `y` this time:

```bash
python 03_hpc/generate_hpc_jobs.py      # regenerating the job files is harmless
```

or run the two commands it printed in step 1. Both routes do the same thing —
re-running the generator is simply the way to avoid retyping a long command.

Each script SSHes in, tarballs the result CSVs on the cluster, pulls one archive,
and extracts it locally. There is no separate retrieve step.

**Where the results land:**

| | Script | Output |
|---|---|---|
| **Single-path regressions** | `compile_nonsp_results.py` | `results/sensitivity_analyses_single_paths/existingresults_manuscript.csv` and `existingresults_manuscript_counterfactual.csv` — two flat tables, one row per model |
| **Mediation models** | `compile_mediation_results.py subset` | `results/{dv}/mediation_models/{model_name}/` — one directory per model (e.g. `results/hppd_binary/…`, `results/caps_vision/…`), holding the summary CSVs plus DHARMa and posterior-predictive PNGs |

Everything downstream reads `results/`. Nothing after this point touches the
cluster.

> **A model that did not converge writes no summary table.** That is deliberate —
> see § Convergence gate in `03_hpc/README.md`. A missing row or a blank
> heatmap cell means "did not converge", not "not tested".

### Step 3 — Main figures and tables

```bash
python 04_visualizations/0X_all_figures.py
```

Produces Figures 2–7, Table 1 and Supplementary Figure S6 — **and the summary CSVs
that Steps 4 and 5 read.** Run it before them. All of the nonparameteric tests are first run here. 

Section toggles at the top of the file (`RUN_HPPD_CAPS_FIGS`, `RUN_TABLES`,
`RUN_DESCRIPTIVE_FIGS`, `RUN_FIGURE_ASSEMBLY`, …) let you re-run one part.
`RUN_DIAGNOSTIC_COMPILATION` shells out to R and is slow; turn it off if you only
want figures.

**Figure 1 is a hand-drawn schematic** with no generating script. It ships as
`results/final_figures/figure_1.svg`.

### Step 4 — Supplementary figures and tables

Order-independent among themselves; all require Step 3.

```bash
cd 04_visualizations/supplement
python beta_sigmoid_creator.py                          # S1
python hgf_ppc_oos_assembly.py                          # S2
python hgf_bms_modified.py                              #      panel source for S3
python hgf_param_recovery_assembly.py                   #      panel source for S3
python hgf_ppc_assembly.py                              #      panel source for S3
python hgf_6panel_assembly.py                           # S3   (run the three above first)
python sensitivity_analyses.py                          # S4
python sensitivity_analyses_mediation.py all            # S5
python vch_beta_qc_scatter_supplement.py                # S7
python hardware_keydown_check.py                        # S8
python consort_diagram.py                               # S9
python caps_vision_confounds_spearman.py                # Table S1
python mann_whitney_table_hppd_binary.py                # Table S2
python regression_results_table_nominal_sensitivity.py  # Table S3 a, b
python mediation_results_table_nominal_sensitivity.py   # Table S3 c, d
python regression_results_table.py                      # Table S4
python mediation_results_table.py                       # Table S5
python spearman_table_caps_vision.py                    # Table S6
```

`fdr_correction.py` is imported by the two nonparametric tables, not run directly —
it is the single place Benjamini-Hochberg is applied, so both tables cannot drift.

Every command above is local and offline; none of them contacts the cluster.
`sensitivity_analyses_mediation.py all` aggregates the fitted mediation models
in `results/` and then draws; pass `heatmap` instead to redraw from the
aggregated CSV without re-reading the model directories.

### Step 5 — The Results narrative

```bash
python 05_results_narrative/results_narrative.py
```

Writes `results_narrative_output_editted.txt`: the manuscript's Results section,
every number computed from `results/` rather than transcribed.

---

## What produces what

Outputs are named for their place in the manuscript.

### Figures — `results/final_figures/`

| | Script | Source panel |
|---|---|---|
| **Figure 1** | *(hand-drawn)* | `figure_1.svg` |
| **Figure 2** | `0X_all_figures.py` | PPA history distributions |
| **Figure 3** | `0X_all_figures.py` + `caps_item_distributions_hppd_split.py` | CAPS item distributions by PPA |
| **Figure 4** | `0X_all_figures.py` | SP predictors → PPA / CAPS |
| **Figure 5** | `0X_all_figures.py` | VCH behaviour + mediation |
| **Figure 6** | `0X_all_figures.py` | VCH computations (HGF) + mediation |
| **Figure 7** | `0X_all_figures.py` | β, SDT, detection curves |

### Table — `results/descriptive/tables/`

| | Script |
|---|---|
| **Table 1** | `0X_all_figures.py` (`table_1.docx`) |

### Supplementary figures — `results/supplement/`

| | Script | Path |
|---|---|---|
| **S1** | `beta_sigmoid_creator.py` | `hgf_figures/supplementary_figure_s1.*` |
| **S2** | `hgf_ppc_oos_assembly.py` | `hgf_figures/supplementary_figure_s2.*` |
| **S3** | `hgf_6panel_assembly.py` | `hgf_figures/supplementary_figure_s3.*` |
| **S4** | `sensitivity_analyses.py` | `sensitivity_analyses/supplementary_figure_s4.*` |
| **S5** | `sensitivity_analyses_mediation.py` | `sensitivity_analyses_mediation/supplementary_figure_s5.*` |
| **S6** | `0X_all_figures.py` | `supplementary_figure_s6.*` |
| **S7** | `vch_beta_qc_scatter_supplement.py` | `vch_beta_qc_scatter_supplement/supplementary_figure_s7.*` |
| **S8** | `hardware_keydown_check.py` | `hardware_keydown_check/supplementary_figure_s8.*` |
| **S9** | `consort_diagram.py` | `consort_diagram/supplementary_figure_s9.*` |

### Supplementary tables — `results/supplement/tables/`

| | Script | File |
|---|---|---|
| **S1** | `caps_vision_confounds_spearman.py` | `supplementary_table_s1.docx` |
| **S2** | `mann_whitney_table_hppd_binary.py` | `supplementary_table_s2.docx` |
| **S3** | `regression_results_table_nominal_sensitivity.py` | `supplementary_table_s3ab.docx` |
| **S3** | `mediation_results_table_nominal_sensitivity.py` | `supplementary_table_s3cd.docx` |
| **S4** | `regression_results_table.py` | `supplementary_table_s4.docx` |
| **S5** | `mediation_results_table.py` | `supplementary_table_s5.docx` |
| **S6** | `spearman_table_caps_vision.py` | `supplementary_table_s6.docx` |

Every figure is written as `.png`, `.tiff` and (where vector) `.svg`; every table as
`.docx`, `.csv` and `.png`.

---

## Conventions worth knowing before you change anything

- **Gelman normalization divides by 2 SD**, not 1: `(x - mean) / (2 * sd)`.
- **94% HDI is the reported interval**, everywhere. Coefficient tables carry it
  as `hdi_lower_94` / `hdi_upper_94`, mediation tables as `lower_94_hdi` /
  `upper_94_hdi`. There is no toggle: reporting an equal-tailed 95% CI instead
  would mean changing the readers, not flipping a constant.
- **SP-user filter is always `psycheduse_yn == "Yes"`.**
- **Categorical covariates are never normalized in Python.** They are passed raw
  and factored in R. The list is derived from `master_config.py`, which generates
  `categorical_factor_vars.R` for both model scripts — do not declare factors
  anywhere else.
- **Model-type names are HPC directory names.** Two entries sharing a name
  overwrite each other's results on the cluster, silently.
- **brms strips underscores** from variable names in coefficient columns:
  `vch_bl_yes_0_normalized` → `vchblyes0normalized`.

`results/` and `data/` are gitignored. `results/` is created by the pipeline;
`data/final/` ships with the repository.

---

## Directory map

```
02_hgf_modeling/     Julia HGF fitting, model comparison, recovery, PPCs
03_hpc/             SLURM job generation, the two brms model scripts, result compilation
04_visualizations/  0X_all_figures.py (main figures + Table 1) and supplement/
05_results_narrative/  results_narrative.py
modules/            master_config.py, shared plotting and normalization helpers
data/final/         the de-identified analysis dataframes
results/            all output (created by the pipeline)
```
