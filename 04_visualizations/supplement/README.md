# 04_visualizations/supplement/

Scripts for generating figures and analyses associated with the supplementary
materials of the HPPD manuscript.

All outputs go to `results/supplement/` (gitignored — regenerate by running the
scripts). Subdirectories within `results/supplement/` mirror the script names.

**Output files are named for their place in the manuscript**, not for the script
that made them: `supplementary_figure_s{n}.*` and `supplementary_table_s{n}.*`.
The subdirectory still names the analysis, so `supplementary_figure_s8.png` lives
in `hardware_keydown_check/`. Renaming a figure means editing one output stem in
one script; the mapping is tabulated below.

**All figure scripts** (those that write to `results/supplement/` but NOT to
`results/supplement/tables/`) produce both a `.png` and a `.tiff` at the same
DPI alongside each other. The `.tiff` files are for journal submission; the
`.png` files are for manuscript drafting and Google Docs preview.

---


## Redacted inputs

`qc_notes` (RA free text) and `phone_number` are **not** in the shipped
recruitment CSV. The scripts here that used them now read frozen values:

- `modules/qc_redacted_categories.py` holds the QC-failure category
  memberships as record_id sets, snapshotted from the original notes.
  `consort_diagram.py` looks them up rather than pattern-matching, and uses
  them only to fill the quality-control box of the diagram — it writes no QC
  table and no participant-level QC artifact. Published counts are unchanged
  (39 / 12 / 24 / 6 of 69).
- `phone_cc_234` replaces `phone_number`; it is the +234 country-code boolean,
  the only thing the audit ever derived from the number.
- `task_data_prltask_present`, `continue_date_passed` and
  `timesincesurveystart_gt14` replace the blob, the date and the timestamp.

That classification is therefore a **replay of a frozen result**, not a live
computation; every other count in the participant flow still computes normally.
A data refresh therefore requires the frozen module to be regenerated alongside
the recruitment CSV.


## Editable table requirement (Scientific Reports)

**Scientific Reports requires editable tables (Word or TeX), not images.**
All table-generating scripts in this directory produce three output formats:

| Format | Purpose |
|--------|---------|
| `.docx` | **Primary submission artifact.** Formatted Word table with booktabs-style borders, bold/italic section headers, amber-highlighted P-cells, and proper subscript rendering (e.g. *r*<sub>rb</sub>). Imports cleanly into Google Docs. |
| `.csv` | Machine-readable sidecar with raw numeric values. Encoded as UTF-8 with BOM for Excel/Google Sheets compatibility. |
| `.png` | Visual reference for manuscript drafting. Not submitted. |

No table script writes `.html` or `.tex`. `results/supplement/tables/` should
contain exactly the six supplementary tables, in those three formats each.

If you add a new table script that outputs to `results/supplement/tables/`,
it **must** also produce `.docx` and `.csv` companions. Use the shared
`_docx_helper.py` module for Word table generation.

**`_docx_helper.py`** — shared module exporting `save_docx_tables()`. Takes
the same `display_df` / `row_meta` / `prob_rows` data that drives the PNG
renderer. Supports subscript rendering via `subscript_map` parameter and
landscape orientation via `landscape=True`.

### Which script feeds which supplementary table

| Supp. table | Panels | Script |
|---|---|---|
| S1 | — | `caps_vision_confounds_spearman.py` |
| S2 | — | `mann_whitney_table_hppd_binary.py` |
| S3 | a, b | `regression_results_table_nominal_sensitivity.py` |
| S3 | c, d | `mediation_results_table_nominal_sensitivity.py` |
| S4 | a, b | `regression_results_table.py` |
| S5 | a, b | `mediation_results_table.py` |
| S6 | — | `spearman_table_caps_vision.py` |

### Which script feeds which supplementary figure

| Supp. figure | Script | Output |
|---|---|---|
| S1 | `beta_sigmoid_creator.py` | `hgf_figures/supplementary_figure_s1.*` |
| S2 | `hgf_ppc_oos_assembly.py` | `hgf_figures/supplementary_figure_s2.*` |
| S3 | `hgf_6panel_assembly.py` | `hgf_figures/supplementary_figure_s3.*` |
| S4 | `sensitivity_analyses.py` | `sensitivity_analyses/supplementary_figure_s4.*` |
| S5 | `sensitivity_analyses_mediation.py` | `sensitivity_analyses_mediation/supplementary_figure_s5.*` |
| S6 | `04_visualizations/0X_all_figures.py` | `results/supplement/supplementary_figure_s6.*` |
| S7 | `vch_beta_qc_scatter_supplement.py` | `vch_beta_qc_scatter_supplement/supplementary_figure_s7.*` |
| S8 | `hardware_keydown_check.py` | `hardware_keydown_check/supplementary_figure_s8.*` |
| S9 | `consort_diagram.py` | `consort_diagram/supplementary_figure_s9.*` |

**S6 is written by `0X_all_figures.py`, not by a script in this directory.** It
is the SP use-characteristics figure, assembled from
`results/descriptive/sp_table_distributions.png` in that script's
linked-figures block. It is the one entry in `linked_figures` routed to
`results/supplement/` instead of `results/final_figures/` — see
`_SUPPLEMENT_FIGS` there. Everything else in `results/supplement/` comes from
this directory.

The `.csv` output of each table script is the authoritative numeric record; the
`.docx` and `.png` are renderings of it.

---

## Files

> **Parallel-maintenance warning — sensitivity analysis scripts**
>
> Two sensitivity pipeline scripts live in this directory and both produce outputs under
> `results/supplement/`:
>
> | Script | Analyses covered |
> |--------|-----------------|
> | `sensitivity_analyses.py` | Single-path (nonsp predictor) models |
> | `sensitivity_analyses_mediation.py` | Mediation models (A/B/C′/NIE paths) |
>
> These scripts share the same visual conventions (color maps, cell format, flagging
> logic, compound-figure layout) and were written to be parallel.  **If you ever change
> any display logic in one — color thresholds, heatmap layout, flagging criteria, cell
> text format — you MUST replicate the change in the other.** Do not let them diverge.

| Script | Output dir | Description |
|--------|-----------|-------------|
| `consort_diagram.py` | `results/supplement/consort_diagram/` | **Supplementary Figure S9.** CONSORT-style participant flow diagram: screening → final analytic sample → analytic subsamples. Counts every participant at exactly one stage so each arrow subtracts, and asserts chain closure at runtime. Owns both sets of exclusion-reason bullets — the screening loop and the QC classification both run here, off masks it builds itself, with no sibling script and no intermediate CSV. Writes the figure only. See details below. |
| `caps_vision_confounds_spearman.py` | `results/supplement/tables/` | **Supplementary Table S1.** Spearman correlations between `caps_vision` and all demographic/clinical table variables (continuous, ordinal, and binary), plus four 6-month substance-use count variables. See details below. |
| `beta_sigmoid_creator.py` | `results/supplement/hgf_figures/` | Plots the unit-square sigmoid function P(u_t) = b_t^β / (b_t^β + (1−b_t)^β) for β ∈ {0.25, 1.75}, illustrating how action precision controls belief-to-response mapping. Pure math — reads nothing. See details below. |
| `hgf_param_recovery_assembly.py` | `results/supplement/hgf_figures/` | Assembles four prior-based HGF parameter recovery scatter plots into a stacked 4-panel figure with pair titles and row labels. See details below. |
| `hgf_ppc_assembly.py` | `results/supplement/hgf_figures/` | Stacks the two PPC participant-spread plots (conditions + blocks), removing the legend from the bottom panel and replacing the x-axis label "Condition (signal strength)" with "Contrast Intensity". See details below. |
| `hgf_ppc_oos_assembly.py` | `results/supplement/hgf_figures/` | **Supplementary Figure S2.** Three-panel figure (a, b, c): out-of-set empirical condition-by-condition detection rates (SP users), stacked PPC conditions (empiric top / nominal bottom) with shared y-axis label, and out-of-set non-hallucinator group per-block detection rates at 75% contrast with 94% CI. Writes `supplementary_figure_s2.*`. The four-panel variant, whose panel d showed two exemplar QUEST staircases, is archived at `06_submission/superseded_figures/`. |
| `hgf_bms_modified.py` | `results/supplement/hgf_figures/` | Re-renders the RFX-BMS Ef + PXP bar chart from `bms_summary.csv` with the suptitle removed and x-axis labels capitalised. Intermediate figure consumed by `hgf_6panel_assembly.py`. See details below. |
| `hgf_6panel_assembly.py` | `results/supplement/hgf_figures/` | Assembles six HGF validation panels into a 2-row × 3-column publication figure (`supplementary_figure_s3.png`). Reads five source PNGs from `julia_hgf_ch/` plus the modified BMS figure. See details below. |
| `sensitivity_analyses.py` | `results/supplement/sensitivity_analyses/` | **Supplementary Figure S4.** Draws the compound single-path sensitivity heatmap (both DVs, 12 model types) from `existingresults_manuscript_counterfactual.csv`. Figure only — it submits no jobs and pulls nothing from the cluster; both are done by `03_hpc/generate_hpc_jobs.py`. See details below. |
| `sensitivity_analyses_mediation.py` | `results/supplement/sensitivity_analyses_mediation/` | **Supplementary Figure S5.** Aggregates the fitted mediation models into `compiled_sensitivity_mediation.csv`, then draws the NIE compound heatmap (both DVs across the sensitivity covariate sets). Submits no jobs — `03_hpc/generate_hpc_jobs.py` fits every type in its `CUSTOM_MED_TYPES`. `PATHS_TO_PLOT` is `['NIE']`; add `'A path'`, `'B path'`, `"C' path"` back to render those too. See details below. |
| `regression_results_table.py` | `results/supplement/tables/` | Two-panel compound table of Bayesian regression results for all forest-plot predictors (sp_predictors, vch_behavior, vch_computations) and both primary outcomes (hppd_binary, caps_vision). See details below. |
| `mediation_results_table.py` | `results/supplement/tables/` | Two-panel publication table of all mediation results (4 mediators × 2 outcomes, `nice_covariates_spusers`). Rows: a/b/c′ paths + NIE/NDE/Total/PMed per mediator. Columns: β 94%HDI P(β≠0) Δ 94%HDI P(Δ≠0). See details below. |
| `regression_results_table_nominal_sensitivity.py` | `results/supplement/tables/` | Two-panel companion to `regression_results_table.py`, restricted to the two HGF parameters estimated with **nominal, QUEST-derived** detection probabilities (`vch_nu_nominal`, `vch_beta_nominal`). Supplies panels **a** and **b** of Supplementary Table S5. |
| `mediation_results_table_nominal_sensitivity.py` | `results/supplement/tables/` | Two-panel companion to `mediation_results_table.py` for the same two nominal-prior mediators. Supplies panels **c** and **d** of Supplementary Table S5. |
| `mann_whitney_table_hppd_binary.py` | `results/supplement/tables/` | Single-panel publication table of Mann-Whitney U test statistics for `hppd_binary` (PPA History), grouped by IV type. Columns: PPA(−) Mdn [Q1,Q3] · PPA(+) Mdn [Q1,Q3] · U · r_rb · p · p_FDR+ · p_FDR− · N. Reads the canonical `boxplot_grid.csv` sidecars (SP users only, n≈186) generated by `0X_all_figures.py`; the two FDR columns come from `fdr_correction.py`. See details below. |
| `spearman_table_caps_vision.py` | `results/supplement/tables/` | Single-panel publication table of Spearman correlations between `caps_vision` and every forest-plot predictor, grouped by IV type. Columns: ρ · p · p_FDR+ · p_FDR− · N. The `sp_predictors` block is the AGE-CONTROLLED partial ρ; the VCH blocks are zero-order, and the table carries no marker saying so. Companion to `mann_whitney_table_hppd_binary.py` — keep the two parallel. See details below. |
| `fdr_correction.py` | `results/supplement/fdr_correction/` | Benjamini-Hochberg FDR correction of every frequentist test in the manuscript for `sp_predictors` / `vch_behavior` / `vch_computations`. Two families (9 `caps_vision` Spearman correlations; 9 `hppd_binary` Mann-Whitney U tests), each corrected across all nine tests (`p_FDR+`) and within each predictor block (`p_FDR-`). Single source for the FDR columns in both supplement tables. CSV only. See details below. |
| `vch_beta_qc_scatter_supplement.py` | `results/supplement/vch_beta_qc_scatter_supplement/` | **Supplementary Figure S7.** 5-panel scatter figure (2×3 grid, 1 slot empty): vch_beta vs. curated QC metrics (SP users). Points colored by vch_beta (lower = electric blue). Y-axis labels from `VARIABLE_REGISTRY.plot_label`. Likert y-ticks on effort_qc and distraction_qc. No titles. Composite scores excluded (experimental; reported nowhere). See details below. |
| `hardware_keydown_check.py` | `results/supplement/hardware_keydown_check/` | **Supplementary Figure S8.** Does the participant's display covary with `d_prime_overall` or `vch_threshold`? Display hardware is hand-coded from free text into a 3-level display class (`monitor_check_operationalized_final`), the same coding used as the HPC hardware-control covariate — which ships as a column of `data/final/df_public_*.csv`, so this script no longer writes a lookup CSV for the job generators. One **stacked 2×1 figure at half a Nature page width** — **a** d′ over **b** threshold — Kruskal-Wallis per panel, per-group n printed above every box, every pairwise post-hoc in the CSV. Defines its own `hardware_boxplot_grid()` rather than reusing `multipanel_boxplot_grid`. See details below. |
| `diagnostics/create_mediation_diagnostic_compilation.py` | `{model_dir}/diagnostic_compilation.png` | Per-model diagnostic compilation figure (MCMC traces, PP checks, DHARMa). Calls `_compile_diagnostics_helper.R` via subprocess. See details below. |
| `diagnostics/_compile_diagnostics_helper.R` | `{model_dir}/compiled_*.png` | Companion R script that generates 13 cached diagnostic PNGs from a `.RData` fit object. Called by `create_mediation_diagnostic_compilation.py` — do not run directly. |
| `diagnostics/compile_mediation_diagnostic_pdfs.py` | `results/supplement/diagnostics/mediation_diagnostics_{dv}.pdf` | Batch-generates `diagnostic_compilation.png` for all parseable mediation models, then assembles per-DV multi-page PDFs. See details below. |
| `diagnostics/compile_single_path_diagnostic_pdfs.py` | `results/supplement/diagnostics/single_path_diagnostics_{dv}.pdf` | Assembles per-DV diagnostic PDFs for nonsp single-path brms models from HPC-generated compilation PNGs. See details below. |

## `beta_sigmoid_creator.py`

### Purpose

Plots the unit-square sigmoid function used in the HGF response model:

    P(u_t) = b_t^β / (b_t^β + (1 − b_t)^β)

Two curves are drawn for β ∈ {0.25, 1.75} (β=1 identity line excluded), illustrating
how the action precision parameter β controls how deterministically the model's belief
maps to a response probability. Higher β pushes probabilities toward 0 and 1; lower β
compresses them toward 0.5.

Staggered dotted reference lines are drawn at x = 0.3 and x = 0.7 for each curve:
vertical lines from y = 0 to the curve height, and horizontal lines from x = 0 to the
x reference position, in the respective curve color. Axis labels "P(Target | Cue)" (x)
and "P(Reporting 'Yes')" (y) are always shown. All text uses Arial font.

**Config block at top of file** (all easy to change):

| Variable | Default | Purpose |
|----------|---------|---------|
| `BETAS` | `[0.25, 1.75]` | β values to plot (one curve each); β=1 omitted |
| `COLORS` | blue → dark | Curve colors; index matches `BETAS` order |
| `X_REFS` | `[0.3, 0.7]` | x positions for staggered reference lines |
| `REF_LW_THICK` | `3.5` | Dotted linewidth for light-blue (low β) reference lines |
| `REF_LW_THIN` | `1.2` | Dotted linewidth for dark-gray (high β) reference lines |
| `SHOW_LEGEND` | `False` | Toggle per-β legend |
| `SHOW_TICKS` | `False` | Toggle x/y tick marks and labels |
| `DPI` | `600` | Output resolution |

The "increasing β" annotation uses a gradient arrow (electric blue → dark grey) drawn
via `LineCollection` in axes-fraction coordinates.

### Reads

Nothing — pure math (no data files required).

### Output

```
results/supplement/hgf_figures/supplementary_figure_s1.png    (manuscript drafting)
results/supplement/hgf_figures/supplementary_figure_s1.tiff   (journal submission)
```

Output directory is gitignored — regenerate by running the script.

### How to run

```bash
/usr/local/bin/python3.12 04_visualizations/supplement/beta_sigmoid_creator.py
```

Run from any directory (paths resolve relative to this file).

### Common things to change

- `MODEL_VARIANTS`: the sensitivity model types rendered as columns, in display
  order. Names must match the HPC results directory names, which are the entries
  of `BASE_MODELS` in `03_hpc/generate_hpc_jobs.py` — keep the two lists in step.
- `HEATMAP_SOURCE_CSV`: the single CSV every heatmap column is read from.
- `CANONICAL_MODEL_TYPE` / `SECOND_CANONICAL_TYPE` / `THIRD_CANONICAL_TYPE`: the
  three leftmost columns. They select rows out of `HEATMAP_SOURCE_CSV` rather
  than being read from a second file.
- `DHARMA_ALPHA` (default 0.05) and the flagging logic in
  `_diagnostic_status_cf()`: keep in sync with `sensitivity_analyses_mediation.py`.

The per-DV heatmaps (`hppd_binary_sensitivity_heatmap`,
`caps_vision_sensitivity_heatmap`) and the single-path diagnostic PDFs this
script used to assemble are gone. The diagnostic PDFs the supplement actually
uses come from `diagnostics/compile_single_path_diagnostic_pdfs.py`, driven by
`0X_all_figures.py`, and land in `results/supplement/diagnostics/`.

---

## `sensitivity_analyses_mediation.py`

> **Parallel-maintenance:** this script mirrors `sensitivity_analyses.py`
> (single-path). If you change display logic in one, replicate it in the other.
> See the warning block at the top of this file.

### Purpose

Runs sensitivity analyses for all Bayesian **mediation models**, asking whether
the key mediation effects (A path: SP→mediator; B path: mediator→DV; C′: direct
SP→DV; NIE: natural indirect effect) are robust to covariate set choice.

Produces four compound figures — one per path type — each with one panel per DV
(`hppd_binary` top, `caps_vision` bottom).  Each panel is a heatmap with:

- **Rows** = (spvar, mediator) combinations
- **Columns** = 4 sensitivity covariate sets + 1 canonical reference column
- **Cell text** = `mean ± hdi_half` on the response scale (posterior mean — see
  `POINT_ESTIMATE_COL` in `modules/master_config.py`)
- **Red border / asterisk** = DHARMa diagnostic flag OR MCMC flag (Rhat > 1.01,
  ESS < 400, or divergent transitions > 0)

The canonical reference column (`nice_covariates_spusers`) is read from the
primary results tree (`results/{dv}/mediation_models/`) — no SSH required for it.

### Stages

The script has two stages, toggled by `RUN_*` flags in the CONFIG block at the
top. Job submission is **not** one of them — `03_hpc/generate_hpc_jobs.py` fits
every type in its `CUSTOM_MED_TYPES`, which must stay in step with
`SENSITIVITY_MED_TYPES` here.

| Stage | Flag | CLI arg | What it does |
|-------|------|---------|-------------|
| 1. Compile | `RUN_COMPILE` | `compile` | Reads the per-model CSVs already in `results/{dv}/mediation_models/` and writes the long-format `compiled_sensitivity_mediation.csv`. Local only — it reports any model missing from `results/` and carries on without it. |
| 2. Heatmap | `RUN_HEATMAP` | `heatmap` | Reads that CSV and draws one compound heatmap per entry in `PATHS_TO_PLOT` — currently `['NIE']`, the supplementary figure. |

Neither stage contacts the cluster. Pulling fitted models back is
`03_hpc/compile_mediation_results.py subset`; run that first if the compile
stage reports models missing from `results/`.

### CONFIG block (key variables)

```python
# Sensitivity covariate types (9 total, all different from the primary
# pipeline's nice_covariates_spusers).  Keep in sync with MODEL_VARIANTS in
# sensitivity_analyses.py — the two scripts are maintained in parallel.
SENSITIVITY_MED_TYPES = [
    'empirical_covariates_spusers',
    'nice_covariates',                    # full sample  -> SECOND_CANONICAL slot
    'nice_covariates_spusers_iqr',        # IQR fence    -> THIRD_CANONICAL slot
    'age_control_spusers',
    'true_univariate_spusers',
    'nice_covariates_beta_spusers',
    'nice_covariates_spusers_nonan_caps',
    'drugs_month_spusers',                # + 6 past-month drug binaries
    'drugs_trimmed_month_spusers',        # + depressants/cannabis/stimulants
]

# Primary-pipeline covariate type — read from results/{dv}/mediation_models/ (no HPC pull needed)
CANONICAL_MED_TYPE = 'nice_covariates_spusers'

# Analyses: one dict per (spvar, mediator, dv) combination
MED_ANALYSES = [...]  # auto-built from _MED_MEDIATORS × spvar/dv pairs

# Which mediation paths to heatmap — all 4 by default
PATHS_TO_PLOT = ['A path', 'B path', "C' path", 'NIE']
```

### Reads

- `results/{dv}/mediation_models/{canonical_model_name}/` — canonical reference
  (4 CSV files per model, read locally)
- `results/{dv}/mediation_models/{model_name}/` — the sensitivity model types,
  read locally. Anything absent is reported and skipped; pull it with
  `03_hpc/compile_mediation_results.py subset`.

### Outputs

| File | Description |
|------|-------------|
| `results/sensitivity_analyses_mediation/compiled_sensitivity_mediation.csv` | Long-format compiled results (one row per model × path). Columns: `dv, spvar, mediator, cov_type, path, mean, hdi_low, hdi_high, p_above_0, p_below_0, p_direction, dharma_flagged_dv, dharma_flagged_med, dharma_flagged_any, mcmc_flagged, max_rhat, min_ess, num_divergents, model_name`. The point-estimate column is named for `POINT_ESTIMATE_COL` (`modules/master_config.py`), currently `mean`. |
| `results/supplement/sensitivity_analyses_mediation/a_path_compound_heatmap.png` | A path (SP → mediator) compound figure — hppd_binary + caps_vision panels |
| `results/supplement/sensitivity_analyses_mediation/a_path_compound_heatmap.tiff` | A path — TIFF for journal submission |
| `results/supplement/sensitivity_analyses_mediation/b_path_compound_heatmap.png` | B path (mediator → DV) compound figure |
| `results/supplement/sensitivity_analyses_mediation/b_path_compound_heatmap.tiff` | B path — TIFF for journal submission |
| `results/supplement/sensitivity_analyses_mediation/c_prime_path_compound_heatmap.png` | C′ path (direct SP → DV) compound figure |
| `results/supplement/sensitivity_analyses_mediation/c_prime_path_compound_heatmap.tiff` | C′ path — TIFF for journal submission |
| `results/supplement/sensitivity_analyses_mediation/supplementary_figure_s5.png` | NIE (natural indirect effect) compound figure |
| `results/supplement/sensitivity_analyses_mediation/supplementary_figure_s5.tiff` | NIE — TIFF for journal submission |

All output directories are gitignored — regenerate by running the script.

### How to run

```bash
cd hppd_manuscript_public

# Aggregate the fitted models, then draw
/usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses_mediation.py all

# Redraw from the aggregated CSV without re-reading the model directories
/usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses_mediation.py heatmap
```

### Common things to change

- `SENSITIVITY_MED_TYPES`: add/remove sensitivity covariate sets.  A set added
  here must also exist in `BASE_COVARIATE_SETS` **and** have a matching branch in
  `get_covs()` in `03_hpc/generate_hpc_jobs.py` — that function resolves
  covariate strings through a hardcoded `if/elif` chain, not a dict lookup, and
  raises `ValueError: Unknown cov_type` otherwise.  Names must
  match the `cov_type` keys used in `generate_hpc_jobs.py`.
- `CANONICAL_MED_TYPE`: the reference covariate set (leftmost column).  Must
  already exist in `results/{dv}/mediation_models/`.
- `MED_ANALYSES`: list of `{spvar, mediator, dv}` dicts.  Add new combinations
  here; the generate stage will add corresponding HPC jobs automatically.
- `PATHS_TO_PLOT`: subset to `['NIE']` if you only want indirect-effect heatmaps.
- `HPC_MED_RESULTS` (in CONFIG): update if the HPC directory structure changes.
- `_CMAP`, `_NEUTRAL_COLOR`, flagging thresholds in `_plot_med_heatmap_panel()`:
  keep in sync with `sensitivity_analyses.py` if visual conventions change.
- DHARMa column names in `_read_equation_diagnostics()`: mediation summary CSVs
  use `dharma_ks_pval` (not `dharma_uniformity_pval`).  Do not confuse with the
  nonsp-predictor pipeline.

---

## `regression_results_table.py`

Produces a two-panel supplementary compound table (`regression_results_table.png`)
that compiles Bayesian regression results for the same predictor groups and DVs
shown in the primary forest plots of `0X_all_figures.py`.

### Panel layout

| | Description |
|---|---|
| **Panel a** | 9 predictor rows (grouped by SP Use Patterns / VCH Task Behavior / HGF Estimates). Columns: standardized β section (PPA History + CAPS Vision, each with β̂, HDI₉₄–, HDI₉₄+, P(β̂≠0)) and marginal contrast Δ section (PPA History + CAPS Vision, each with Δ̂, P(Δ̂≠0)). N on far right. |
| **Panel b** | Identical to Panel a, with a third column group in the β section for `hu_{predictor}_normalized`: the hurdle (zero-inflation) component of the `caps_vision` hurdle-negbinom model. |

### Data sources

| Section | Source CSV | Filter | Key columns |
|---------|-----------|--------|-------------|
| β̂, HDI₉₄, P(β̂≠0) | `existingresults_manuscript.csv` | `covariates == CANONICAL_COV` | `Estimate`, `hdi_lower_94`, `hdi_upper_94`, `prob_above_0`, `prob_below_0` |
| Δ̂, P(Δ̂≠0) | `existingresults_manuscript_counterfactual.csv` | `cov == CANONICAL_COV` | `mean` (via `point_estimate()`), `hdi_lower_94`, `hdi_upper_94`, `prob_above_0`, `prob_below_0` |

`CANONICAL_COV` is `nice_covariates_spusers` (`regression_results_table.py:63`);
the two CSV paths are set at lines 55–56. Note the filter column differs between
the files — `covariates` in the beta CSV, `cov` in the counterfactual one. The
Δ̂ point estimate goes through `point_estimate()` from `modules/master_config.py`
rather than reading `mean` directly, so both halves of a table row summarize the
posterior the same way.

### Symbol conventions

- **β̂** = posterior mean standardized coefficient (Gelman-normalized predictor: divided by 2 SD)
- **Δ̂** = marginal response-scale counterfactual contrast: E[Y|X=mean+1 SD] − E[Y|X=mean]
- **P(·≠0)** = max(P(·>0), P(·<0)) — posterior probability that the effect is nonzero
- Amber-highlighted cells: P ≥ 0.90 (matching the forest-plot highlight threshold)

### Future columns under consideration

The following columns are available in `existingresults_manuscript.csv` and
may be added in a later revision without changing the data source:

- **Estimate Error** (`Est.Error` column) — posterior SD of the coefficient
- **Effective Sample Size** (`Bulk_ESS` column) — MCMC chain quality indicator

### Common things to change

- `CANONICAL_COV` (default `'nice_covariates_spusers'` — the primary manuscript
  model): swap to any of the 10 model types present in `existingresults_manuscript.csv`
  to show a different specification. Confirm the type exists in **both** the beta and
  counterfactual CSVs first — they do not carry identical model-type sets.
- `IVTYPES`: subset to fewer predictor groups if the table becomes too tall.
- `PROB_THRESHOLD` (default 0.90): amber-highlight cutoff.
- `BETA_COLS` / `DELTA_COLS`: add or rename statistic columns.
- `OUTPUT_DPI` (default 400): reduce to 150 for faster iteration during editing.

### How to run

```bash
cd 04_visualizations
/usr/local/bin/python3.12 supplement/regression_results_table.py
```

Output:
- `results/supplement/tables/regression_results_table.png`
- `results/supplement/tables/regression_results_table.csv` — raw numeric values (panel, dv, section, predictor, beta, HDI bounds, probabilities, delta, N; hu columns for caps_vision)

---

## `mediation_results_table.py`

Publication-style compound table summarising **all mediation analysis results**
for the manuscript's canonical covariate set (`nice_covariates_spusers`).

### Panel layout

| Panel | DV | SP predictor |
|-------|----|-------------|
| **a** | `hppd_binary` (PPA History) | `psychedelic_age` (`spage`) |
| **b** | `caps_vision` (CAPS Vision Score) | `avg_life_dose` (`avgdose`) |

Both panels share the same column structure:

```
Path | β | β 94% HDI | P(β≠0) | Δ/Δ_med | Δ 94% HDI | P(Δ≠0)
```

### Mediators (4 per panel, each a section header)

| Key | Display label |
|-----|--------------|
| `vchbeta` | VCH Decision Noise (β) |
| `vchnu` | VCH Top-down Bias (ν) |
| `vchrate` | VCH Response Rate |
| `vchthreshold` | VCH 75% Detection Threshold |

### Rows per mediator section

| Row label | Source | β? | Δ? |
|-----------|--------|----|----|
| a | `path_coefficients_summary.csv` path=`a` + `path_counterfactual_summary.csv` "A path" | ✓ | ✓ |
| b | same, path=`b` + "B path" | ✓ | ✓ |
| c′ | same, path=`c_prime` + "C' path" | ✓ | ✓ |
| b (hu) | `path_coefficients_summary.csv` path=`b_hu` | ✓ | — (Panel b only) |
| c′ (hu) | `path_coefficients_summary.csv` path=`c_hu` | ✓ | — (Panel b only) |
| NIE | `mc_mediation_summary.csv` prefix "NIE" | — | ✓ |
| NDE | `mc_mediation_summary.csv` prefix "NDE" | — | ✓ |
| Total | `mc_mediation_summary.csv` prefix "TE" | — | ✓ |
| PMed | `mc_mediation_summary.csv` prefix "PMed" | — | ✓ (proportion, 0–1 scale) |

`b (hu)` and `c′ (hu)` appear only in Panel b (`caps_vision`, hurdle-negbinom);
there are no counterfactual Δ values for the hurdle paths.
`NIE`/`NDE`/`Total`/`PMed` are response-scale quantities with no normalised β.

### Data sources

Model directories: `results/{dv}/mediation_models/{dv}_{spvar}_{mediator}_nice_covariates_spusers/`

| File | Columns used | Notes |
|------|-------------|-------|
| `path_coefficients_summary.csv` | `mean`, `hdi_lower_94`, `hdi_upper_94`, `prob_above_0`, `prob_below_0` | β and 94% HDI for each path |
| `path_counterfactual_summary.csv` | `mean`, `hdi_low`, `hdi_high`, `p_above_0`, `p_below_0` | Δ for A/B/C′ paths (note: column names differ from path_coefficients — `hdi_low`/`hdi_high`, not `hdi_lower_94`/`hdi_upper_94`) |
| `mc_mediation_summary.csv` | `mean`, `hdi_low`, `hdi_high`, `p_above_0`, `p_below_0` | NIE/NDE/TE/PMed; matched by `effect.str.startswith()` |

### Two posterior summaries share the Δ column — read this before changing a number

The Δ column carries **two different summaries**, decided by which CSV the row
came from:

| Rows | Source CSV | Reported point estimate | Constant |
|---|---|---|---|
| a, b, c′ | `path_counterfactual_summary.csv` | posterior **mean** | `POINT_ESTIMATE_COL` |
| NIE, NDE, Total, PMed | `mc_mediation_summary.csv` | posterior **median** | `MC_EFFECT_POINT_ESTIMATE_COL` |

The four MC effects are produced by `hpc_mediation.R`'s Monte-Carlo integration
over the mediator's posterior predictive. For a `caps_vision` DV that predictive
passes through `hurdle_negbinomial`'s log link and is unbounded in the mediator,
so a handful of extreme draws can supply essentially the whole mean — before this
was fixed, Supplementary Table S7 b reported `NDE = +63646555.766`, and NIE and
Total rendered as `—` because their means were non-finite. The median, the 94% HDI
and the direction probabilities are all quantile-based and stable. Full rationale,
including the cost (means are exactly additive, NIE + NDE = TE; medians are not),
is at `modules/master_config.py:60-111`.

Consequences for this script:

- `get_cf_effect(df, prefix, mc_integrated=...)` — the caller decides, per call
  site, never the function by sniffing `source`. The MC call sites (display rows
  **and** the CSV builder) pass `mc_integrated=True`.
- The header of the Δ point-estimate column is **`Δ/Δ_med`** — matplotlib mathtext
  in the PNG, a real Word subscript in the `.docx` via `_docx_helper`'s
  `subscript_map` (`DELTA_SUBSCRIPT_MAP` at the top of the script). Only that
  column is marked; the HDI and P columns are identical under either summary and
  carry no caveat.
- The `.csv` sidecar has a **`Δ summary`** column holding `mean` / `median` per
  row, so the machine-readable output says on its face which is which.

Both `median` and `mean` are written to every CSV by the R pipeline, so the choice
can be revisited from results already on disk without refitting. Readers call
`point_estimate()` rather than naming a column, so a CSV written before the mean
columns existed raises a named error identifying the file and the required refit
instead of silently reporting the other summary.

`mediation_results_table_nominal_sensitivity.py` is maintained in parallel with
this script and has the identical treatment. So does the NIE row in
`sensitivity_analyses_mediation.py`, whose heatmap renders `p_direction` only —
the point estimate there reaches `compiled_sensitivity_mediation.csv` and nothing
else, but it is a median for the same reason.

### Styling

Matches `regression_results_table.py` / `generate_publication_table.py`:
Arial font, bold-italic section headers, TB structural lines only, amber P-cells ≥ 0.90.

### Common things to change

- `COV_SET` (default `"nice_covariates_spusers"`): covariate set suffix in directory names.
- `MEDIATORS` / `MEDIATOR_LABELS`: add/remove mediators or update display labels.
- `PROB_THRESHOLD` (default 0.90): amber-highlight cutoff.
- `pred_frac` inside `draw_panel()` (default 0.13): fraction of axis width for the
  Path column. Increase if any label is truncated.
- `DAT_H`, `SEC_H`, `COL_H` (in `main()`): row heights in inches; reduce to compress
  the figure if it becomes too tall for a journal supplement page limit.

### How to run

```bash
cd 04_visualizations
/usr/local/bin/python3.12 supplement/mediation_results_table.py
```

Output:
- `results/supplement/tables/mediation_results_table.png`
- `results/supplement/tables/mediation_results_table.csv` — raw numeric values (panel, dv, sp_predictor, mediator, path, beta, HDI bounds, probabilities, delta)

---

## `mann_whitney_table_hppd_binary.py`

Publication-style single-panel table of Mann-Whitney U test statistics for
the `hppd_binary` (PPA History) outcome, grouped by IV type — the frequentist
complement to the Bayesian forest-plot panels in Figures 1–2.

### Panel layout

One panel ("a") with 9 predictor rows grouped under three bold-italic section
headers:

| Section header | Predictors |
|----------------|-----------|
| SP Use Patterns | Age of First Use, Lifetime SP Uses (Count), Avg. Dose Used |
| VCH Task Behavior | 75% Detection Threshold, Hit Rate (75% Contrast), VCH Rate |
| HGF Estimates | Prior Weighting (ν), Decision Noise (β), Learning Rate (ω) |

### Columns

```
Predictor | PPA(−) Mdn [Q1, Q3] | PPA(+) Mdn [Q1, Q3] | U | r_rb | p
          | p_FDR+ | p_FDR- | N
```

- `Mdn [Q1, Q3]`: combined string, 3 decimal places each value
- `U`: Mann-Whitney U statistic, 1 decimal place
- `r_rb` (rendered as italic-r subscript-rb via mathtext): rank-biserial
  correlation, signed (+/−), 3 decimal places
- `p`: uncorrected p-value, 3 decimal places; "< .001" for very small values
- `p_FDR+` (bold upright p, subscript FDR+): Benjamini-Hochberg across all 9
  frequentist tests for `hppd_binary`
- `p_FDR-` (bold upright p, subscript FDR−): Benjamini-Hochberg within the predictor's
  own `iv_type_dict` block only (m = 3)
- `N`: total sample size for that predictor's test

The two FDR headers are **bold upright**, not italic, matching the plain `p`,
`U` and `N` headers. Matplotlib mathtext ignores
`cell.set_text_props(weight="bold")`, so the weight has to live inside the
mathtext itself — hence `$\mathbf{p}_{\mathbf{FDR}+}$` rather than
`$p_{\mathrm{FDR}+}$`. The `.docx` segments carry no `italic` key for the same
reason; bold is inherited from the header-row call in `_docx_helper._write_cell`.
`r_rb` is deliberately left italic and unbolded — that is standard notation for
the rank-biserial correlation and predates these columns.

Both FDR columns are imported from `fdr_correction.py` via
`family_frame("hppd_binary")`; this script never runs BH itself. See the
`fdr_correction.py` section below for the family definitions.

**`p_FDR-` is not uniformly smaller than `p_FDR+`.** BH is a step-up procedure,
so a test can be rescued by a smaller p-value ranked below it in the larger
family. `vch_beta` is the live example: p = .0229 gives p_FDR- = .0687 within
`vch_computations` (m=3), but p_FDR+ = .0626 across all nine, because
`vch_threshold` at p = .0278 pulls it down under monotonicity enforcement.

### Styling

Identical to `regression_results_table.py` panel a: Arial font,
bold column headers, bold-italic section-header rows with TB lines only,
right-divider on Predictor column, thick bottom closing line, no highlighting.

### Data sources

```
results/hppd_binary/{ivtype}/data_visualization/summary_results/boxplot_grid.csv
```
for `ivtype ∈ {sp_predictors, vch_behavior, vch_computations}`.

These CSVs are the sidecar files written by `multipanel_boxplot_grid()` when
`0X_all_figures.py` runs (requires `RUN_HPPD_CAPS_FIGS = True`). They contain
one row per predictor with test = `mann_whitney`, drawn from `df_sp_plot`
(SP users only, n ≈ 186). Re-run `0X_all_figures.py` to regenerate them if
the underlying data changes.

### Common things to change

- `IVTYPES`: subset or reorder the three IV-type groups.
- Formatting functions (`fmt_mdn_iqr`, `fmt_u`, `fmt_r`, `fmt_p`): adjust
  decimal places or p-value threshold for "< .001" notation.
- `_col_widths` inside `draw_panel()`: adjust per-column width fractions if
  any column content is clipped at a different figure width.
- `fig_w` in `main()` (default 17.5 in): increase if Mdn [Q1,Q3] columns are
  truncated; decrease to match a specific journal column width. It is tied to
  `_col_widths`: at `mdn_frac = 0.18` each Mdn column needs
  `fig_w ≥ 3 / (0.18 × 0.97) ≈ 17.2` to hold 23 characters.
- `OUTPUT_DPI` (default 400): reduce to 150 for faster iteration.

### How to run

```bash
cd 04_visualizations
/usr/local/bin/python3.12 supplement/mann_whitney_table_hppd_binary.py
```

Output:
- `results/supplement/tables/mann_whitney_table_hppd_binary.png`
- `results/supplement/tables/mann_whitney_table_hppd_binary.csv` — raw numeric values (section, predictor, PPA(-)/(+) medians and quartiles, U, r_rb, p, p_FDR+, p_FDR-, N)
- `results/supplement/tables/mann_whitney_table_hppd_binary.docx`

The `.docx` renders 13x9 — the two FDR columns (`p_FDR+`, `p_FDR-`) are extra
relative to a 13x7 layout. The `.csv` is the authoritative numeric output.

---

## `spearman_table_caps_vision.py`

Publication-style single-panel table of Spearman correlations between
`caps_vision` (CAPS visual items endorsed) and every forest-plot predictor,
grouped by IV type. Direct companion to `mann_whitney_table_hppd_binary.py`,
which does the same job for `hppd_binary`.

> **Parallel-maintenance warning.** These two table scripts share renderer
> style, formatting helpers, column conventions and FDR source. If you change
> display logic in one — column widths, cell format, section styling, the FDR
> columns — replicate it in the other. Do not let them diverge.

### Panel layout

One panel with 9 predictor rows under three bold-italic section headers
(`SP Use Patterns`, `VCH Task Behavior`, `HGF Estimates`), in `iv_type_dict`
order.

### Columns

```
Predictor | ρ | p | p_FDR+ | p_FDR- | N
```

- `ρ`: Spearman correlation, signed (+/−), 3 decimal places. The header is
  rendered as an italic Greek rho in the .docx via `SUBSCRIPT_MAP`.
- `p`: uncorrected p-value, 3 decimal places; "< .001" for very small values
- `p_FDR+` (bold upright p, subscript FDR+): Benjamini-Hochberg across all 9
  frequentist tests for `caps_vision`
- `p_FDR-` (bold upright p, subscript FDR−): Benjamini-Hochberg within the predictor's
  own `iv_type_dict` block only (m = 3)
- `N`: pairwise-complete sample size for that correlation

### Correlation type by block

The three `sp_predictors` rows are **partial** Spearman correlations
controlling for `age_v2` (`correlation_grid_age_control.csv`). The six VCH rows
are **zero-order** (`correlation_grid.csv`), because no age-controlled variant
is generated for those blocks.

This matches the manuscript — Fig. 1b shows the age-controlled panel and
`05_results_narrative/results_narrative.py` quotes the age-controlled values.
The table legend states which rows are which, and the sidecar `.csv` carries a
`Correlation Type` column per row.

### Data sources

Everything — ρ, p, N, and both FDR columns — comes from
`fdr_correction.family_frame("caps_vision")`. This script reads no results CSV
and runs no BH itself, so the grid mapping above lives only in
`CAPS_GRID_BY_IVTYPE` in `fdr_correction.py`.

The underlying grids are sidecars written by `correlation_matrix_plot()` when
`0X_all_figures.py` runs with `RUN_HPPD_CAPS_FIGS = True`. Re-run it to
regenerate them if the underlying data changes.

### Styling

Identical to `mann_whitney_table_hppd_binary.py`: Arial, bold column headers,
bold-italic section-header rows with TB lines only, right-divider on the
Predictor column, thick bottom closing line, no highlighting.

### Common things to change

- `_col_widths` inside `draw_panel()` (default `[0.40] + [0.12] * 5`): must sum
  to 1.0.
- `fig_w` in `main()` (default 11.0 in). At this width the longest label,
  "Contingency Belief Evolution Rate (ω)" (37 chars), overhangs its cell
  slightly — the same behaviour the Mann-Whitney table already has. Raise to
  ~12.1 to fit it fully. The `.docx` autofits and wraps regardless.
- `fmt_rho` / `fmt_p`: decimal places and the "< .001" threshold. These are
  deliberately byte-identical to the Mann-Whitney script's helpers.
- `OUTPUT_DPI` (default 400).

### Common pitfalls

- **`ρ` values print with variable decimal places.** `_strip()` removes trailing
  zeros, so −0.040 renders as "-0.04" while −0.156 renders in full. This is
  inherited from the Mann-Whitney formatter (`fmt_r`) and is consistent across
  both tables; change both or neither.
- **`p_FDR+` and `p_FDR-` can round to the same displayed value.** `vch_beta`
  is .0497 in both, printing as "0.050" in each column — visually at the α=.05
  line while actually just under it. If that matters for the reviewer response,
  raise the decimal places in `fmt_p`.

### How to run

```bash
cd 04_visualizations
/usr/local/bin/python3.12 supplement/spearman_table_caps_vision.py
```

Output:
- `results/supplement/tables/spearman_table_caps_vision.png`
- `results/supplement/tables/spearman_table_caps_vision.csv` — raw numeric values (section, predictor, correlation type, rho, p, p_FDR+, p_FDR-, N)
- `results/supplement/tables/spearman_table_caps_vision.docx`

---

## `fdr_correction.py`

### Purpose

Benjamini-Hochberg FDR correction of every frequentist test the manuscript
reports for the three canonical predictor blocks (`sp_predictors`,
`vch_behavior`, `vch_computations`).

Two independent families, one per outcome, nine predictors each:

| Family | Tests | n |
|---|---|---|
| `caps_vision` | 9 Spearman correlations | 130 (sp), 113 (vch) |
| `hppd_binary` | 9 Mann-Whitney U tests | 186 (sp), 160 (vch) |

Each family is corrected **two ways**:

| Pass | m | Contents | Shown in tables as |
|---|---|---|---|
| `_m9` | 9 | All nine frequentist tests for that outcome. | `p_FDR+` |
| `_m3` | 3 | Frequentist tests within one `iv_type_dict` block, run separately per block. | `p_FDR-` |

`FDR_PLUS_COL` and `FDR_MINUS_COL` at the top of the module name the two columns
the tables render, so both table scripts are pinned to the same passes.

BH is run with `statsmodels.stats.multitest.multipletests(method="fdr_bh")`,
which returns step-up adjusted p-values with monotonicity already enforced, so
`p_fdr_m9` / `p_fdr_m3` compare directly against alpha = 0.05.
`fdr_correction_verify.py` checks that against an independent reimplementation,
against R's `p.adjust`, and against a published known-answer dataset.

### `caps_vision` x `sp_predictors` uses the AGE-CONTROLLED partial Spearman

Two grids exist for that block. `correlation_grid_age_control.csv` is the one
the manuscript figure shows (Fig. 1b) and the one
`05_results_narrative/results_narrative.py` quotes. The zero-order
`correlation_grid.csv` is **not** used for `sp_predictors`. The two `vch_*`
blocks have no age-controlled variant, so they use `correlation_grid.csv`.
`CAPS_GRID_BY_IVTYPE` at the top of the script is the single place this mapping
lives.

### `p_FDR-` is not uniformly smaller than `p_FDR+`

This reads as a bug and is not. BH adjusts by **m/i** — family size over rank —
not by m, so a smaller family does not imply a smaller adjusted p; that
intuition holds for Bonferroni and fails for a step-up procedure. See the
`fdr_correction_verify.py` section below for the two mechanisms, worked through
on the live numbers and confirmed independently in R.

### Reads

```
results/caps_vision/sp_predictors/data_visualization/summary_results/
    correlation_grid_age_control.csv
results/caps_vision/{vch_behavior,vch_computations}/data_visualization/
    summary_results/correlation_grid.csv
results/hppd_binary/{ivtype}/data_visualization/summary_results/
    boxplot_grid.csv                                    (test == "mann_whitney")
```

All are looked up through `_require()` / `_one_row()`, which raise rather than
substituting a near-miss file or silently taking the first of several matching
rows.

### Imported by the table scripts

`mann_whitney_table_hppd_binary.py` and `spearman_table_caps_vision.py` both
call `family_frame("<outcome>")` and read `FDR_PLUS_COL` / `FDR_MINUS_COL`
rather than re-deriving anything. That keeps the family definitions, the
age-controlled-grid mapping, and the BH settings in exactly one place, and
guarantees the uncorrected `p` printed in a table is the same number the
correction was computed from.

Importing this module runs `os.chdir(_HERE)` as a side effect. All callers sit
in this directory and chdir there themselves, so it is a no-op — but keep that
in mind before importing it from elsewhere.

### Output

```
results/supplement/fdr_correction/fdr_correction.csv
```

18 rows (2 families x 9 tests). Columns:

| Column | Meaning |
|---|---|
| `family` | `caps_vision` or `hppd_binary` |
| `test_type` | `spearman` or `mann_whitney` |
| `test_detail` | e.g. `partial Spearman \| age_v2` |
| `ivtype`, `ivtype_label` | predictor block |
| `predictor`, `predictor_label` | variable and its `dv_to_lab_short` label |
| `n` | sample size for that test |
| `effect_name`, `effect` | `rho` / `U` and its value |
| `p_raw`, `p_tails` | uncorrected p and its tail convention |
| `sig_raw` | `p_raw < 0.05` |
| `m_m9`, `rank_m9`, `p_fdr_m9`, `sig_fdr_m9` | `p_FDR+` |
| `m_m3`, `rank_m3`, `p_fdr_m3`, `sig_fdr_m3` | `p_FDR-`; ranks restart at 1 in every block |
| `source_file` | path the row came from, relative to `results/` |

Written UTF-8 with BOM so Excel / Google Sheets render the Greek glyphs.
Output directory is gitignored — regenerate by running the script.

### How to run

```bash
cd 04_visualizations
/usr/local/bin/python3.12 supplement/fdr_correction.py
```

Prints a per-family table and the surviving-test counts, then writes the CSV.

### Common things to change

- `FDR_METHOD`: `"fdr_bh"` (Benjamini-Hochberg). Switch to `"fdr_by"` for
  Benjamini-Yekutieli if arbitrary-dependence control is wanted; BY is the more
  conservative choice given that predictors within a block are correlated.
- `ALPHA`: 0.05.
- `CAPS_GRID_BY_IVTYPE`: which Spearman grid each block reads (see above).
- `IVTYPES`: adding a fourth block changes both family sizes; the `_m9` / `_m3`
  column names are then misleading and should be renamed in the same commit.

### Common pitfalls

- **Changing `IVTYPES` or a family definition silently changes two submission
  tables.** Both table scripts import from here. Re-run all three scripts, plus
  `fdr_correction_verify.py`, after any edit.
- **No `.docx` / `.png` companion.** This script writes to
  `results/supplement/fdr_correction/`, not `results/supplement/tables/`, so the
  editable-table requirement at the top of this README does not apply. If the
  output is ever promoted to a supplementary table, it must gain `.docx` and
  `.png` outputs via `_docx_helper.py` at that point.
- **Re-run after `0X_all_figures.py`.** The Spearman and Mann-Whitney p-values
  come from the `data_visualization/summary_results/` sidecars, which
  `0X_all_figures.py` regenerates.

## `vch_beta_qc_scatter_supplement.py`

### Purpose

Publication-quality supplementary scatter figure testing whether vch_beta
(decision precision) correlates continuously with task quality-control metrics.
If low-beta participants were disengaged, QC metrics should be worse for them —
this figure shows they are not.

### Scientific context

The vch_beta finding (low beta → higher HPPD/CAPS outcomes) could in principle
be an artifact of task disengagement rather than a genuine perceptual trait.
This figure presents seven scatter plots examining that possibility.  The
direction of associations — where significant — runs opposite to the disengagement
concern: high-beta participants show more same-response-train behavior, worse
threshold calibration, and faster RT during streaks, not low-beta participants.

### Sample

SP users only (`psycheduse_yn == "Yes"`), n = 186 (n = 160 with JSON-decoded
task metrics; n = 151 with self-report QC items).

### Panels (2 × 3 grid; 1 slot empty)

| Row | Col | Variable | n | ρ | p | sig |
|-----|-----|----------|---|---|---|-----|
| 0 | 0 | `effort_qc` | 151 | −0.159 | 0.051 | ~ |
| 0 | 1 | `distraction_qc` | 151 | −0.027 | 0.740 | ns |
| 0 | 2 | `n_timeouts_total` | 160 | −0.098 | 0.216 | ns |
| 1 | 0 | `n_timeout_trials` | 160 | −0.105 | 0.186 | ns |
| 1 | 1 | `threshold_empiric_v_nominal` | 160 | +0.435 | <0.001 | *** |
| (1, 2) | — | *(empty)* | — | — | — | — |

**Note:** The streak-based composite scores (`z_composite_rt`, `z_composite_accuracy`)
are **experimental** — they use a multiplicative + directional formula still under
development. They are reported nowhere and are intentionally excluded from this
supplement figure; the columns are retained only so their definition is documented.

### Design

- **No titles.** Y-axis label = `VARIABLE_REGISTRY[col]['plot_label']`, wrapped
  at 28 characters via `textwrap.fill()`.
- **Likert y-ticks** on `effort_qc` and `distraction_qc` (1–5 scale):
  Strongly Disagree / Disagree / Neither Agree nor Disagree / Agree / Strongly Agree.
- **Point color** = vch_beta, via `electric_blue_palette` (reversed): lower = electric blue.
- **Spearman ρ and p** annotated upper-right of each panel.
- **Shared colorbar** placed outside the grid on the right edge.
- `VARIABLE_REGISTRY` plot_labels for all 5 variables are defined in
  `modules/master_config.py` (added alongside this script).

### Reads

- `most_recent_public_df()` (`modules/data_prep.py`)
- `task_data_vch_short_psychedelic_bl` (JSON embedded in df, decoded per-participant)
- `VARIABLE_REGISTRY`, `electric_blue_palette` (`modules/master_config.py`)

### Output

```
results/supplement/vch_beta_qc_scatter_supplement/supplementary_figure_s7.png
results/supplement/vch_beta_qc_scatter_supplement/supplementary_figure_s7.tiff
```

Output directory is gitignored (covered by global `results/` entry).

### Common things to change

- **`PANELS`** list in the script: add, remove, or reorder panels here. The
  plot_label is always fetched from `VARIABLE_REGISTRY` at render time —
  update the label there, not in the script.
- **`YLABEL_WRAP`** (default 28): increase to allow longer single-line y-labels.
- **`LIKERT_COLS`**: set of column names that receive Likert y-ticks.  Add any
  new 1–5 Likert variable here; the Likert label dict is `LIKERT_TICKS`.
- **`figsize`** (default `(N_COLS * 4.8, N_ROWS * 4.2)`): scale up for higher-
  resolution publication output.
- **`dpi`** (default 200): increase to 300 for final journal submission.

### How to run

```bash
/usr/local/bin/python3.12 04_visualizations/supplement/vch_beta_qc_scatter_supplement.py
```

Run from repo root.  Output: `results/supplement/vch_beta_qc_scatter_supplement/supplementary_figure_s7.{png,tiff}`

---

## `hardware_keydown_check.py`

### Purpose

Reviewer response. A reviewer asked whether the participant's **display** —
uncontrolled in an unsupervised online task — covaries with the two behavioural
measures we report, `d_prime_overall` and `vch_threshold`.

Display hardware was collected as a REDCap **free-text** field
(`monitor_check`), which cannot be tested directly. It is hand-coded into a
3-level display class, `monitor_check_operationalized_final`:

| Level | What it is |
|---|---|
| **Mac** | integrated Apple displays — consistently LED-backlit IPS LCD panels |
| **External Monitor** | a display separate from the machine driving it, and so the class with the most variable viewing distance, panel type and setup |
| **Windows/Other Laptop** | everything else |

This is the only grouping of the free text that corresponds to a physical
property plausibly able to produce a display difference, and the only one with
enough participants in every cell to test. It is also the coding carried into
the models as the hardware-control covariate (model type
`nice_covariates_spusers_hardware_control`), so the figure and the covariate
describe the same variable.

**Analysis rule.** Kruskal-Wallis across the three display classes, per DV,
with pairwise Mann-Whitney post-hoc contrasts (Bonferroni-adjusted within
panel). Whether a display difference is carried into the reported analyses is
decided on that omnibus test and its post-hoc contrasts.

```
row a   d_prime_overall  ×  monitor_check_operationalized_final
row b   vch_threshold    ×  monitor_check_operationalized_final
```

### Derivation

The coding collapses a finer-grained hand coding of the free-text field: `Mac Laptop` + `Mac Desktop/Monitor`
→ **Mac**; `Windows/Other Laptop` → **Windows/Other Laptop**;
`Standalone Monitor` + `TV used as monitor` → **External Monitor**;
`No Response` → NaN.

Five responses are then corrected by name, enforcing one rule the finer coding
does not: **the level describes the display, not the machine.** A response
naming a separate external display codes as External Monitor whatever machine
drives it; an Apple desktop codes as Mac:

| Raw string | finer coding | `_final` |
|---|---|---|
| `Mac 22" desktop` | Windows/Other Laptop | **Mac** |
| `apple desktop` | Windows/Other Laptop | **Mac** |
| `touchscreen laptop` | Standalone Monitor | **Windows/Other Laptop** |
| `Macbook Pro with a 27" 4K monitor` | Mac Laptop | **External Monitor** |
| `Macbook pro … Acer Nitro XF273 …` | Mac Laptop | **External Monitor** |

`Apple Studio Monitor` and `4K 27" Monitor on macOS 26` stay **External
Monitor** by the same rule — Apple machines, separate displays.

This is provenance, not live code: the script reads the finished
`monitor_check_operationalized_final` column and derives nothing. Recorded here
so the coding can be audited without opening the export that built it.

#### Where the column is read from

`monitor_check_operationalized_final` is **a materialised column of
`data/final/df_public_*.csv`**, read from there by every consumer — this script
and the HPC job generator alike. Nothing in this repository derives it, so the
figure and the `nice_covariates_spusers_hardware_control` covariate are the same
variable by construction rather than by two derivations that have to agree.

**To change the coding:** change it where the dataframe is built and re-export.

The script raises on a missing column, a non-numeric DV, or the two DVs no
longer sharing the same missing rows — no silent adaptation. It also raises if a
participant has a blank `monitor_check` but a non-null display class, which
would mean the coding and the raw field had come apart.

### Figure format — journal sizing, fonts, panel letters

The figure is built **at its final printed size**: it is never fed through
`modules/figure_assembly.py`, so nothing rescales it after the fact.

**Layout.** A stacked single column — `d_prime_overall` on top (**a**),
`vch_threshold` below (**b**). Both rows share the grouping variable, so only
panel **b** carries x-tick labels and the x-axis title
(`xaxis_bottom_row_only=True`). Repeating them on panel **a** would say the
same thing twice and cost vertical space the half-page width cannot spare.

**Width.** `FIG_WIDTH_IN` = **half** `JOURNAL_DOUBLE_COL_MM` (183 mm, the
Scientific Reports / Nature full-page measure defined in `0X_all_figures.py`)
= 91.5 mm ≈ 3.602″. `PANEL_HEIGHT_IN` = 2.6″ per row, so the saved figure is
3.6″ × 5.2″ (1080 × 1560 px at 300 dpi).

> **`bbox_inches='tight'` is off for this figure** (`bbox_tight=False`).
> Tight bbox crops to the artists, so the saved width is *not* the figsize
> width — which would defeat sizing to a journal measure. `tight_layout()` is
> what keeps the labels on the canvas instead. Verify after any layout change:
> the PNG must be 1080 px wide at 300 dpi.

**Fonts.** Arial, declared as a sans-serif *stack*
(`font.family='sans-serif'` + `font.sans-serif=['Arial','Helvetica','DejaVu
Sans']`) exactly as in `0X_all_figures.py`, so a machine without Arial degrades
predictably. Point sizes are the manuscript targets from that file, applied
directly rather than back-computed through `compute_source_fontsize()` (there
is no assembly scaling to cancel):

| Constant | pt | Source constant in `0X_all_figures.py` | Used for |
|---|---|---|---|
| `FONT_AXIS_LABEL` | 9.2 | `TARGET_AXIS_LABEL` | x- and y-axis labels |
| `FONT_TICK_LABEL` | 7.0 | `TARGET_TICK_LABEL` | tick labels, both axes |
| `FONT_SIG_MARKER` | 24.0 | `TARGET_SIG_MARKER` | post-hoc asterisks (drawn at 0.5×) |
| `FONT_ANNOT` | 7.0 | — (small-text tier, = tick size) | omnibus stat block |
| `FONT_GROUP_N` | 7.0 | — (small-text tier, = tick size) | the `n = …` row |
| `FONT_PANEL_LABEL` | 20.0 | `FIGURE_LABEL_FONTSIZE` | the bold **a** / **b** |

**Panel letters.** Drawn by `figure_assembly._add_panel_label()` — the same
helper every assembled manuscript figure uses (Arial bold, black, `va='top'`,
upper-left of the axes at `x=0.01, y=0.99`) — imported rather than
reimplemented so these labels cannot drift from the rest of the paper's.

**Point-unit geometry was halved** with the width. Marker sizes and line widths
are in points and therefore do *not* shrink with the figure: the values
inherited from the 14″-wide exploratory layout buried the boxes at journal
width. `strip_size` 4 → 2.5, `line_width` 3 → 1.5, `mean_marker_size` 14 → 7;
the mean-marker edge and the post-hoc bracket are now derived as
`line_width / 2` rather than pinned at 1.5.

> **Restoring the wide exploratory schemes.** The journal sizing assumes 3
> groups. Schemes 1–2 (11–17 levels) need `fig_width=None` (re-enables the
> group-count-driven auto-widening), `xaxis_bottom_row_only=False` (their rows
> use *different* grouping variables), `tick_fontsize=None` (re-enables the
> auto-shrinking tick size), and the original point-unit geometry above.

### Statistics

- **Kruskal-Wallis H** across the three display classes, per DV.
- **Effect size.** ε² = H / (N − 1).
- **Post-hoc.** Every pairwise Mann-Whitney is computed and written to the
  summary CSV with a raw and a **Bonferroni-adjusted** p-value (adjusted within
  panel, over all three pairs). Brackets are drawn on the figure only when the
  omnibus test is significant.
- The summary CSV carries a `test` column on every row, so which test produced
  a given row is never implicit.

### Why a new function instead of `multipanel_boxplot_grid`

`hardware_boxplot_grid()` is defined **inside this script** (not in
`modules/visualization_helpers_parts/`) and is visually identical to
`multipanel_boxplot_grid` — translucent box, per-group strip plot, mean marker,
manually drawn quartile rules, no spines, no ticks. It differs in five ways the
hardware layout requires:

1. **Every panel keeps its own y-axis and y-label.** The parent hides the
   y-axis for every column but the first, which assumes all panels in a row
   share a DV. Here the DV varies **by panel** (d′ vs. threshold), so hiding
   one panel's axis would make it unreadable.
2. **Each panel can keep its own x-label.** The parent labels only the bottom
   row, which assumes a shared grouping variable — true of the current stacked
   layout (`xaxis_bottom_row_only=True`), but *not* of the exploratory schemes,
   where the grouping variable could vary by row, and
   every row needs its own label. Both are supported; the flag picks one.
3. **X-tick labels rotate and the figure widens with group count.** The parent's
   horizontal, space-wrapped labels collide well before 17 levels.
   `_wrap_label()` breaks on spaces **and after `/`**, because these names come
   from free text: plain `textwrap` produces `Other/Uncategori` + `zed`.
4. **Effect sizes and all post-hoc tests always reach the CSV** (the parent
   writes pairwise rows only when it draws them).
5. **It writes `.png`, `.tiff` and `.svg`** (the parent writes `.png` + `.svg`).
   TIFFs use `compression='tiff_lzw'` — lossless, and the difference between a
   3 MB file and a 120 MB one for the scheme-1 figure.

> If you ever want this reusable elsewhere, promote it to
> `modules/visualization_helpers_parts/`.

### Reads

```
data/final/df_public_*.csv     (most recent, via most_recent_public_df)
```

Columns used: `monitor_check`, `monitor_check_operationalized_final`,
`d_prime_overall`, `vch_threshold`.

The script raises on a missing column, a non-numeric DV, or the two DVs no
longer sharing the same missing rows — no silent adaptation.

### Sample

228 rows in the analysis dataframe. Two exclusions, in order:

1. **Non-responders** on the relevant hardware field (see caveats), applied per
   grouping column so answering one question but not the other still counts
   toward the question that *was* answered.
2. **Listwise deletion on the DV** — the 34 VCH-excluded participants. Both DVs
   share exactly the same missing rows (asserted at runtime).

The two exclusions overlap almost completely, so every panel ends at
**N = 193**. **No SP-user filter** — the reviewer's question is about the
measurement device, which is independent of drug exposure.

### Data caveats — read before interpreting any panel

1. **Non-response is not a hardware category.** 15 people left the hardware
   question blank and are dropped (`DROP_NON_RESPONDERS = True`). The exclusion
   is keyed on the **raw free-text field being blank**, not on the
   `"No Response"` label, so it does not depend on how non-response was coded.
   The script asserts afterwards that no `"No Response"` label survives the
   raw-field mask; if the two disagree it raises rather than proceeding.

2. **Non-response is almost perfectly confounded with VCH exclusion.** 14 of
   the 15 non-responders already have missing d′/threshold, so only **one** was
   ever in the analysis sample. Dropping moves per-panel N from 194 to 193 and
   removes the `"No Response"` level (n = 1) entirely. People who skipped the
   hardware question are essentially the people who did not produce usable task
   data — consistent with partial-completion dropout, not with anything about
   hardware. Set `DROP_NON_RESPONDERS = False` to test the coding as it stands.

3. **Single-participant sensitivity.** With one non-responder in the analysis
   sample, a result at the margin can turn on that participant. Per-panel N and
   per-group n are in the summary CSV (`sample_size`, `group_ns`) and printed
   above every box, so a caption can state the cell sizes honestly.

### Outputs

```
results/supplement/hardware_keydown_check/
  supplementary_figure_s8.png     manuscript drafting / Docs preview
  supplementary_figure_s8.tiff    journal submission, LZW, same DPI
  supplementary_figure_s8.svg     true-vector figure assembly
  summary_results/
    supplementary_figure_s8.csv        this scheme
    hardware_operationalization_all_schemes.csv  every active panel, stacked
                                                 (112 rows at current scope)
```

Summary CSV columns: `scheme`, `panel`, `dv`, `group_var`, `test`, `n_groups`,
`sample_size`, `groups`, `group_ns`, `row_type` (`overall` |
`pairwise_posthoc`), `comparison`, `df`, `statistic`, `p_value`,
`p_value_bonferroni`, `effect_size_name`, `effect_size`, and the
`group_{1,2}_{n,median,iqr_q1,iqr_q3}` block (populated for Mann-Whitney rows).

### How to run

```bash
cd hppd_manuscript_public
/usr/local/bin/python3.12 04_visualizations/supplement/hardware_keydown_check.py
```

Runs in a few seconds. Prints the per-panel test results and an omnibus table
suitable for pasting into a reviewer response.

### Common things to change

- `DROP_NON_RESPONDERS` (default `True`) — exclude vs. retain the people who
  left a hardware question blank (caveat 2). Keyed on the raw free-text field,
  not the `"No Response"` label; `RAW_FIELD_FOR` maps each family of
  operationalisation columns to its backing raw field.
- `PAIRWISE_LINE_MAX_GROUPS` (5) — how many groups a panel may have before
  post-hoc brackets are suppressed from the figure. The tests reach the CSV
  regardless.
- `ACTIVE_SCHEME` / `ACTIVE_STEM_INDEX` — which grouping column is rendered
  (`'final'`, monitor). The journal sizing assumes 3 groups in a single stacked
  column, so changing this also means revisiting **Figure format** above.
- `show_group_n` / `group_n_fontsize` (args to `hardware_boxplot_grid`) — the
  per-group "n = …" row above the boxes. On by default.
- `GROUP_LABELS` — the x-axis title. The rendered entry reads **"Display
  Hardware Reported"**, a reader-facing description of the variable. What the
  three levels mean is documented under **Purpose** above.
- `DV_COLS` / `GROUP_STEM` — the column and row variables. `DV_COLS` order is
  the column order; `GROUP_STEM` order is the row order.
- Panel geometry: `FIG_WIDTH_IN` (3.602″ = half a 183 mm page) and
  `PANEL_HEIGHT_IN` (2.6″ per row) set the printed size. `width_per_group`
  (0.85), `min_panel_width` (7.0) and `max_panel_width` (13.0) apply only when
  `fig_width=None`; the cap exists because 17 levels at an uncapped width
  produce a 32-inch-wide figure no supplement page can use.
- Font sizes: the `FONT_*` constants (see **Figure format** above). Change them
  there, not at the call site, so the whole figure stays on one scale.
- `panel_labels` (passed from `main()` as `['a', 'b']`, generated from
  `DV_COLS`) and `panel_label_fontsize`. Reordering `DV_COLS` reorders both the
  rows and the letters.
- Level order defaults to **descending group n** (a display choice only — both
  tests are invariant to group order). Pass `'order'` in a panel spec to fix it.

---

## `diagnostics/` — Diagnostic PDF compilation

This subdirectory contains scripts that compile brms model diagnostic figures
into multi-page PDFs for supplement submission. All scripts resolve paths
relative to their own location via `_PROJECT_ROOT`, so they can be run from
any working directory.

`0X_all_figures.py` adds `supplement/diagnostics/` to `sys.path` so that
the existing `from compile_mediation_diagnostic_pdfs import ...` and
`from compile_single_path_diagnostic_pdfs import ...` imports work.

### Staleness detection

Both pipelines use **mtime-based staleness checks** so that re-syncing
updated HPC results automatically triggers diagnostic regeneration:

- **Mediation:** `generate_compiled_figures()` compares the `fit_*.RData`
  mtime against the oldest `compiled_*.png`. If the fit is newer, R is re-run
  even if all compiled PNGs exist.
- **Single-path:** `nonsp_diagnostic_worker.R` compares the compilation PNG
  mtime against the fit `.RData` and all source PNGs (pp_check, traces,
  heteroscedasticity). If any input is newer, the compilation is regenerated.

To force full regeneration regardless of timestamps, pass `--force` (mediation)
or delete the output PNGs (single-path, on HPC).

---

### Pipeline: regenerating mediation diagnostic PDFs

Produces `results/supplement/diagnostics/mediation_diagnostics_{dv}.pdf` for
`dv` in `[caps_vision, hppd_binary]`.

**Configurable values** (change these to target different models):
- `model_type` — covariate-set suffix. Default: `nice_covariates_spusers`.
  Passed via `--model-type` CLI flag or the `model_type=` keyword argument.
- `DVS` — DV list. Hardcoded in `compile_mediation_diagnostic_pdfs.py`
  (`DV_DIRS = ['caps_vision', 'hppd_binary']`).

```bash
# ── Step 0: Open SSH socket (one-time, lasts until you close it) ─────────
ssh -MNf bouchet
# Approve the DUO push on your phone.

# ── Step 1: Pull latest mediation results from HPC ──────────────────────
#   Pulls CSVs, DHARMa PNGs, PP-check PNGs, and .RData fit objects.
#   Results land in results/{dv}/mediation_models/{model_name}/.
#   To change which models are pulled, edit CONFIG in 03_hpc/generate_hpc_jobs.py.
cd 03_hpc
python compile_mediation_results.py subset

# ── Step 2: Generate compiled diagnostic PNGs + assemble PDFs ────────────
#   Scans results/{dv}/mediation_models/ for directories ending in the
#   model_type, calls _compile_diagnostics_helper.R per model (parallel),
#   then assembles multi-page PDFs.
#
#   The staleness check auto-regenerates any model whose fit_*.RData is
#   newer than its compiled_*.png files (i.e. after Step 1 syncs new fits).
#
#   CONFIG: --model-type controls which covariate-set directories to include.
#           --workers sets R parallelism (default 4).
cd ../04_visualizations/supplement/diagnostics
python compile_mediation_diagnostic_pdfs.py --model-type nice_covariates_spusers

# Variants:
#   python compile_mediation_diagnostic_pdfs.py --workers 8   # faster
#   python compile_mediation_diagnostic_pdfs.py --force        # re-run R even if cached
#   python compile_mediation_diagnostic_pdfs.py --pdf-only     # skip R; PDFs from cached PNGs
```

**Output:**
```
results/supplement/diagnostics/mediation_diagnostics_caps_vision.pdf
results/supplement/diagnostics/mediation_diagnostics_hppd_binary.pdf
```

---

### Pipeline: regenerating single-path diagnostic PDFs

Produces `results/supplement/diagnostics/single_path_diagnostics_{dv}.pdf` for
`dv` in `[caps_vision, hppd_binary]`.

**Configurable values** (change these to target different models):
- `MODEL_TYPE` — covariate-set suffix. Default: `nice_covariates_spusers`.
  Hardcoded at top of `compile_single_path_diagnostic_pdfs.py` (line ~67).
- `DVS` — DV list. Hardcoded at top of `compile_single_path_diagnostic_pdfs.py`
  (line ~62).
- `DEFAULT_MODEL_TYPE` — same value, in `03_hpc/generate_nonsp_diagnostic_jobs.py`
  (line ~67). Must match `MODEL_TYPE` above.
- `IVTYPES` — IV type groups. Hardcoded in `compile_single_path_diagnostic_pdfs.py`
  (line ~59).

```bash
# ── Step 0: Open SSH socket (one-time) ───────────────────────────────────
ssh -MNf bouchet

# ── Step 1: Generate HPC diagnostic jobs (if not already done) ───────────
#   Reads combined_all_analyses.txt (from generate_hpc_jobs.py) and writes
#   nonsp_diagnostic_jobs.txt — one dSQ job per (predictor, model_type, dv).
#
#   CONFIG: --model-type controls which covariate-set models get diagnostic jobs.
cd 03_hpc
python generate_nonsp_diagnostic_jobs.py --model-type nice_covariates_spusers

# ── Step 2: Transfer worker script + job file to HPC and submit ──────────
#   The script prints exact scp + dSQ commands. Copy-paste them:
scp nonsp_diagnostic_worker.R transfer-bouchet:/nfs/roberts/scratch/pi_arp29/msg74/aim1_baseline_final/nonsp_predictor_analyses/
scp nonsp_diagnostic_jobs.txt  transfer-bouchet:/nfs/roberts/scratch/pi_arp29/msg74/aim1_baseline_final/nonsp_predictor_analyses/
ssh bouchet
  cd /nfs/roberts/scratch/pi_arp29/msg74/aim1_baseline_final/nonsp_predictor_analyses
  module load dSQ
  dsq --job-file nonsp_diagnostic_jobs.txt --mem-per-cpu 8g -t 15:00 --mail-type ALL
  # → run the sbatch command that dSQ prints
  exit

# ── Step 3: Wait for HPC jobs to finish (check with squeue -u msg74) ────

# ── Step 4: Retrieve compilation PNGs from HPC ──────────────────────────
#   Option A (programmatic — recommended):
python generate_nonsp_diagnostic_jobs.py --retrieve

#   Option B (manual tarball — the script prints exact commands if the
#   hpc_mirror directory does not exist):
cd ../04_visualizations/supplement/diagnostics
python compile_single_path_diagnostic_pdfs.py
#   ↑ If hpc_mirror/ is missing, this prints tarball retrieval commands.

# ── Step 5: Assemble per-DV PDFs ────────────────────────────────────────
#   Scans data/final/nonsp_predictor_analyses/hpc_mirror/ for
#   {predictor}/{model_type}/results/diagnostics/{dv}_diagnostic_compilation.png
#   and assembles multi-page PDFs.
cd ../04_visualizations/supplement/diagnostics
python compile_single_path_diagnostic_pdfs.py
```

**Output:**
```
results/supplement/diagnostics/single_path_diagnostics_caps_vision.pdf
results/supplement/diagnostics/single_path_diagnostics_hppd_binary.pdf
```

> **Note on the single-path staleness check:** The mtime check runs on HPC
> (inside `nonsp_diagnostic_worker.R`), not locally. If you re-run
> `nonsp_predictors.R` for a model and its fit/PNGs are now newer than the
> existing compilation PNG, the diagnostic worker will automatically regenerate
> it on the next dSQ submission. Locally, the staleness boundary is the
> `--retrieve` step — re-pulling always overwrites the local mirror.

---

### `diagnostics/create_mediation_diagnostic_compilation.py` + `diagnostics/_compile_diagnostics_helper.R`

Generates a supplementary diagnostic compilation figure (`diagnostic_compilation.png`)
for a single brms mediation model results directory. All diagnostic figures are regenerated
directly from the `.RData` fit object (not assembled from pipeline PNGs), ensuring white
backgrounds and properly labelled DHARMa axes.

**Layout** — dynamic based on DV, saved to `{model_dir}/diagnostic_compilation.png`:

*caps_vision (5 rows x 2-3 columns):*

| Row | Columns | Contents |
|-----|---------|---------|
| 0 | 2 | Posterior predictive checks (DV left, mediator right) |
| 1 | 3 | MCMC traces, **mu submodel** — spvar->DV . med->DV . spvar->med |
| 2 | 3 | MCMC traces, **hu submodel** — spvar->DV . med->DV . spvar->med (col 3 blank) |
| 3 | 2 | DHARMa comprehensive — QQ + residuals vs. fitted (DV left, mediator right) |
| 4 | 3 | DHARMa residuals vs. spvar (DV) . vs. mediator (DV) . vs. spvar (mediator) |

*hppd_binary and all other DVs (4 rows x 2-3 columns):*

| Row | Columns | Contents |
|-----|---------|---------|
| 0 | 2 | Posterior predictive checks |
| 1 | 3 | MCMC traces, mu submodel — spvar->DV . med->DV . spvar->med |
| 2 | 2 | DHARMa comprehensive |
| 3 | 3 | DHARMa residuals vs. predictors |

The mu trace row is generated at half height (H3/2 = 2.05") for caps_vision and full
height (H3 = 4.1") for other DVs, so aspect ratios are preserved in both layouts.
The hu trace row is always H3/2.  For DVs with no hu submodel (e.g. hppd_binary),
blank white placeholder PNGs are saved for the hu files to satisfy the cache check.

Column headers above row 0 use DAG notation:
`<spvar_label>  ->  <dv_label>` and `<mediator_label>  ->  <dv_label>`.

**Workflow:**
1. Python (`create_mediation_diagnostic_compilation.py`) parses the model directory name to
   extract DV / spvar / mediator, resolves labels from `VARIABLE_REGISTRY`, and calls
   `_compile_diagnostics_helper.R` via `subprocess`.
2. R script generates 13 `compiled_*.png` files in the model directory and exits.
3. Python assembles the PNGs into the final figure with `matplotlib.gridspec`.

Generated `compiled_*.png` files are cached — subsequent calls skip the R step unless
the fit `.RData` is newer than the oldest compiled PNG (mtime staleness check), or
`--force` / `force_regenerate=True` is passed.  **Important:** existing models have the
old `compiled_traces_spvar_dv.png` / `compiled_traces_med_dv.png` / `compiled_traces_spvar_med.png`
files on disk.  These are no longer used; the new `compiled_traces_mu_*` and
`compiled_traces_hu_*` files are required.  Re-run with `--force` to regenerate.

**CLI:**
```bash
cd 04_visualizations/supplement/diagnostics
python create_mediation_diagnostic_compilation.py \
    ../../../results/caps_vision/mediation_models/caps_vision_avgdose_vchrate_nice_covariates_spusers
# Add --force to re-run R even if compiled_ files already exist.
```

**API (from `0X_all_figures.py` or other scripts):**
```python
from create_mediation_diagnostic_compilation import make_diagnostic_compilation
make_diagnostic_compilation(
    '../results/caps_vision/mediation_models/caps_vision_avgdose_vchrate_nice_covariates_spusers',
    force_regenerate=False,   # True to re-run R helper
    dpi=150,
)
```

**Integration in `0X_all_figures.py`:** Set `RUN_DIAGNOSTIC_COMPILATION = True` (inside
`RUN_DIAGNOSTICS` block) to batch-generate compilation figures for all four
manuscript mediation panels (Figs 2e, 2f, 4c, 4d). Defaults to `False` because the R
step takes several minutes per model.

**Known spvar / mediator shorthands** (parsed from the directory name):

| Shorthand | Full column name |
|-----------|-----------------|
| `spage` | `psychedelic_age` |
| `avgdose` | `avg_life_dose` |
| `lifenomic` | `psycheduse_life_nomic` |
| `vchrate` | `vch_bl_yes_0` |
| `vchthreshold` | `vch_threshold` |
| `vchnu` | `vch_nu` |
| `vchbeta` | `vch_beta` |

Add new shorthands to `SPVAR_TO_FULL` / `MEDIATOR_TO_FULL` at the top of
`create_mediation_diagnostic_compilation.py` when adding new mediation models.

**`_compile_diagnostics_helper.R`** — companion R script called as a subprocess.
Takes 15 positional CLI arguments (fit path, output dir, variable names and labels,
figure dimensions; arg 15 is `H3_HALF`). Saves 13 `compiled_*.png` files:
`compiled_pp_check_{dv,med}.png`, `compiled_dharma_{dv,med}.png`,
`compiled_traces_mu_{spvar_dv,med_dv,spvar_med}.png`,
`compiled_traces_hu_{spvar_dv,med_dv,spvar_med}.png`,
`compiled_resid_dv_sp.png`, `compiled_resid_dv_med.png`, `compiled_resid_med_sp.png`.
Do not call this script directly — use `create_mediation_diagnostic_compilation.py`.

---

### `diagnostics/compile_mediation_diagnostic_pdfs.py`

Batch generates `diagnostic_compilation.png` for every parseable brms mediation
model directory (those with a `fit_*.RData` file), then assembles per-DV multi-page PDFs.

**Output:**
```
results/supplement/diagnostics/mediation_diagnostics_caps_vision.pdf
results/supplement/diagnostics/mediation_diagnostics_hppd_binary.pdf
```

Each PDF contains one page per model, ordered by `(spvar, mediator, model_name)`.
Separator pages are inserted whenever the `(spvar, mediator)` group changes, so
all covariate-set variants for a given spvar x mediator pair appear together.

**Usage:**
```bash
cd 04_visualizations/supplement/diagnostics
python compile_mediation_diagnostic_pdfs.py              # 4 parallel R workers
python compile_mediation_diagnostic_pdfs.py --workers 8  # faster
python compile_mediation_diagnostic_pdfs.py --force      # re-run R even if cached
python compile_mediation_diagnostic_pdfs.py --pdf-only   # skip R; assemble from
                                                          # existing compiled_*.png only
```

**Runtime:** R generation (the bottleneck) runs in parallel. With 6 workers and ~280 models,
expect ~2 hours total. Already-cached models (all `compiled_*.png` present) are skipped.
After first run, subsequent runs are fast — `--pdf-only` assembles PDFs from cached PNGs in
a few minutes.

**Models skipped:** directories without a `fit_*.RData` file, or whose name uses spvar/mediator
shorthands not in `SPVAR_TO_FULL` / `MEDIATOR_TO_FULL` in `create_mediation_diagnostic_compilation.py`.
A count of skipped models is printed at the end.

---

### `diagnostics/compile_single_path_diagnostic_pdfs.py`

Assembles per-DV diagnostic PDFs for the nonsp single-path brms models
whose forest plots appear in the manuscript figures.

Reads 2x2 diagnostic compilation PNGs produced on HPC by
`nonsp_diagnostic_worker.R` and assembles them into two multi-page PDFs
(one per DV), organized by IV type column:
  col 1 — sp_predictors
  col 2 — vch_behavior
  col 3 — vch_computations

**Output:**
```
results/supplement/diagnostics/single_path_diagnostics_hppd_binary.pdf
results/supplement/diagnostics/single_path_diagnostics_caps_vision.pdf
```

> **Status: NOT YET TESTED end-to-end.** The PDF stage was rewritten to include both the
> 6 original DHARMa panel PNGs (`results/diagnostics/`) and the pp_check /
> convergence PNGs added in June 2026. The local `hpc_mirror` only contains
> `pp_checks/` and `convergence_tests/` from a prior partial pull; the
> `results/diagnostics/` PNGs have never been pulled.

**Prerequisites:** Diagnostic compilation PNGs must be available locally in
`data/final/nonsp_predictor_analyses/hpc_mirror/`. After running
`nonsp_diagnostic_worker.R` jobs on HPC, retrieve with `--retrieve` or rsync
(the script prints exact tarball commands if the mirror directory is missing).
See the "Pipeline: regenerating single-path diagnostic PDFs" section above for
the full step-by-step workflow.

**HPC staleness check:** `nonsp_diagnostic_worker.R` compares the compilation
PNG mtime against its inputs (fit `.RData`, pp_check, traces, heteroscedasticity
PNGs). If any input is newer, the compilation is regenerated automatically on
the next dSQ submission — no manual deletion required.

**Usage:**
```bash
cd 04_visualizations/supplement/diagnostics
python compile_single_path_diagnostic_pdfs.py
```

**CONFIG** (top of file) mirrors the relevant config from `0X_all_figures.py`.
When changing the manuscript model type or IV groups, update both files.
