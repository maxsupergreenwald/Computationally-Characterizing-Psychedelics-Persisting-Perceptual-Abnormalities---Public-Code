# `04_visualizations/`

Everything the manuscript shows, except the supplement's standalone figures and
tables (those live in `supplement/`, which has its own README).

| | |
|---|---|
| `0X_all_figures.py` | Figures 2, 4, 5, 6, 7, Table 1, Supplementary Figure S6 |
| `caps_item_distributions_hppd_split.py` | Figure 3 |
| `supplement/` | Supplementary Figures S1–S5, S7–S9 and Tables S1–S6 — see `supplement/README.md` |
| `linked_figures/` | Copies of the assembled figures, downsampled to <1 MB, tracked in git so a Google Apps Script can sync them into the manuscript document |

Figure 1 is a hand-drawn schematic with no generating script. It ships as
`results/final_figures/figure_1.svg`.

---

## Running it

```bash
python 04_visualizations/0X_all_figures.py
```

Run it from anywhere — the script `chdir`s to its own directory at startup, so
all relative paths resolve regardless of where you invoke it.

It reads the compiled model output under `results/` and the shipped analysis
dataframe, and writes both the figures and **the summary CSVs that
`supplement/` and `05_results_narrative/` read**. Run it before either of those.

`caps_item_distributions_hppd_split.py` is independent and can run in any order.

---

## The data it reads

The analysis dataframe ships fully prepared. Read it directly:

```python
from data_prep import most_recent_public_df
DF_PUBLIC = most_recent_public_df(DATA_DIR)      # newest data/final/df_public_*.csv
df = pd.read_csv(DF_PUBLIC, low_memory=False)
```

There is no preparation step to run and no `load_and_prepare_data()` in this
repository — every derived column (`hppd_binary`, `subtle`, `baggot_total`, the
SDT/metacognition block, the recalculated VCH hit rates, corrected
`avg_life_dose`) is already a column of the CSV. The script prints which file it
picked.

**`df_sp` is the SP-user subsample.** The filter is always
`psycheduse_yn == "Yes"` — not `!= "No"`, and not `psycheduse_life_nomic > 0`.
Those give different Ns.

---

## Effect scale and analysis sample

Reported effects are **counterfactual marginal contrasts on the response
scale**: E[Y | X = mean + 1 SD] − E[Y | X = mean], from `posterior_epred()` on
the fitted brms model.

| DV family | What the contrast is |
|---|---|
| `bernoulli` (`hppd_binary`) | risk difference, probability scale |
| `hurdle_negbinom` / `negbinom` (`caps_vision`) | difference in expected count |
| `student_t` / `gaussian` | raw score difference |

These are reported instead of standardized log-scale coefficients because they
are in the units of the outcome. Forest plots therefore read
`results/sensitivity_analyses_single_paths/existingresults_manuscript_counterfactual.csv`
and are drawn by `counterfactual_forest_plot`, with the reference line at **0**
— no `exp()`, no signed-reciprocal axis, no OR/IRR tick remapping.

Mediation and forest-plot models are fitted in **SP users only** (n ≈ 186): the
mediation question is about variation in exposure among the exposed, so SP-naive
participants carry no information about it. Descriptive figures use the full
QC-passing sample (n ≈ 228).

---

## Configuration

All of it is at the top of `0X_all_figures.py`.

| Constant | Default | Purpose |
|---|---|---|
| `MODEL_TYPE` | `"nice_covariates"` | Covariate set whose results the figures read. Also the results directory name. |
| `IVTYPES` | `sp_predictors`, `vch_behavior`, `vch_computations`, `sdt_hppd` | IV groups that get a data-visualisation figure and a forest plot |
| `RESULTS_BASE` | `"../results"` | Root of the compiled model output |
| `FIGURE_FORMAT` / `FIGURE_DPI` | `"png"` / `600` | Intermediate panel format |
| `FINAL_FIGURE_FORMATS` | png, tiff, svg | Journal-submission copies written to `results/final_figures/` |

In Table 1, Fisher's exact test is substituted automatically wherever chi-square
is invalid (any expected cell count < 5) and the contingency table is 2×2;
`statistic_type` records which test produced each p-value in the sidecar CSV.

### Section toggles

Each runs independently; set any to `False` to skip it.

| Toggle | What it does |
|---|---|
| `RUN_HPPD_CAPS_FIGS` | Boxplot grids, correlation grids, forest plots, HGF state trajectories |
| `RUN_TABLES` | Table 1 and the split demographic/clinical tables |
| `RUN_DESCRIPTIVE_FIGS` | SP exposure distribution summaries |
| `RUN_MEDIATION_DIAGRAMS` | Mediation path diagrams |
| `RUN_FIGURE_ASSEMBLY` | Composites the panels into the numbered manuscript figures |
| `RUN_DIAGNOSTICS` | DHARMa review of the mediation models behind Figures 5 and 6 |
| `RUN_DIAGNOSTIC_COMPILATION` | Per-model `diagnostic_compilation.png`. **Shells out to R and is slow** — turn this off if you only want figures. |
| `RUN_SINGLE_PATH_DIAGNOSTICS` | Assembles single-path diagnostic PDFs from the pulled HPC results. No R needed. |
| `RUN_LINKED_FIGURES` | Copies assembled figures into `linked_figures/` |

---

## What produces what

Figures are written to `results/final_figures/` in every format in
`FINAL_FIGURE_FORMATS`; Table 1 goes to `results/descriptive/tables/`.

| | Script | Content |
|---|---|---|
| **Figure 2** | `0X_all_figures.py` | PPA history distributions |
| **Figure 3** | `caps_item_distributions_hppd_split.py` | CAPS item distributions split by PPA history |
| **Figure 4** | `0X_all_figures.py` | SP predictors → PPA / CAPS forest plots |
| **Figure 5** | `0X_all_figures.py` | VCH behaviour + mediation |
| **Figure 6** | `0X_all_figures.py` | VCH computations (HGF) + mediation |
| **Figure 7** | `0X_all_figures.py` | β, signal detection, detection curves |
| **Table 1** | `0X_all_figures.py` | `table_1.docx` |
| **Figure S6** | `0X_all_figures.py` | Routed to `results/supplement/`, not `final_figures/` |

### Figure 3

`caps_item_distributions_hppd_split.py` draws a 2×2 panel over the six binary
CAPS vision items (`caps_bl_{4,26,31,23,19,22}`):

| Panel | Content |
|---|---|
| a | Binary endorsement count and %, PPA(−) vs PPA(+) side by side |
| b | Frequency (`caps_bl_{x}c`) |
| c | Intrusiveness (`caps_bl_{x}b`) |
| d | Distress (`caps_bl_{x}a`) |

Panels b–d are split violins **conditioned on endorsement** — structural zeros
are excluded, so they describe severity among endorsers, not the whole sample.
Rows null on any of the six binary items are dropped first. Sample: SP users
with valid CAPS, n = 130 (PPA(−) = 38, PPA(+) = 92).

---

## The summary-CSV sidecar

When a helper is called with `savepath`, it writes two things:

```
{savepath}.png
{dirname(savepath)}/summary_results/{basename(savepath)}.csv
```

The CSV holds the numbers actually plotted. This is the contract that lets
`05_results_narrative/results_narrative.py` report every statistic without
transcribing it. **If you move a figure, its sidecar moves with it** — the
narrative resolves paths, not figures.

---

## Helper functions

Imported from `modules/visualization_helpers`, which auto-loads
`modules/visualization_helpers_parts/`. See that directory's README for the full
index. Every helper there is called by this pipeline.

Figure labels and palettes come from `modules/master_config.py` —
`dv_to_lab_short`, `binary_palette`, `caps_vision_palette`,
`electric_blue_palette`. Do not define a local label dict or palette.

---

## Pitfalls

1. **Run `0X_all_figures.py` before `supplement/` or `05_results_narrative/`.**
   Both read summary CSVs this script writes. Running them first gives you
   stale numbers or a missing-file error, not a warning.

2. **All figures must use Arial.** Global `rcParams` are set immediately after
   the matplotlib import (search `GLOBAL FONT`). `sns.set_style()` resets
   `font.family`, so any call to it must re-enforce Arial straight afterward.
   The PIL-assembled figures (3 and 6) load `Arial Bold.ttf` explicitly.

3. **`MODEL_TYPE` must exist in `results/`.** It is the results directory name,
   not a label. Pointing it at a covariate set the cluster never ran gives an
   empty figure rather than an error.

4. **`RUN_DIAGNOSTIC_COMPILATION` calls R.** It needs `brms`, `DHARMa` and
   `DHARMa.helpers` installed and the fitted `.RData` objects present under
   `results/`. It is by far the slowest section.

5. **`linked_figures/` is deliberately tracked in git** and deliberately
   downsampled to under 1 MB per image. It exists so an Apps Script can pull
   the current figures into the manuscript document. Do not gitignore it, and
   do not treat its contents as publication-resolution — `results/final_figures/`
   holds those.
