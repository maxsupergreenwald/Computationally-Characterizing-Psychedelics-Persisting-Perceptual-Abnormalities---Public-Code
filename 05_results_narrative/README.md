# Results Narrative — `05_results_narrative/`

One script. It recomputes every statistic reported in the manuscript's Results
section from the shipped analysis dataframe and the model-result CSVs, drops
them into the manuscript prose, and writes the finished section to a plain
`.txt` file.

The output is the Results section **verbatim** — no banner, no separators, no
markup. What the file contains is what the submitted manuscript says.

---

## Quick Start

```bash
cd hppd_manuscript_public/05_results_narrative
python results_narrative.py
```

Writes `results_narrative_output_editted.txt` (the `OUTPUT_TXT` constant) next
to the script, and prints the same text to stdout along with a running log of
every intermediate quantity.

The script anchors all in-repo paths to its own location, so it runs correctly
from any working directory.

### Files

| File | Description |
|---|---|
| `results_narrative.py` | The whole pipeline: loads data, computes statistics, assembles prose, writes the `.txt` |
| `results_narrative_output_editted.txt` | Generated output — the manuscript Results section |
| `README.md` | This file |

---

## `results_narrative.py`

### Structure

| Section | Contents |
|---------|----------|
| **1 — Imports & Setup** | Imports, `_REPO_ROOT` anchoring, `master_config` imports |
| **2 — Data Loading & Config** | Loads the analysis dataframe, the recruitment export and the counterfactual CSV; defines every figure/table label, the citation numerals and the output filename |
| **3 — Utility Functions** | Formatting helpers, row selection, and the result-string builders |
| **4 — Variable Calculations** | Every named quantity that appears in the prose, blocked by manuscript section |
| **5 — Narrative Text Blocks** | The prose f-strings, one per manuscript section |
| **6 — Output** | `build_results_text()` / `write_results_txt()` |

### Data files loaded

| What | Where from |
|---|---|
| Analysis dataframe | `most_recent_public_df(data/final)` → newest `df_public_*.csv`; prints which file it chose. Already fully prepared — every derived column is in the CSV, so there is no preparation step to run |
| SP-user subsample | `df[df[SP_USER_COL] == SP_USER_VALUE]` — the canonical filter, from `master_config` |
| Recruitment export | `RECRUIT_CSV` in `modules/master_config.py` → `data/final/df_recruit_public_<date>.csv`, the de-identified screening export |
| Trial-level belief states | `data/final/vch_master_public.csv` — named directly; raises `FileNotFoundError` if absent |
| Single-path model results | `results/sensitivity_analyses_single_paths/existingresults_manuscript_counterfactual.csv` (`COUNTERFACTUAL_CSV`) |
| Mediation A/B/C′ paths | `results/{dv}/mediation_models/{model}/path_counterfactual_summary.csv` |
| Mediation indirect effects | `results/{dv}/mediation_models/{model}/mc_mediation_summary.csv` |
| Nonparametric group tests | `results/{dv}/{iv_type}/data_visualization/summary_results/boxplot_grid.csv` |
| Nonparametric correlations | `results/{dv}/{iv_type}/data_visualization/summary_results/correlation_grid.csv` |

**Nothing is silently substituted.** A missing counterfactual row raises a
`ValueError` that names the dv × spvar and lists the covariate sets that *are*
present. A missing narrative block raises from `build_results_text()`. The one
deliberate exception is `_cf_result_string_or_fallback` (see below), which falls
back loudly.

### Effect scale

Every reported effect is a **response-scale marginal contrast**,
Δ = E[Y | X = μ + 1 SD] − E[Y | X = μ], with the reference at 0. There is no
log-scale / OR / IRR mode: the counterfactual scale is the only one the script
reports.

- **Single paths** come from `COUNTERFACTUAL_CSV`.
- **Mediation A/B/C′ paths** come from `path_counterfactual_summary.csv`.
- **Indirect effects (NIE)** come from `mc_mediation_summary.csv`.
- **CAPS analyses use the SP-user sample** (`nice_covariates_spusers`, n ≈ 130
  with VCH data).

**Format:** `Δ = X, P(Δ>0) = Y%, 94% HDI [A, B]` — the sign gives the direction,
and `P` is whichever tail lies in the direction of the estimate.

For `hppd_binary` outcomes the contrast is a probability, so `_cf_ppa_pct_parts`
scales Δ and the HDI by ×100 and appends `%`. Its inline "delta only" string
trims a trailing `.0` (8.0 → `8%`) to match the manuscript; the full string and
the HDI keep one decimal.

**Point estimates never read a column name directly.** Every Δ goes through
`master_config.point_estimate()`, which reads `POINT_ESTIMATE_COL` and raises
`PointEstimateColumnMissing` — naming the file and the refit required — if that
column is absent from the CSV.

### `Δ` vs `Δmed`

Indirect effects out of `mc_mediation_summary.csv` are the posterior **median**
(`MC_EFFECT_POINT_ESTIMATE_COL`), not the mean: the Monte-Carlo integration over
the mediator's posterior predictive has tails heavy enough that a handful of
draws can supply the whole mean. Every other Δ in the file is a posterior mean.
Rationale in the `Reported point estimates` block of `modules/master_config.py`.

`_cf_nie_string()` therefore labels its output `Δmed = …`. The HDI and the
direction probability are quantile-based and unaffected, so they keep the plain
`Δ`.

**All eight indirect effects reported in the manuscript carry the `Δmed` label** —
the four in the VCH-behaviour section (threshold and VCH rate, for both outcomes)
and the four in the HGF section (β and ν). Four come straight from
`_cf_nie_string()`; the other four are assembled by hand in Section 4
(`MHA_allstats`, `MHV_allstats`, `MDT_allstats`, `MDV_allstats`) because the
sentence needs the magnitude and the interval in separate places. If you add a
ninth, label it the same way — a bare `Δ` on an indirect effect would read as a
posterior mean.

⚠️ The manuscript-sync tool in the master repository
(by a tool in the private master repo; not shipped here)
matches these strings with a regex. Its `STAT_PAT` must keep the `Δmed = `
alternative **ahead of** the bare `Δ = ` one — alternation is first-match, so the
other order makes every `Δmed` statistic silently stop syncing into the
manuscript.

### Known lookup quirks

- **HPPD outcome key.** The counterfactual CSV stores the HPPD outcome under
  `dv = 'hppd_binary'`. Older CSV versions used `'persist_vis_yn'`; a lookup with
  that string returns zero rows and raises.
- **`vch_threshold` × `caps_vision` × `nice_covariates_spusers`** was not run at
  the counterfactual scale. `_cf_result_string_or_fallback` falls back to
  `nice_covariates` (full sample) with a printed warning and a
  `[full-sample fallback — … not available]` note appended to the string. As of
  the current result set the fallback does **not** fire; if it ever does, the
  note appears in the output text and the run log.
- **`nice_covariates_nopsychosis_iqr`** does not exist in the counterfactual
  CSV. No sentence depends on it; the ν sensitivity analysis reported in the
  text uses `nice_covariates_spusers_iqr` and
  `nice_covariates_beta_spusers_iqr`, both of which are present.
- **hu (hurdle) components are not reported.** They live in
  `hu_paths_summary.csv` inside each mediation model directory. All reported
  effects are mu (count/continuous) components.

---

## Output format

`build_results_text()` assembles the file as:

```
<recruitment paragraph>          ← the lead paragraph; no heading
                                  ← blank line
<section heading>
<paragraph>
                                  ← blank line
<paragraph>
                                  ← blank line
<section heading>
...
```

Three rules make this exact, and all three are easy to break:

1. **Blocks are emitted without stripping.** Each narrative literal starts and
   ends with a newline; those newlines *are* the blank lines in the output.
   `beta_sdt_results` opens with a blank line because the manuscript does.
2. **Trailing spaces at the end of some paragraphs are intentional** — they come
   from the manuscript. An editor configured to trim trailing whitespace on save
   will silently change the output.
3. **A missing or empty block raises** rather than being skipped. An absent
   section means an upstream result failed to load; dropping it silently would
   hide that.

`RESULT_SECTION_ORDER` is the list of `(heading, variable)` pairs, in manuscript
order. `LEAD_PARAGRAPH_VAR` names the unheaded opening paragraph.

### Narrative blocks

| Variable | Manuscript section |
|---|---|
| `recruitment_results_text` | Participant flow (lead paragraph, no heading) |
| `clinical_demographic_results` | Demographics, clinical, and SP use history |
| `ppa_history_results` | PPA History & Current PPAs |
| `ppa_hx_sp_results` | Earlier age at first SP use → lifetime PPA risk |
| `caps_sp_results` | Higher average SP doses → more current PPAs |
| `vch_behavior_results` | Detection threshold and VCH rate (incl. their mediation) |
| `vch_computations_results` | HGF parameters (incl. their mediation) |
| `beta_sdt_results` | β, SDT and confidence exploratory analyses |

---

## Figure, table and citation numbering

Everything the narrative can reference is a named constant in the
`FIGURES AND TABLES LISTED HERE` block in Section 2. Renumber there and nowhere
else.

**Main figures** derive from `FIRST_FIGURE_NUM` (currently `1`):

| Constant | Value | Contents |
|---|---|---|
| `PPA_FIG` | 2 | PPA history distributions |
| `CAPS_FIG` | 3 | CAPS item distributions |
| `FIG_SP_PREDICTORS` | 4 | SP predictors of PPA history & CAPS vision |
| `FIG_VCH_BEHAVIOR` | 5 | VCH behavioural results + their mediation diagrams |
| `FIG_VCH_COMPUTATIONS` | 6 | HGF parameter results |
| `FIG_MEDIATION` | 6 | Same figure as above — the two were merged |
| `FIG_BETA` | 7 | β / SDT exploratory panel |

**Panel letters** are written into the per-section reference constants at the
bottom of each Section 4 block (e.g. `vch_behavior_ppa_hx_regression`,
`mediation_caps_vchrate_panel`). Figure 5's mediation panels are laid out
row-wise:

| Panel | Predictor → Outcome | Mediator |
|---|---|---|
| `5e` | `psychedelic_age` → PPA history | `vch_threshold` |
| `5f` | `avg_life_dose` → CAPS vision | `vch_threshold` |
| `5g` | `psychedelic_age` → PPA history | `vch_bl_yes_0` (VCH rate) |
| `5h` | `avg_life_dose` → CAPS vision | `vch_bl_yes_0` (VCH rate) |

Figure 6 carries the HGF mediation panels: `6n` (age → β → PPA history) and
`6o` (dose → β → CAPS vision).

**Supplementary references actually reaching the text:**

| Constant | Value |
|---|---|
| `supp_fig_consort` | Supplementary Figure S9 |
| `SP_FIG` / `SP_FIG_SHORT` | Supplementary Figure S6 / Supplementary Fig. S6 |
| `sensitivity_analysis_heatmap_single_path` | Supplementary Figure S4 |
| `task_engagement_fig_threshold_error` | Supplementary Figure S7 |
| `mann_whitney_u_table` | Supplementary Table S2 |
| `supp_table_mediation_sensitivity` | Supplementary Table S3 |

Only references a sentence actually cites are defined; there is no longer a
registry of unreferenced supplementary labels to drift out of date.

**Citations** are superscript numerals keyed to the manuscript reference list,
held in the `CITE` dict:

| Key | Numeral |
|---|---|
| `izmi_2024` | 3 |
| `carhart_nutt_2010` | 6 |
| `baggott_2011` | 7 |
| `zhou_2025` | 8 |
| `muller_2022` | 9 |
| `kvam_2023` | 10 |
| `kessler_2005` | 56 |
| `hirschfeld_2023` | 57 |
| `sdt_low_sensitivity` | 58 |

Grouped citations are assembled in the prose — `,`-joined for a list
(`6,7,10`), en-dash-joined for a contiguous run (`6–9`). Renumbering the
reference list means editing this dict only.

---

## Data flow and integrity

- **`df`** (N = 228) — the shipped analysis sample. Every derived column is
  already present; there is no preparation step in this repo.
- **`df_sp`** (N = 186) — SP users, `psycheduse_yn == "Yes"` via `SP_USER_COL` /
  `SP_USER_VALUE`.
- **`df_hppd`** — PPA-positive subset (`persist_vis_yn > 0`), defined once.

**Neither `df` nor `df_sp` is mutated anywhere after loading.** The QC audit loop
in Section 4 checks `df_recruit` for participants failing 2+
attention/effort/distraction checks and emits a **warning only** — it does not
filter `df`. `RFL = 0` in the current data, so the branch never fires. If it ever
does, the exclusion belongs upstream in the data export, not here.

All reported statistics come from `df` or `df_sp`; `df_recruit` is used only for
the recruitment counts.

### Recruitment counts

Two sets of counts coexist and they mean different things.

**Marginal counts** (`RS`, `RE`, `RSP`, …) are each computed over everyone in
`df_recruit`, so they do not subtract into a chain.

**Stage-conditioned counts** (`RE_INELIGIBLE`, `RE_MINMEASURES`, `RE_AFTERMIN`,
`RE_QCFAIL`, `RE_COMPLETERS`, `RE_SALVAGED`, `RE_RETAINED`) place every
participant at exactly one stage, so the flow closes exactly. `_flow_checks`
asserts every step at runtime and prints the closure; a failure prints a
`WARNING: recruitment flow no longer closes` block naming the broken identity.

**`RE` (eligible = 492)** is `RS - len(df_eligibility[si_2_v2.isna()])`, not
`len(df_recruit[si_2_v2.notna()])`. Early-cohort participants were enrolled
before the `geo_crit` REDCap field existed, so the ineligibility loop's
`elif pd.isna(row["geo_crit"])` branch mislabels ~133 of them as
fraud-associated. Restricting `df_eligibility` to `si_2_v2.isna()` isolates the
757 who are both loop-ineligible *and* never found eligible by REDCap:
`1249 − 757 = 492`.

**Two orthogonal descriptions of the same flow.** The recruitment paragraph and
the CONSORT diagram (`04_visualizations/supplement/consort_diagram.py`) describe
the same participants by different routes. Both close on 228, and the assertions
tie them together.

| | Route | Closes as |
|---|---|---|
| Diagram | subtract exclusions from the eligible pool | `492 − 195 = 297`, `297 − 69 = 228` |
| Paragraph | count completers, remove QC failures, add salvage back | `270 − 69 = 201`, `201 + 27 = 228` |

Every quality-control failure sits inside the completers and none is a salvage
record (asserted), so the two routes never double-subtract anyone.

**The `"Completed"` row (347) of the `df_rec_count` audit table is neither of
these.** It counts participants who reached the end-of-survey QC items
(`honesty_qc` not null) or completed as a student — regardless of eligibility or
whether the record was analysable — and is not a reliable completion flag
(records 2369 and 2438 have `honesty_qc` filled but did not finish the
questionnaire). Treat it as an audit figure, not a flow count; no sentence uses
it.

**`RSTF` (student task-fail records).** `RSTF_MASK` identifies non-SP students
(`student_yn == 1`, `psycheduse_yn == 2`) who failed QC but passed the RAVEN
screen and were not flagged by the RA as bad data (`qc_bad_data`). They failed
the behavioural tasks, were awarded credit, and their data cannot be used — they
are not genuine QC failures. Records 1989, 1997, 2578, 2666, validated at runtime
by the RSTF print. They are excluded from `RFQ` and from the "Failed QC" audit
row, and appear instead among the eligible non-completers.

**Sentence-initial counts are spelled out** by `spell_sentence_initial()`,
derived from the data rather than hardcoded, so the word tracks the number. It
covers 1–99; anything outside that falls through to digits — rephrase the
sentence rather than extending it blindly.

---

## Utility function reference

**Formatting**

| Function | Returns |
|---|---|
| `_fmt_num(x, digits=3)` | Number as string; a value that rounds to an integer prints without decimals |
| `_fmt_trim_num(x, digits=2)` | Same, but strips trailing zeros (`8.0` → `"8"`) |
| `_fmt_p(p)` | `"< 0.001"` or `"= 0.xxx"` |
| `_fmt_prob(prob)` | Posterior probability as `"xx.x%"`, trailing zero trimmed |
| `_num_pct_from_mask(mask, denom, semicolon=None)` | `"n (xx.x%)"`, or `"n; xx.x%"` with `semicolon=True` |
| `_pct_str(df, mask, denom_mask=None)` | Whole-number percentage as a string |
| `_select_single_row(df, mask, context)` | The single matching row; raises naming `context` if none, prints a note and takes the largest-N row if several |

**Statistics**

| Function | Reads | Returns |
|---|---|---|
| `_mw_result_string(nonparam_df, var_name, verbose=None)` | boxplot grid CSV | `"U = …, p …"`, plus rrb and group medians when `MANN_WHITNEY_VERBOSE` |
| `_spearman_string(nonparam_df, row_var, column_var)` | correlation grid CSV | `"ρ = …, p …"` |
| `_partial_spearman_str(df_sub, x, y, covariate='age_v2')` | dataframe | Partial Spearman as `"ρ = …, p …"` — ranks all three variables, then partial Pearson on the ranks (matches `correlation_matrix_plot`) |
| `_perm_mean_diff(plus, minus, n_boot=10000, seed=42)` | arrays | `"Δ = …, 94% CI [.., ..]"` — bootstrap CI on the group-mean difference, matching the trajectory figures |

**Counterfactual result builders**

| Function | Reads | Returns |
|---|---|---|
| `_cf_result_string(cf_df, spvar, cov, dv)` | `COUNTERFACTUAL_CSV` | `(str, n_obs, estimate)`; raises listing available covariate sets if the row is missing |
| `_cf_result_string_or_fallback(cf_df, spvar, preferred, fallback, dv)` | as above | Tries `preferred`, else falls back with a printed warning and an appended note |
| `_cf_ppa_pct_parts(cf_df, spvar, cov, dv)` | as above | `(full_str, delta_abs_str, otherstats_str, n_obs, est_raw)` with Δ and HDI ×100 for probability-scale outcomes |
| `_cf_path_string(cf_path_df, effect_label)` | `path_counterfactual_summary.csv` | `"Δ = …"` for `'A path'`, `'B path'` or `"C' path"` |
| `_cf_nie_string(mc_df)` | `mc_mediation_summary.csv` | `"Δmed = …"` for the `NIE` row |

**Output**

| Function | Description |
|---|---|
| `build_results_text()` | Assembles the full narrative; raises on a missing or empty block |
| `write_results_txt(output_path=OUTPUT_TXT)` | Writes that text to disk |

---

## Toggles

| Constant | Effect |
|---|---|
| `FIRST_FIGURE_NUM` | Shifts every main figure number |
| `OUTPUT_TXT` | Output path (anchored to this directory) |
| `MANN_WHITNEY_VERBOSE` | `True` adds rrb and group medians/IQRs to every Mann-Whitney string. Default `False` — the manuscript reports `U` and `p` only |

Supplementary figure and table numbers are literals, not derived from an offset:
the supplement is not ordered by the same sequence as the reference list, so an
offset constant only obscured which label was which.

**There is no interval-type toggle.** Every reported interval is a 94% HDI, read
from `hdi_lower_94` / `hdi_upper_94` by the `_cf_*` builders. Reporting an
equal-tailed 95% CI instead would mean changing those builders, not flipping a
constant.

---

## CAPS descriptive section — calculation notes

**Denominator:** SP users with a non-missing `caps_bl_1` (`_df_caps`, N ≈ 130).
All percentages in the block use it unless noted.

**CAPS vision items:** `[4, 26, 31, 23, 19, 22]` — binary `caps_bl_{x}`,
frequency `caps_bl_{x}c`, distress `caps_bl_{x}a`.

| Variable | Definition |
|---|---|
| `caps_pos_majority` | `"majority"` / `"minority"`, on whether > 50% of `_df_caps` have `caps_vision > 0` |
| `caps_pos_pct` | % of `_df_caps` with `caps_vision > 0` |
| `caps_ppa_pos_pct` | % of `caps_vision > 0` participants with `persist_vis_yn == 1` |
| `caps_freq_fives` | % of `caps_vision > 0` participants with any `caps_bl_{x}c > 4` ("All the Time") |
| `pct_caps_vision_distress_0` | % where the **max** distress across **endorsed** items is < 2 (all "Not at all") |
| `pct_caps_vision_distress_over3` | % where **any** endorsed item is rated > 3 ("Firmly" or "Very") |
| `chi_square_caps_ppa` | Permutation χ² string — computed and printed for reference; no sentence currently reports it |

**Distress** is evaluated only over items the participant endorsed
(`caps_bl_x == 1`); non-endorsed items (0 = absent) are excluded. Scale:
0 = absent, 1 = Not at all, 2 = Slightly, 3 = Somewhat, 4 = Firmly, 5 = Very.

**Permutation chi-square.** Tests whether the distribution of endorsements
across the six CAPS items differs between PPA-history+ and PPA-history−
participants, restricted to `caps_vision > 0`. Classical χ² is unreliable here
(~9 PPA− participants have current symptoms, so expected cell counts < 5). The
permutation approach shuffles `persist_vis_yn` 10,000 times at the person level,
preserving within-person item correlations, and takes `(B+1)/(N_perm+1)`
(Phipson & Smyth 2010). Seed 42. Current value: χ²(5) = 2.52,
p_permutation = 0.509.

---

## Style guidance for editing the prose

- Mimic the author's naming style and tone; do not impose a template.
- Statistics appear inline immediately after the thing they describe:
  `lifetime SP use ({CB_lifetime}) and average SP dose ({CB_dose})`.
- Variables hold statistics only — no interpretation in the variable name.
- Greek symbols are written directly: `ρ =`, not `rho =`.
- Order from strongest effect to weakest where the argument allows.
- Keep sensitivity text short: if outlier exclusion is the only change, one
  clause is enough.
- **Never hardcode a number into the prose.** Every count, percentage and
  interval in the text comes from a variable computed in Section 4. If a
  sentence needs a number that does not exist yet, compute it there.

**Effect-size reference for Mann–Whitney rank-biserial correlations:**

| rrb | Interpretation |
|---|---|
| ≈ 0.10 | Small |
| ≈ 0.30 | Medium |
| ≥ 0.50 | Large |

If a variable's meaning is unclear, check `../modules/master_config.py` and the
repository `README.md` before writing prose.
