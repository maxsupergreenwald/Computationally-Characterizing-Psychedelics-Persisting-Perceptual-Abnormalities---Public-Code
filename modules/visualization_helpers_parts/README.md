# `visualization_helpers_parts/`

The plotting and table helpers used by `04_visualizations/0X_all_figures.py`.
One helper (or one tightly-related group) per module.

Import through `visualization_helpers`, never from this package directly:

```python
import sys; sys.path.insert(0, '../modules')
from visualization_helpers import multipanel_boxplot_grid, correlation_matrix_plot
```

`../visualization_helpers.py` imports every `.py` file here whose name does not
start with `_`, and re-exports the names in that file's `__all__`. Adding a file
is enough to register it.

Every module in this directory is called by the figure pipeline. If you add one
that isn't, it does not belong here.

---

## File index

| Module | Exports | Produces |
|---|---|---|
| `correlation_matrix_plot.py` | `correlation_matrix_plot` | Scatter grids with per-panel Spearman ρ — the `correlation_grid.png` behind Figures 3–7 |
| `multipanel_boxplot_grid.py` | `multipanel_boxplot_grid` | Box + strip grids split by a binary group — the `boxplot_grid.png` behind Figures 2, 5, 6 |
| `create_mediation_diagram.py` | `create_mediation_diagram` | Mediation path diagrams (Figures 5, 6) |
| `counterfactual_forest_plot.py` | `counterfactual_forest_plot` | Forest plots of response-scale marginal contrasts (Figure 4) |
| `state_trajectories.py` | `plot_state_trajectories`, `compute_state_stats`, `get_state_label`, `DV_GROUPING_MAP` | Block-level HGF state trajectories (Figure 6 d–g) |
| `generate_publication_table.py` | `generate_publication_table`, `generate_publication_table_thickonly`, `generate_combined_split_table_thickonly` | Table 1 (`table_1.docx`) and the split demographic/clinical tables |
| `get_field_label_dict.py` | `get_field_label_dict` | Reads REDCap choice labels from `modules/redcap_data_dictionary.csv` |
| `load_most_recent_csv.py` | `load_most_recent_csv` | Newest file matching a prefix, by mtime |
| `binary_palette.py` | `binary_palette` | Re-export of the canonical two-colour palette from `master_config` |

---

## Conventions

**Labels come from `master_config`.** `create_mediation_diagram` and
`multipanel_boxplot_grid` import `dv_to_lab_short` directly. Do not pass a
different label dictionary unless the figure genuinely needs one.

**Colours come from `master_config`.** `binary_palette`, `caps_vision_palette`
and `electric_blue_palette` are the only palettes used by manuscript figures.

**The `savepath` sidecar.** When a helper is called with `savepath`, it writes
the figure to `{savepath}.png` and the plotted numbers to
`{dirname(savepath)}/summary_results/{basename(savepath)}.csv`. The CSV is what
`05_results_narrative/results_narrative.py` reads — if a figure moves, its
sidecar moves with it.

**Adding a helper.** Define one public function, declare
`__all__ = ["function_name"]` at the bottom, and give the module its own
imports. Do not re-export another module's names; the loader warns on
collisions and the winner is decided by alphabetical order, not by intent.

---

## Pitfalls

1. **`master_config` must be on `sys.path`.** Several modules import from it at
   module level, so `modules/` has to be importable before
   `visualization_helpers` is imported. Every shipped script does
   `sys.path.insert(0, '../modules')` first.

2. **Stale `__pycache__` after a rename.** The loader enumerates `.py` files on
   disk, but Python will happily import a cached module whose source is gone.
   Clear `__pycache__/` after deleting or renaming a helper.

3. **`correlation_matrix_plot` ranks before correlating.** Its ρ is a Spearman
   computed by rank-transforming and then correlating, matching
   `_partial_spearman_str` in `results_narrative.py`. Do not swap in a
   Pearson-on-raw-values implementation; the two would silently disagree.
