"""Counterfactual forest plot: marginal response-scale effects (no exp transform).

This function does NOT apply any exp()
transformation or standardization. Effects from the counterfactual analysis are
already on the response scale:
  - bernoulli (e.g. hppd_binary)        → risk difference (probability scale)
  - hurdle_negbinom / negbinom           → expected count difference
  - student_t / gaussian                 → raw score difference

The reference line is drawn at 0 (not 1). There is no signed-reciprocal axis
transformation, no OR/IRR tick-label remapping, and no subsample standardization.

Required columns in results_df (after renaming counterfactual CSV columns via
the loader in 0X_all_figures.py):
    dv, var, covariates, Estimate, hdi_lower_94, hdi_upper_94

Optional columns: prob_above_0, prob_below_0

Column renaming expected upstream:
    mean     → Estimate   (POINT_ESTIMATE_COL, modules/master_config.py)
    spvar    → var
    cov      → covariates
    settings → model
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def counterfactual_forest_plot(
    results_df: pd.DataFrame,
    dv: str,
    covariates: str,
    x_label: str,
    predictors: list[str] | None = None,
    savepath: str | None = None,
    binary_palette: list[str] | tuple[str, str] = ("#faeee2", "#bc511c"),
    positive_palette_index: int = 1,
    negative_palette_index: int = 0,
    prob_threshold: float = 0.9,
    alpha_high: float = 0.85,
    alpha_low: float = 0.40,
    use_alpha_by_significance: bool = False,
    figsize: tuple[float, float] = (9, 4.5),
    row_spacing: float = 0.6,
    title: str | None = None,
    ylabel: str = "",
    ci_linewidth: float = 9.0,
    ci_linestyle: str = "-",
    estimate_marker: str = "o",
    estimate_size: float = 420,
    vertical_line_alpha: float = 0.5,
    null_line_color: str = "#000000",
    null_line_style: str = "-",
    null_line_width: float = 1.8,
    tick_fontsize: float = 16,
    label_fontsize: float = 18,
    xlabel_fontsize: float | None = None,
    ylabel_fontsize: float | None = None,
    xlabel_labelpad: float = 10,
    xlabel_x_position: float | None = None,
    xlabel_y_position: float | None = None,
    title_fontsize: float = 20,
    y_label_map: dict | None = None,
    predictor_order: list[str] | None = None,
    xlim: tuple[float, float] | None = None,
    seaborn_context: str | None = "poster",
    estimate_edgecolor: str = "none",
    estimate_linewidth: float = 0.0,
    significance_fontsize: float = 28,
    significance_y_offset: float = 0.02,
    savetype: str = "png",
    show_plot: bool = True,
):
    """
    Forest plot for counterfactual (marginal response-scale) effects.

    Parameters
    ----------
    results_df : pd.DataFrame
        Results with columns: dv, var, covariates, Estimate, hdi_lower_94,
        hdi_upper_94. Column names are those produced by renaming the
        counterfactual CSV upstream (median→Estimate, spvar→var, etc.).
    dv : str
        Outcome variable to plot. Rows with results_df["dv"] != dv are dropped.
    covariates : str
        Covariate-set string to plot (e.g. "nice_covariates_spusers"). Rows
        with results_df["covariates"] != covariates are dropped.
    x_label : str
        X-axis label describing the effect scale (e.g. "Change in PPA Probability",
        "Change in Symptom Count"). Comes from cfg["counterfactual_x_label"] in
        DV_CONFIGS.
    predictors : list[str] | None
        If provided, keep only these predictors (after stripping _normalized).
    savepath : str | None
        If provided, figure is saved here (extension added automatically).

    Returns
    -------
    pd.DataFrame
        The filtered/annotated work DataFrame used for plotting.
    """
    req = {"dv", "var", "covariates", "Estimate", "hdi_lower_94", "hdi_upper_94"}
    missing = req - set(results_df.columns)
    if missing:
        raise ValueError(f"Missing required columns in results_df: {sorted(missing)}")

    work = results_df.copy()
    work = work[(work["dv"] == dv) & (work["covariates"] == covariates)].copy()
    work = work[work["var"] != "Intercept"].copy()

    if predictors is not None:
        allowed_canonical = {p.replace("_normalized", "") for p in predictors}
        work = work[work["var"].str.replace("_normalized", "", regex=False).isin(allowed_canonical)].copy()

    if len(work) == 0:
        print(f"\033[1mWARNING:\033[0m No rows found for dv='{dv}' and covariates='{covariates}'.")
        return work

    # ── Probability of direction ─────────────────────────────────────────────
    p_above = pd.to_numeric(work["prob_above_0"], errors="coerce") if "prob_above_0" in work.columns else pd.Series(np.nan, index=work.index)
    p_below = pd.to_numeric(work["prob_below_0"], errors="coerce") if "prob_below_0" in work.columns else pd.Series(np.nan, index=work.index)
    prob = pd.concat([p_above, p_below], axis=1).max(axis=1)
    if prob.dropna().max() > 1.5:  # stored as 0–100
        prob = prob / 100.0
    work["prob_proportion"] = prob

    work["Estimate"]      = pd.to_numeric(work["Estimate"],      errors="coerce")
    work["hdi_lower_94"]  = pd.to_numeric(work["hdi_lower_94"],  errors="coerce")
    work["hdi_upper_94"]  = pd.to_numeric(work["hdi_upper_94"],  errors="coerce")

    # ── Predictor ordering ───────────────────────────────────────────────────
    work["predictor_raw"] = work["var"].str.replace("_normalized", "", regex=False)

    if predictor_order is not None:
        order = [p.replace("_normalized", "") for p in predictor_order]
    elif predictors is not None:
        order = [p.replace("_normalized", "") for p in predictors]
    else:
        order = list(work["predictor_raw"])

    order_map = {k: i for i, k in enumerate(order)}
    work["order"] = work["predictor_raw"].map(order_map)
    fallback = pd.Series(len(order_map) + np.arange(len(work)), index=work.index, dtype=float)
    work["order"] = work["order"].where(work["order"].notna(), fallback)
    work = work.sort_values("order").reset_index(drop=True)

    # ── Labels ───────────────────────────────────────────────────────────────
    if y_label_map is None:
        y_label_map = {}
    work["y_label"] = work["predictor_raw"].map(lambda x: str(y_label_map.get(x, x)))

    # ── Color by sign of estimate (0 = boundary) ─────────────────────────────
    pos_color = binary_palette[positive_palette_index]
    neg_color = binary_palette[negative_palette_index]
    work["color"] = [pos_color if v >= 0 else neg_color for v in work["Estimate"]]

    plot_df = work.dropna(subset=["Estimate", "hdi_lower_94", "hdi_upper_94"]).copy()
    if len(plot_df) == 0:
        print("\033[1mWARNING:\033[0m No rows left after dropping NaN estimates; skipping plot.")
        return work

    # ── Plot ─────────────────────────────────────────────────────────────────
    if seaborn_context:
        sns.set_context(seaborn_context)
        # sns.set_context can mutate global rcParams; re-enforce Arial
        import matplotlib as _mpl
        _mpl.rcParams['font.family'] = 'sans-serif'
        _mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

    y = np.arange(len(plot_df)) * row_spacing
    fig, ax = plt.subplots(figsize=figsize)

    for i, (_, row) in enumerate(plot_df.iterrows()):
        yi = y[i]
        pprob = row["prob_proportion"]
        ci_alpha  = (alpha_high if pd.notna(pprob) and pprob >= prob_threshold else alpha_low
                     ) if use_alpha_by_significance else alpha_high
        ax.hlines(
            y=yi,
            xmin=row["hdi_lower_94"],
            xmax=row["hdi_upper_94"],
            color=row["color"],
            alpha=ci_alpha,
            linewidth=ci_linewidth,
            linestyle=ci_linestyle,
            zorder=2,
        )
        ax.scatter(
            row["Estimate"],
            yi,
            color=row["color"],
            alpha=1.0,
            s=estimate_size,
            marker=estimate_marker,
            edgecolors=estimate_edgecolor,
            linewidths=estimate_linewidth,
            zorder=3,
        )
        star = ""
        if pd.notna(pprob):
            if pprob >= 0.999:
                star = "***"
            elif pprob >= 0.99:
                star = "**"
            elif pprob >= 0.95:
                star = "*"
        if star:
            ax.text(
                row["Estimate"],
                yi - significance_y_offset,
                star,
                ha="center", va="bottom",
                fontsize=significance_fontsize,
                fontweight="bold",
                color="#000000",
                zorder=4,
            )

    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["y_label"], fontsize=tick_fontsize, rotation=0)
    ax.invert_yaxis()

    # Reference line at 0 (not 1 — counterfactual effects are centered on 0)
    ax.axvline(
        0.0,
        color=null_line_color,
        linestyle=null_line_style,
        linewidth=null_line_width,
        alpha=vertical_line_alpha,
        zorder=1,
    )

    if xlim is not None:
        ax.set_xlim(*xlim)

    ax.margins(y=0.24)
    ax.tick_params(axis="both", which="both", length=0)
    ax.tick_params(axis="x", labelsize=tick_fontsize)

    xlab_fs = label_fontsize if xlabel_fontsize is None else xlabel_fontsize
    ylab_fs = label_fontsize if ylabel_fontsize is None else ylabel_fontsize

    ax.set_xlabel(x_label, fontsize=xlab_fs, labelpad=xlabel_labelpad)
    if xlabel_x_position is not None or xlabel_y_position is not None:
        x_pos = xlabel_x_position if xlabel_x_position is not None else 0.5
        if xlabel_y_position is not None:
            y_pos = xlabel_y_position
        else:
            fig.canvas.draw()
            ax_h_pts = ax.get_position().height * fig.get_figheight() * 72
            tick_h_pts = (
                ax.xaxis.get_ticklabels()[0].get_fontsize()
                if ax.xaxis.get_ticklabels() else 12
            )
            y_pos = -(tick_h_pts * 1.5 + xlabel_labelpad) / ax_h_pts
        ax.xaxis.set_label_coords(x_pos, y_pos)

    ax.set_ylabel(ylabel, fontsize=ylab_fs)
    if title:
        ax.set_title(title, fontsize=title_fontsize)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    fig.tight_layout()

    if savepath:
        base = Path(savepath)
        if base.suffix.lower() in {".png", ".svg", ".pdf", ".jpg", ".jpeg"}:
            base = base.with_suffix("")
        out_dir = base.parent
        out_dir.mkdir(parents=True, exist_ok=True)
        fig_path = out_dir / f"{base.name}.{savetype}"
        fig.savefig(fig_path, dpi=600 if savetype == "png" else None, bbox_inches="tight")
        # Also save vector SVG for true-vector figure assembly
        svg_path = out_dir / f"{base.name}.svg"
        fig.savefig(svg_path, format='svg', bbox_inches="tight")
        print(f"Saved counterfactual plot: {fig_path}")

    if show_plot:
        plt.show()
    else:
        plt.close(fig)

    return work


__all__ = ["counterfactual_forest_plot"]
