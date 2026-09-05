"""
caps_item_distributions_hppd_split.py

Generates the 2×2 descriptive figure of CAPS vision item distributions split
by PPA history (persist_vis_yn) — manuscript Figure 3.

  caps_item_distributions_hppd_split.png
      Panel a is a binary endorsement count + %.  Panels b–d use vertical split
      violins conditioned on endorsement (structural zeros excluded).  Rows
      where any of the 6 binary CAPS items is null are dropped before
      plotting.

Sample (as of Jan 2026 data):
    SP users with valid CAPS: n = 130  (PPA(−) = 38, PPA(+) = 92)

Panel layout:
    a (top-left)     Binary endorsement count + % for
                     caps_bl_{4,26,31,23,19,22}, PPA(−)/PPA(+) grouped
                     bars side by side per item.

    b (top-right)    Frequency  (caps_bl_{x}c)
    c (bottom-left)  Intrusiveness (caps_bl_{x}b)
    d (bottom-right) Distress    (caps_bl_{x}a)

X-tick labels:
    Panel a:     plot_label from VARIABLE_REGISTRY.
    Panels b–d:  wrapped plot_label + "\\n{dimension}"
                 (e.g. "Shapes, Lights, Colors\\nFrequency").

Violin drawing notes:
    • Values are filtered to > 0 (conditional on item endorsement).
    • At each item x-position: left half = PPA(−), right half = PPA(+).
    • Groups with n ≥ 5 endorsed get a KDE half-violin (bw ≥ 0.5).
    • Groups with 1 ≤ n < 5 get individual scatter dots instead.
    • Groups with n = 0 are left blank — visually communicates rarity.
    • No quartile/IQR marks in this version.

Color scheme:
    binary_palette from master_config: [0] = PPA(−), [1] = PPA(+).

Usage (standalone):
    cd 04_visualizations
    /usr/local/bin/python3.12 caps_item_distributions_hppd_split.py

Usage (called from 0X_all_figures.py):
    from caps_item_distributions_hppd_split import (
        make_caps_item_distributions_hppd_split,
    )
    make_caps_item_distributions_hppd_split(
        df_sp=df_sp,
        binary_palette=binary_palette,
        hppd_neg_label=HPPD_NEG_LABEL,
        hppd_pos_label=HPPD_POS_LABEL,
        save_path=FIG_DIR / "caps_item_distributions_hppd_split.png",
    )
"""

import sys
import textwrap
from pathlib import Path

# Resolve modules/ relative to this file so the script works when invoked
# from any working directory (0X_all_figures.py os.chdir's to 04_visualizations/).
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent / "modules"))

import numpy as np
import scipy.stats as sci_stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

from master_config import VARIABLE_REGISTRY

###############################################################################
# Shared configuration
###############################################################################

# CAPS vision items to display, in the requested order.
CAPS_ITEMS = [4, 26, 31, 23, 19, 22]

# Binary endorsement columns (value 0 = No, 1 = Yes).
BINARY_COLS = [f"caps_bl_{x}" for x in CAPS_ITEMS]

# Sub-dimension panels: (column_suffix, display_label).
# c = Frequency, b = Distraction, a = Distress.
SUBDIM_PANELS = [
    ("c", "Frequency"),
    ("b", "Distraction"),
    ("a", "Distress"),
]

# Ordinal y-tick labels for violin panels.
# Frequency (c): REDCap scale 1–5 for how often the experience occurs.
# Distraction/Distress (b/a): REDCap scale 1–5 for how bothersome/distressing.
_YTICK_LABELS_FREQ = [
    "Hardly\nat all",
    "Not\nOften",
    "Sometimes",
    "Regularly",
    "All the\nTime",
]
_YTICK_LABELS_DIST = [
    "Not at\nall",
    "Slightly",
    "Somewhat",
    "Firmly",
    "Very",
]

# Fallback group labels; overridden by 0X_all_figures.py via the sobp_poster toggle.
_DEFAULT_NEG_LABEL = "PPA (-)"
_DEFAULT_POS_LABEL = "PPA (+)"

# Minimum KDE violin threshold.  Groups below this use scatter dots instead.
_KDE_MIN_N = 5

# Shared font sizes and spacing — applied identically to all panels.
_FS_AXIS  = 19.0   # y-axis label
_FS_TICK  = 9.8    # tick labels (both x and y) — 8.5 × 1.15
_FS_BAR   = 7.2    # bar count/% annotations
_YLABEL_PAD = 8    # labelpad for all y-axis labels

# Default output paths (relative to project root).
_DEFAULT_SAVE = (
    _HERE.parent / "results" / "descriptive"
    / "caps_item_distributions_hppd_split.png"
)

###############################################################################
# Shared helpers
###############################################################################

def _darken(color, factor=0.60):
    """Return a darkened version of *color* (used for IQR/median line colour)."""
    r, g, b, *_ = mcolors.to_rgba(color)
    return (r * factor, g * factor, b * factor)


def _prep_data(df_sp):
    """
    Drop rows with any null in the six binary CAPS columns and split by
    persist_vis_yn (0 = PPA−, 1 = PPA+).  Returns (df, df_neg, df_pos).

    Uses persist_vis_yn rather than hppd_binary because hppd_binary may have
    been remapped to display strings by 0X_all_figures.py before this function
    is called.
    """
    df = df_sp.dropna(subset=BINARY_COLS).copy()
    return df, df[df["persist_vis_yn"] == 0], df[df["persist_vis_yn"] == 1]


def _build_figure_skeleton(hppd_neg_label, hppd_pos_label, col_neg, col_pos,
                            n_neg, n_pos):
    """Shared 2×2 figure + GridSpec with shared legend at top-right."""
    plt.rcParams.update({
        "font.family":     "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
    })
    fig = plt.figure(figsize=(15, 11))
    gs  = gridspec.GridSpec(
        2, 2, figure=fig,
        left=0.08, right=0.97, top=0.91, bottom=0.22,
        wspace=0.38, hspace=0.62,
    )
    fig.legend(
        handles=[
            mpatches.Patch(facecolor=col_neg, alpha=0.86,
                           label=f"{hppd_neg_label} (n={n_neg})"),
            mpatches.Patch(facecolor=col_pos, alpha=0.86,
                           label=f"{hppd_pos_label} (n={n_pos})"),
        ],
        loc="upper right", fontsize=20, frameon=False, ncol=2,
        bbox_to_anchor=(0.97, 0.995),
    )
    return fig, gs


def _draw_panel_a(ax, df_neg, df_pos, n_neg, n_pos, plot_labels,
                  col_neg, col_pos, hppd_neg_label, hppd_pos_label):
    """Binary endorsement grouped bar chart (shared by both figure variants)."""
    xs, w = np.arange(len(CAPS_ITEMS)), 0.38

    def _clean(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(length=0)

    counts_neg = np.array([(df_neg[c] == 1).sum() for c in BINARY_COLS])
    counts_pos = np.array([(df_pos[c] == 1).sum() for c in BINARY_COLS])
    top_a = max(int(counts_neg.max()), int(counts_pos.max()), 1)

    bars_neg = ax.bar(xs - w / 2, counts_neg, w,
                      color=col_neg, alpha=0.86, edgecolor="white", lw=0.6,
                      label=hppd_neg_label, zorder=3)
    bars_pos = ax.bar(xs + w / 2, counts_pos, w,
                      color=col_pos, alpha=0.86, edgecolor="white", lw=0.6,
                      label=hppd_pos_label, zorder=3)

    for bars, counts, tot in [(bars_neg, counts_neg, n_neg),
                               (bars_pos, counts_pos, n_pos)]:
        for bar, cnt in zip(bars, counts):
            if cnt > 0:
                xc   = bar.get_x() + bar.get_width() / 2
                y_pct = bar.get_height() + top_a * 0.022
                # Percentage — bold, just above bar
                ax.text(xc, y_pct,
                        f"({cnt / tot * 100:.0f}%)",
                        ha="center", va="bottom", fontsize=_FS_BAR,
                        color="#333", fontweight="bold", clip_on=False)
                # Count — normal weight, above the percentage
                ax.text(xc, y_pct + top_a * 0.060,
                        f"{cnt}",
                        ha="center", va="bottom", fontsize=_FS_BAR,
                        color="#333", clip_on=False)

    ax.set_xticks(xs)
    ax.set_xticklabels(
        ["\n".join(textwrap.wrap(lbl, 32)) for lbl in plot_labels],
        rotation=35, ha="right", fontsize=_FS_TICK,
    )
    ax.set_ylabel("Endorsement", fontsize=_FS_AXIS, labelpad=_YLABEL_PAD)
    ax.set_ylim(0, top_a * 1.55)
    ax.tick_params(axis="y", labelsize=_FS_TICK)
    _clean(ax)


###############################################################################
# PRIMARY figure — vertical split violins (conditional on endorsement)
###############################################################################

def make_caps_item_distributions_hppd_split(
    df_sp,
    binary_palette,
    hppd_neg_label=_DEFAULT_NEG_LABEL,
    hppd_pos_label=_DEFAULT_POS_LABEL,
    save_path=None,
):
    """
    Create the PRIMARY 2×2 CAPS item distribution figure (violin version).

    Panels b–d show vertical split half-violins for Frequency, Intrusiveness,
    and Distress sub-dimensions among participants who *endorsed* each item
    (structural zeros excluded).

    Groups with n < 5 endorsers draw individual scatter dots instead of a
    violin — this is typically the PPA(−) group, and the sparsity itself is
    informative.

    Parameters
    ----------
    df_sp : pd.DataFrame
        SP-user subsample (psycheduse_yn == "Yes").
    binary_palette : list
        [neg_colour, pos_colour] from master_config.binary_palette.
    hppd_neg_label : str
        Display label for PPA(−) group.
    hppd_pos_label : str
        Display label for PPA(+) group.
    save_path : str or Path, optional
        Output file.  Defaults to
        results/descriptive/caps_item_distributions_hppd_split.png.
    """
    if save_path is None:
        save_path = _DEFAULT_SAVE
    save_path = Path(save_path)

    df, df_neg, df_pos = _prep_data(df_sp)
    n_neg, n_pos = len(df_neg), len(df_pos)
    print(f"  CAPS item distributions (violin): {len(df)} SP users with valid CAPS "
          f"({hppd_neg_label}: {n_neg}, {hppd_pos_label}: {n_pos})")

    plot_labels = [VARIABLE_REGISTRY[col]["plot_label"] for col in BINARY_COLS]
    col_neg = binary_palette[0]
    col_pos = binary_palette[1]

    fig, gs = _build_figure_skeleton(
        hppd_neg_label, hppd_pos_label, col_neg, col_pos, n_neg, n_pos,
    )

    lbl_kw = dict(fontsize=26, fontweight="bold", fontfamily="Arial",
                  va="bottom", ha="left", clip_on=False)
    panel_letters = iter("abcd")

    # ── Panel a ───────────────────────────────────────────────────────────
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.text(-0.13, 1.04, f"{next(panel_letters)}", transform=ax_a.transAxes, **lbl_kw)
    _draw_panel_a(ax_a, df_neg, df_pos, n_neg, n_pos, plot_labels,
                  col_neg, col_pos, hppd_neg_label, hppd_pos_label)

    # ── Panels b, c, d — vertical split violins ──────────────────────────
    panel_positions = [(0, 1), (1, 0), (1, 1)]
    xs = np.arange(len(CAPS_ITEMS))
    half_w = 0.40     # violin half-width in data units; gap ≈ 0.20 between items
    alpha  = 0.78

    def _vert_half_violin(ax, x_ctr, data, color, flip_left=False):
        """
        Draw a vertical half-violin at x_ctr.

        flip_left=True  → density extends LEFT  (PPA(−) side)
        flip_left=False → density extends RIGHT (PPA(+) side)

        Falls back to jittered scatter dots when n < _KDE_MIN_N.
        Does nothing when n == 0.
        """
        data = np.asarray(data, float)
        data = data[np.isfinite(data) & (data > 0)]  # conditional on endorsement
        n = len(data)

        if n == 0:
            return

        if n < _KDE_MIN_N:
            # Too few for a meaningful KDE — show individual observations as dots.
            # Horizontal jitter uses a fixed seed so the figure is reproducible.
            rng = np.random.default_rng(seed=int(x_ctr * 100 + (0 if flip_left else 1)))
            jitter = rng.uniform(-0.06, 0.06, n)
            # Offset dots toward their respective side so they don't straddle center.
            side_offset = -0.10 if flip_left else 0.10
            ax.scatter(x_ctr + side_offset + jitter, data,
                       color=color, alpha=0.75, s=22, zorder=5,
                       edgecolors="white", linewidths=0.3)
            return

        # KDE — use a minimum bandwidth of 0.5 so ordinal 1–5 values smooth well.
        bw = max(0.5, sci_stats.gaussian_kde(data).factor)
        kde = sci_stats.gaussian_kde(data, bw_method=bw)
        kde_max = float(kde(np.linspace(data.min(), data.max(), 400)).max())

        # Evaluate over the full scale range [0.8, 5.2] so violin tails look complete.
        y_lo = max(0.8, data.min() - 0.3)
        y_hi = min(5.2, data.max() + 0.3)
        ys      = np.linspace(y_lo, y_hi, 400)
        density = kde(ys)
        scaled  = density / density.max() * half_w

        if flip_left:
            ax.fill_betweenx(ys, x_ctr - scaled, x_ctr,
                             color=color, alpha=alpha, zorder=2)
        else:
            ax.fill_betweenx(ys, x_ctr, x_ctr + scaled,
                             color=color, alpha=alpha, zorder=2)

        # IQR / median lines — horizontal lines at Q25, Q50, Q75.
        # Length proportional to KDE density at that quantile (same normalisation
        # as the violin fill).  Median is solid+wider; quartiles are dashed+thinner.
        dark = _darken(color)
        q25, q50, q75 = np.percentile(data, [25, 50, 75])
        for q, ls, lw in [(q25, "--", 0.9), (q50, "-", 1.7), (q75, "--", 0.9)]:
            reach = float(kde(np.array([q]))[0]) / kde_max * half_w
            if flip_left:
                ax.plot([x_ctr - reach, x_ctr], [q, q],
                        color=dark, ls=ls, lw=lw, alpha=0.9, zorder=5)
            else:
                ax.plot([x_ctr, x_ctr + reach], [q, q],
                        color=dark, ls=ls, lw=lw, alpha=0.9, zorder=5)

    for (row, col_i), (suffix, dim_label) in zip(panel_positions, SUBDIM_PANELS):
        ax = fig.add_subplot(gs[row, col_i])
        ax.text(-0.13, 1.04, f"{next(panel_letters)}",
                transform=ax.transAxes, **lbl_kw)

        for i, (x_ctr, item) in enumerate(zip(xs, CAPS_ITEMS)):
            bin_col = f"caps_bl_{item}"
            sub_col = f"caps_bl_{item}{suffix}"

            # Conditional values: endorsed items only (sub-dim value > 0)
            neg_vals = df_neg.loc[df_neg[bin_col] == 1, sub_col].dropna()
            pos_vals = df_pos.loc[df_pos[bin_col] == 1, sub_col].dropna()

            _vert_half_violin(ax, float(x_ctr), neg_vals, col_neg, flip_left=True)
            _vert_half_violin(ax, float(x_ctr), pos_vals, col_pos, flip_left=False)

            # Thin center divider to visually separate the two halves
            ax.axvline(x_ctr, color="#cccccc", lw=0.6, zorder=1)

        ax.set_xticks(xs)
        ax.set_xticklabels(
            ["\n".join(textwrap.wrap(lbl, 32)) + f"\n{dim_label}"
             for lbl in plot_labels],
            rotation=35, ha="right", fontsize=_FS_TICK,
        )
        ax.set_xlim(-0.6, len(CAPS_ITEMS) - 0.4)

        # Y-axis: ordinal text tick labels + wrapped ylabel (one word per line).
        # labelpad large enough so the ylabel clears the text tick labels.
        ytick_labels = _YTICK_LABELS_FREQ if suffix == "c" else _YTICK_LABELS_DIST
        ax.set_ylim(0.5, 5.5)
        ax.set_yticks([1, 2, 3, 4, 5])
        ax.set_yticklabels(ytick_labels, fontsize=_FS_TICK)
        ylabel_text = "\n".join(f"{dim_label} Rating".split())
        ax.set_ylabel(ylabel_text, fontsize=_FS_AXIS, labelpad=_YLABEL_PAD)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(length=0)

    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {save_path}")


###############################################################################

if __name__ == "__main__":
    import pandas as pd
    from master_config import binary_palette as _bp
    from master_config import SP_USER_COL, SP_USER_VALUE

    # The shipped analysis dataframe is already fully prepared; read it directly.
    from data_prep import most_recent_public_df
    df = pd.read_csv(
        most_recent_public_df(_HERE.parent / "data" / "final"), low_memory=False
    )
    df_sp = df[df[SP_USER_COL] == SP_USER_VALUE].copy()   # canonical SP-user filter

    make_caps_item_distributions_hppd_split(
        df_sp=df_sp, binary_palette=_bp,
    )
