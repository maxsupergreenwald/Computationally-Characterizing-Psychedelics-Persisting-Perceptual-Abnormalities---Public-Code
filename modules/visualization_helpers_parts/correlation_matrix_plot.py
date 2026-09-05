"""Helper function `correlation_matrix_plot`."""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def correlation_matrix_plot(
    row_cols,
    column_cols,
    dataframe,
    palette=None,
    color_axis='x',
    color_var=None,
    line_intensity=80,
    sig_fontsize=40,
    xlab_fontsize=16,
    ylab_fontsize=16,
    figsize=None,
    dv_to_label_dict=None,
    savepath=None,
    show_n=None,
    tick_size=12,
    dpi=400,
    add_x_jitter=None,
    scatter_size=50,
    show_xlabel=True,
    line_lw=2,
    control_vars=None,
    ylab_y=None,
    pad_inches=0.1,
):
    """
    Create a multipanel correlation matrix figure with each subplot showing
    scatter plot + regression line + significance markers.

    Parameters
    ----------
    row_cols : list of str
        Column names for variables plotted on rows (y-axis in each subplot)
    column_cols : list of str
        Column names for variables plotted on columns (x-axis in each subplot)
    dataframe : pd.DataFrame
        DataFrame containing all variables
    palette : str | list | matplotlib colormap, optional
        Color palette for scatter points. If list, can be treated as a discrete palette
        when the colored axis has integer levels (e.g., caps_vision 0-6).
    color_axis : {'x','y'}
        Which axis to use for point coloring.  Ignored when color_var is set.
    color_var : str, optional
        Column name to use for point coloring instead of the x or y axis variable.
        When set, this column is included in the NaN-drop and its values drive the
        scatter color.  Useful when the color dimension is a third variable not on
        either axis (e.g. color by caps_vision while x = vch_beta, y = sdt vars).
        Default None (falls back to color_axis logic).
    line_intensity : int, optional
        Intensity of regression line color (0-100). Default is 40
    sig_fontsize : int, optional
        Font size for significance markers. Default is 25
    figsize : tuple, optional
        Figure size (width, height). If None, auto-calculated based on number of panels
    dv_to_label_dict : dict, optional
        Dictionary mapping variable names to display labels.
        Defaults to dv_to_lab_short if None
    savepath : str, optional
        Path to save figure. Default is None (no save)
    dpi : int, optional
        Resolution for saved figure. Default is 400
    add_x_jitter : float, optional
        Amount of jitter to add to x-axis values. If None, no jitter is added.
        Jitter is added as uniform random noise in range [-add_x_jitter, add_x_jitter].
        Default is None
    show_xlabel : bool, optional
        Whether to show the x-axis label on the bottom row. Set to False when
        compositing into a figure where a lower panel provides the axis label.
        Default is True.
    line_lw : float, optional
        Line width for the regression line. Default is 2.
    control_vars : list of str, optional
        If provided, compute partial Spearman correlation controlling for these
        variables.  All variables (row, col, control_vars) are rank-transformed,
        then each is residualised on the control ranks via OLS.  The scatter
        shows the rank-residual partial regression plot and significance stars
        are derived from the partial Spearman rho.  add_x_jitter is ignored
        when control_vars is set (rank residuals are continuous).  Default None.
    ylab_y : float, optional
        Vertical position of y-axis labels in axes coordinates (0 = bottom,
        0.5 = center, 1 = top).  Default None uses matplotlib's default (0.5).
        Lower values shift labels downward, useful when the top label is
        clipped by the panel boundary after assembly.

    Returns
    -------
    fig, axes
        Matplotlib figure and axes objects
    """
    import os
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap

    # Use dv_to_lab_short if no custom dict provided
    if dv_to_label_dict is None:
        global dv_to_lab_short
        if "dv_to_lab_short" not in globals():
            from modules.master_config import dv_to_lab_short
        dv_to_label_dict = dv_to_lab_short

    def _default_palette():
        # Prefer caps_vision palette when plotting caps_vision on the color axis.
        # Also triggered when color_var='caps_vision' (third-variable coloring).
        try:
            from master_config import caps_vision_palette
        except Exception:
            try:
                from modules.master_config import caps_vision_palette
            except Exception:
                try:
                    from ..variable_labels import caps_vision_palette
                except Exception:
                    caps_vision_palette = None

        if color_var == 'caps_vision' and caps_vision_palette is not None:
            return caps_vision_palette
        if color_axis == 'x' and 'caps_vision' in column_cols and caps_vision_palette is not None:
            return caps_vision_palette
        if color_axis == 'y' and 'caps_vision' in row_cols and caps_vision_palette is not None:
            return caps_vision_palette
        return 'Oranges'

    def _as_cmap(pal):
        if hasattr(pal, '__call__') and hasattr(pal, 'N'):
            return pal
        if isinstance(pal, str):
            return plt.cm.get_cmap(pal)
        if isinstance(pal, (list, tuple)):
            return LinearSegmentedColormap.from_list('custom_list', list(pal))
        return plt.cm.get_cmap('Oranges')

    if palette is None:
        palette = _default_palette()

    nrows = len(row_cols)
    ncols = len(column_cols)

    # Auto-calculate figure size if not provided
    if figsize is None:
        fig_width = max(4 * ncols, 8)
        fig_height = max(3 * nrows, 6)
        figsize = (fig_width, fig_height)

    # Create subplot grid
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)

    # Get colormap for line color
    cmap = _as_cmap(palette)
    line_color = cmap(line_intensity / 100)
    summary_rows = []

    # Loop through each subplot
    for r, row_var in enumerate(row_cols):
        for c, col_var in enumerate(column_cols):
            ax = axes[r, c]

            if control_vars is not None:
                # ── Partial Spearman: raw scatter, partial rho annotation ──
                #
                # Display strategy (option 4 of 4):
                #   - Scatter and regression line use RAW x/y values so the
                #     axes remain on the original scale (e.g. CAPS 0-6), which
                #     is interpretable to readers.
                #   - rho and p-value come from the partial Spearman computed
                #     on ranked data (rank residuals via pingouin.partial_corr),
                #     so the significance stars reflect the age-controlled
                #     association, not the marginal one.
                #
                # Rejected alternatives (kept here for reference):
                #   Option 1 — rank-residual scatter: statistically purest
                #     partial regression plot, but x-axis shows opaque residuals
                #     (~-80 to +80) rather than the original scale. Hard for
                #     readers to interpret without explanation.
                #   Option 2 — rank-residual scatter + clarified axis label:
                #     same as option 1 but label reads "CAPS Vision
                #     (age-residualized ranks)". Still requires reader effort.
                #   Option 3 — raw scatter, regression line from added-variable
                #     logic: scatter raw values, draw a partial regression line
                #     (slope from partial rho rescaled to raw units). More
                #     complex to implement correctly; marginal benefit over
                #     option 4.
                _ctrl_list = list(control_vars)
                _all_vars = [col_var, row_var] + _ctrl_list
                # Include color_var in dropna when it's a third variable
                if color_var is not None and color_var not in _all_vars:
                    _all_vars = _all_vars + [color_var]
                _clean = dataframe[_all_vars].dropna().copy()

                if len(_clean) < len(_ctrl_list) + 3:
                    ax.text(0.5, 0.5, 'N/A', transform=ax.transAxes,
                           ha='center', va='center', fontsize=20, color='gray')
                    ax.set_xticks([])
                    ax.set_yticks([])
                    continue

                from scipy.stats import rankdata as _rankdata
                _ranked = pd.DataFrame(
                    {v: _rankdata(_clean[v]) for v in _all_vars},
                    index=_clean.index,
                )

                # Partial Spearman = partial Pearson on ranked data (pingouin)
                try:
                    import pingouin as _pg
                    _covar_arg = _ctrl_list if len(_ctrl_list) > 1 else _ctrl_list[0]
                    _pc = _pg.partial_corr(
                        data=_ranked, x=col_var, y=row_var,
                        covar=_covar_arg, method='pearson',
                    )
                    rho = float(_pc['r'].iloc[0])
                    p_val = float(_pc['p-val'].iloc[0])
                except Exception:
                    from scipy.stats import pearsonr as _pearsonr
                    rho, p_val = _pearsonr(_ranked[col_var], _ranked[row_var])

                n = len(_ranked)

                # Use raw values for the scatter so axes stay on original scale.
                # Include color_var in temp_df when it is a third column.
                _scatter_cols = [col_var, row_var]
                if color_var is not None and color_var not in _scatter_cols:
                    _scatter_cols = _scatter_cols + [color_var]
                temp_df = _clean[_scatter_cols].copy()
                if add_x_jitter is not None:
                    jitter = np.random.uniform(
                        -add_x_jitter, add_x_jitter, size=len(temp_df)
                    )
                    temp_df[col_var + '_jittered'] = temp_df[col_var] + jitter
                    x_col = col_var + '_jittered'
                else:
                    x_col = col_var

            else:
                # ── Standard Spearman ──────────────────────────────────────
                # Remove NaN values for this pair (include color_var when set)
                _base_cols = [col_var, row_var]
                if color_var is not None and color_var not in _base_cols:
                    _base_cols = [col_var, row_var, color_var]
                temp_df = dataframe[_base_cols].dropna().copy()

                if len(temp_df) < 3:  # Need at least 3 points for correlation
                    # Show empty panel with "N/A" text
                    ax.text(0.5, 0.5, 'N/A', transform=ax.transAxes,
                           ha='center', va='center', fontsize=20, color='gray')
                    ax.set_xticks([])
                    ax.set_yticks([])
                    continue

                # Calculate Spearman correlation
                from scipy.stats import spearmanr
                rho, p_val = spearmanr(temp_df[col_var], temp_df[row_var])
                n = len(temp_df)

                # Add jitter to x-axis if requested
                if add_x_jitter is not None:
                    jitter = np.random.uniform(-add_x_jitter, add_x_jitter, size=len(temp_df))
                    temp_df[col_var + '_jittered'] = temp_df[col_var] + jitter
                    x_col = col_var + '_jittered'
                else:
                    x_col = col_var

            # Determine which variable drives point color.
            # color_var (parameter) takes priority; otherwise fall back to color_axis.
            _color_col = color_var if color_var is not None else (
                row_var if color_axis == 'y' else col_var
            )
            color_values = temp_df[_color_col]

            if isinstance(palette, (list, tuple)):
                # Discrete coloring when values are integer levels within palette range
                int_vals = color_values.round().astype(int)
                is_int = np.allclose(color_values, int_vals, atol=1e-6)
                in_range = int_vals.min() >= 0 and int_vals.max() < len(palette)
                if is_int and in_range:
                    temp_df['dv_discrete'] = int_vals
                    palette_map = {i: palette[i] for i in range(len(palette))}
                    sns.scatterplot(
                        data=temp_df,
                        x=x_col,
                        y=row_var,
                        hue='dv_discrete',
                        palette=palette_map,
                        s=scatter_size,
                        alpha=0.6,
                        ax=ax,
                        legend=False,
                    )
                else:
                    # Fall back to continuous coloring with a colormap
                    if color_values.max() != color_values.min():
                        temp_df['dv_normalized'] = (color_values - color_values.min()) / (
                            color_values.max() - color_values.min()
                        )
                    else:
                        temp_df['dv_normalized'] = 0.5
                    sns.scatterplot(
                        data=temp_df,
                        x=x_col,
                        y=row_var,
                        hue='dv_normalized',
                        palette=_as_cmap(palette),
                        s=scatter_size,
                        alpha=0.6,
                        ax=ax,
                        legend=False,
                    )
            else:
                if color_values.max() != color_values.min():
                    temp_df['dv_normalized'] = (color_values - color_values.min()) / (
                        color_values.max() - color_values.min()
                    )
                else:
                    temp_df['dv_normalized'] = 0.5
                sns.scatterplot(
                    data=temp_df,
                    x=x_col,
                    y=row_var,
                    hue='dv_normalized',
                    palette=palette,
                    s=scatter_size,
                    alpha=0.6,
                    ax=ax,
                    legend=False,
                )

            # Regression line with selected color (always use original x values, not jittered)
            # ci=None removes the shaded confidence band while keeping the regression line.
            sns.regplot(
                data=temp_df,
                x=col_var,
                y=row_var,
                scatter_kws={'alpha': 0},
                line_kws={'color': line_color, 'lw': line_lw},
                ci=None,
                ax=ax,
            )

            # Add significance stars
            if p_val < 0.001:
                sig_text = '***'
            elif p_val < 0.01:
                sig_text = '**'
            elif p_val < 0.05:
                sig_text = '*'
            elif p_val < 0.1:
                sig_text = '~'
            else:
                sig_text = ''

            summary_rows.append(
                {
                    "x_var": col_var,
                    "y_var": row_var,
                    "row_var": row_var,
                    "column_var": col_var,
                    "row_label": dv_to_label_dict.get(row_var, row_var),
                    "column_label": dv_to_label_dict.get(col_var, col_var),
                    "rho": rho,
                    "p_value": p_val,
                    "n": n,
                    "sig_text": sig_text,
                }
            )

            if sig_text:
                ax.text(
                    0.5,
                    0.95,
                    sig_text,
                    transform=ax.transAxes,
                    verticalalignment='top',
                    horizontalalignment='center',
                    fontsize=sig_fontsize,
                    fontweight='bold',
                )

            # Add n in corner if specified
            if show_n:
                ax.text(
                    0.98,
                    0.02,
                    f'n={n}',
                    transform=ax.transAxes,
                    verticalalignment='bottom',
                    horizontalalignment='right',
                    fontsize=8,
                    color='gray',
                )

            # Set labels using dv_to_label_dict
            # Only show x-labels on bottom row, and only when show_xlabel=True
            if r == nrows - 1 and show_xlabel:
                xlabel = dv_to_label_dict.get(col_var, col_var)
                ax.set_xlabel(xlabel, fontsize=xlab_fontsize)
            else:
                ax.set_xlabel('')

            # Only show y-labels on left column
            if c == 0:
                ylabel = dv_to_label_dict.get(row_var, row_var)
                # Wrap long y-axis labels so they don't overflow the figure
                import textwrap as _tw
                ylabel = "\n".join(_tw.wrap(str(ylabel), width=20))
                _ylabel_kw = {'fontsize': ylab_fontsize}
                if ylab_y is not None:
                    _ylabel_kw['y'] = ylab_y
                ax.set_ylabel(ylabel, **_ylabel_kw)
            else:
                ax.set_ylabel('')

            # Remove spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

            # Tick parameters
            ax.tick_params(axis="both", which="both", length=0, labelsize=tick_size)

    plt.tight_layout()

    # Save figure if savepath provided
    if savepath:
        os.makedirs(os.path.dirname(savepath), exist_ok=True)
        fig.savefig(savepath, dpi=dpi, bbox_inches='tight', pad_inches=pad_inches)
        # Also save vector SVG for true-vector figure assembly
        _svg_path = savepath.replace('.png', '.svg').replace('.pdf', '.svg')
        fig.savefig(_svg_path, format='svg', bbox_inches='tight', pad_inches=pad_inches)

    # Save summary statistics if savepath provided
    if savepath:
        summary_dir = os.path.join(os.path.dirname(savepath), 'summary_results')
        os.makedirs(summary_dir, exist_ok=True)
        summary_path = os.path.join(
            summary_dir,
            os.path.basename(savepath).replace('.png', '').replace('.pdf', '') + '.csv',
        )
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(summary_path, index=False)

    return fig, axes
