"""Helper function `multipanel_boxplot_grid`."""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from master_config import dv_to_lab_short

def multipanel_boxplot_grid(
    data,
    panel_specs,  # List of dicts: [{'dv': 'var1', 'group_var': 'group1', 'order': [...], 'palette': [...], 'ylabel': 'Label', 'xtick_labels': [...]}, ...]
    nrows=2,
    ncols=3,
    panel_width=7,
    panel_height=5.5,
    strip_alpha=0.4,
    strip_size=4,
    strip_color="red",

    box_alpha=0.2,
    box_width=0.45,

    show_stats=False,
    show_pairwise=True,
    pairwise_line_style='--',
    ylabel_fontsize=20,
    xlabel_fontsize=22,
    title_fontsize=24,
    sig_marker_fontsize=40,
    pairwise_sig_modifier=0.75,
    mean_marker='o',
    mean_marker_size=20,
    line_width=3,
    ytick_font_multiplier=0.7,
    savepath=None,
    ylim=None,
    show_outliers=False,
    xlabel=True,
    xtick_font_multiplier=0.75,
    showcaps=False,
    sig_marker_x=0.5,
    pairwise_line_height_increment=0.35,
    sig_marker_y=0.9,
    bonferroni_correction=None,
    dpi=300,
    highlight_col=None,
    highlight_color='#CC0000',
    ytick_fontsize=None,
    xtick_fontsize=None,
):
    """
    Create an A x B grid of boxplots where each panel can have a different DV and group variable.
    
    Parameters:
    -----------
    data : DataFrame
        The data to plot
    panel_specs : list of dicts
        Each dict specifies a panel with keys:
        - 'dv': dependent variable name
        - 'group_var': grouping variable name
        - 'order': (optional) order of groups
        - 'palette': (optional) color palette for this panel
        - 'ylabel': (optional) y-axis label
        - 'xtick_labels': (optional) custom x-tick labels
    nrows, ncols : int
        Number of rows and columns in the grid
    savepath : str
        Path to save the figure (will save as PNG)
    """
    import os
    import numpy as np
    from scipy.stats import mannwhitneyu, kruskal
    import itertools
    
    n_panels = len(panel_specs)
    if n_panels == 0:
        raise ValueError("panel_specs must contain at least one panel specification")
    
    if n_panels > nrows * ncols:
        raise ValueError(f"panel_specs has {n_panels} panels but grid only has {nrows * ncols} positions")
    
    def _wrap_label(label):
        text = str(label)
        if "\n" in text:
            return text
        return "\n".join(text.split(" "))

    # Create figure with grid
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel_width * ncols, panel_height * nrows),
                             squeeze=False)
    
    results = {}
    
    for panel_idx, spec in enumerate(panel_specs):
        try:
            print(f"\n{'='*60}")
            print(f"Processing panel {panel_idx}: Row {panel_idx // ncols}, Col {panel_idx % ncols}")
            print(f"DV: {spec['dv']}, Group: {spec['group_var']}")
            
            # Calculate row and column for this panel
            row = panel_idx // ncols
            col = panel_idx % ncols
            ax = axes[row, col]
            
            # Extract panel specifications
            dv = spec['dv']
            group_var = spec['group_var']
            order = spec.get('order', None)
            palette_spec = spec.get('palette', ["#3498db", "#e74c3c"])
            ylabel = spec.get('ylabel', dv)
            xtick_labels = spec.get('xtick_labels', None)
            
            # Determine groups first (needed for palette conversion)
            clean_data_temp = data[[dv, group_var]].dropna()
            if order is None:
                unique_groups = sorted(clean_data_temp[group_var].unique())
            else:
                unique_groups = order
            
            n_groups = len(unique_groups)
            
            # Convert string palette to list of colors
            if isinstance(palette_spec, str):
                cmap = plt.cm.get_cmap(palette_spec)
                # Sample from 0.05 to 0.75 of the colormap (5% to 75%)
                palette = [cmap(0.05 + (i / (n_groups - 1) * 0.70)) for i in range(n_groups)]
            else:
                palette = palette_spec
            
            # Clean data for this DV (include highlight_col if available)
            _extra = [highlight_col] if (highlight_col and highlight_col in data.columns) else []
            clean_data = data[[dv, group_var] + _extra].dropna(subset=[dv, group_var])
            
            if len(clean_data) == 0:
                ax.text(0.5, 0.5, f"No data for {dv}", ha='center', va='center', transform=ax.transAxes)
                continue
            
            # Clip data if ylim provided
            if ylim is not None:
                clean_data = clean_data[(clean_data[dv] >= ylim[0]) & (clean_data[dv] <= ylim[1])].copy()
            
            # Determine groups
            if order is None:
                unique_groups = sorted(clean_data[group_var].unique())
            else:
                unique_groups = order
            
            n_groups = len(unique_groups)

            # seaborn 0.13.2 compatibility fix: convert list palette → dict keyed by
            # string group labels, and map the group_var column to those string labels.
            # Without this, seaborn's internal plot_boxes → _configure_legend raises
            # UnboundLocalError: cannot access local variable 'boxprops'.
            _pal_dict = (
                {g: c for g, c in zip(unique_groups, palette)}
                if isinstance(palette, (list, tuple)) else palette
            )
            _hue_col = f"__hue_label_{group_var}"
            clean_data = clean_data.copy()
            _existing_vals = set(clean_data[group_var].dropna().unique())
            if _existing_vals <= set(unique_groups):
                # Data already uses the target string labels
                clean_data[_hue_col] = clean_data[group_var].astype(str)
            elif clean_data[group_var].dtype.kind in 'iuf':
                # Integer/float indices → map 0 → groups[0], 1 → groups[1], …
                clean_data[_hue_col] = clean_data[group_var].map(
                    dict(enumerate(unique_groups))
                )
            else:
                clean_data[_hue_col] = clean_data[group_var].astype(str)

            # 2-group case: Mann-Whitney
            if n_groups == 2:
                groups = unique_groups

                # Create boxplot (no box outline; we draw quartile lines manually)
                # NOTE: legend=False triggers a seaborn 0.13 bug (boxprops UnboundLocalError).
                # Workaround: omit legend=False and remove the auto-created legend after.
                bp = sns.boxplot(
                    data=clean_data,
                    x=_hue_col,
                    y=dv,
                    hue=_hue_col,
                    palette=_pal_dict,
                    order=groups,
                    ax=ax,
                    showfliers=show_outliers,
                    showcaps=showcaps,
                    width=box_width,
                    whiskerprops={"linewidth": line_width},
                    capprops={"linewidth": line_width},
                    medianprops={"linewidth": line_width},
                )
                if ax.get_legend() is not None:
                    ax.get_legend().remove()

                # Set box alpha and remove edge (replaces the old boxprops linewidth=0)
                for patch in ax.patches:
                    patch.set_alpha(box_alpha)
                    patch.set_edgecolor("none")

                # Add strip plot with per-group colors
                for i, group in enumerate(groups):
                    group_df = clean_data[clean_data[_hue_col] == group]
                    sns.stripplot(
                        data=group_df,
                        x=_hue_col,
                        y=dv,
                        order=groups,
                        color=palette[i] if isinstance(palette, (list, tuple)) else list(_pal_dict.values())[i],
                        alpha=strip_alpha,
                        size=strip_size,
                        jitter=0.2,
                        ax=ax,
                    )
                    # Overplot highlighted points in red
                    if highlight_col and highlight_col in clean_data.columns:
                        hi_df = group_df[group_df[highlight_col] > 0]
                        if len(hi_df) > 0:
                            sns.stripplot(
                                data=hi_df,
                                x=_hue_col,
                                y=dv,
                                order=groups,
                                color=highlight_color,
                                alpha=0.9,
                                size=strip_size + 1.5,
                                jitter=0.2,
                                ax=ax,
                                zorder=5,
                            )

                # Add mean markers and quartile lines (matching box color)
                for i, group in enumerate(groups):
                    group_data = clean_data[clean_data[_hue_col] == group][dv]
                    mean_val = group_data.mean()
                    q1 = group_data.quantile(0.25)
                    q3 = group_data.quantile(0.75)
                    # Mean marker
                    ax.plot(
                        i,
                        mean_val,
                        marker=mean_marker,
                        markersize=mean_marker_size,
                        color=palette[i],
                        markeredgecolor=palette[i],
                        markeredgewidth=1.5,
                        zorder=3,
                        alpha=1,
                    )
                    # Quartile lines (25% and 75%)
                    half_width = box_width / 2
                    ax.hlines(
                        [q1, q3],
                        i - half_width,
                        i + half_width,
                        colors=palette[i],
                        linewidth=line_width,
                        alpha=1,
                        zorder=3,
                    )

                # Recolor whiskers/caps/median/fliers to match box color
                for line in ax.lines:
                    xdata = line.get_xdata()
                    if xdata is None or len(xdata) == 0:
                        continue
                    # Use mean x position to infer group index
                    x_mean = float(np.mean(xdata))
                    idx = int(np.clip(round(x_mean), 0, n_groups - 1))
                    line.set_color(palette[idx])
                    line.set_alpha(1)

                # Perform Mann-Whitney test
                group0_data = clean_data[clean_data[_hue_col] == groups[0]][dv]
                group1_data = clean_data[clean_data[_hue_col] == groups[1]][dv]
                statistic, p_value = mannwhitneyu(group0_data, group1_data, alternative='two-sided')
                n1 = len(group0_data)
                n2 = len(group1_data)
                denom = n1 * n2
                rank_biserial = ((2 * statistic) / denom - 1) if denom > 0 else np.nan

                g0_q1 = group0_data.quantile(0.25)
                g0_med = group0_data.median()
                g0_q3 = group0_data.quantile(0.75)
                g1_q1 = group1_data.quantile(0.25)
                g1_med = group1_data.median()
                g1_q3 = group1_data.quantile(0.75)
                
                results[f"{dv}_{group_var}"] = {
                    'test': 'mann_whitney',
                    'dv': dv,
                    'group_var': group_var,
                    'statistic': statistic,
                    'p_value': p_value,
                    'rank_biserial_correlation': rank_biserial,
                    'n_groups': 2,
                    'sample_size': len(clean_data),
                    'groups': [str(g) for g in groups],
                    'group_1_label': str(groups[0]),
                    'group_1_n': n1,
                    'group_1_median': g0_med,
                    'group_1_iqr_q1': g0_q1,
                    'group_1_iqr_q3': g0_q3,
                    'group_2_label': str(groups[1]),
                    'group_2_n': n2,
                    'group_2_median': g1_med,
                    'group_2_iqr_q1': g1_q1,
                    'group_2_iqr_q3': g1_q3,
                }
                
                # Show significance
                if not show_stats:
                    if p_value < 0.001:
                        sig_text = "***"
                    elif p_value < 0.01:
                        sig_text = "**"
                    elif p_value < 0.05:
                        sig_text = "*"
                    elif p_value < 0.1:
                        sig_text = "~"
                    else:
                        sig_text = ""
                    
                    ax.text(sig_marker_x, sig_marker_y, sig_text, transform=ax.transAxes,
                        verticalalignment='top', horizontalalignment='center',
                        fontsize=sig_marker_fontsize, fontweight='bold')
            
            # 3+ group case: Kruskal-Wallis
            else:
                groups = unique_groups

                # Create boxplot (no box outline; we draw quartile lines manually)
                # NOTE: legend=False triggers a seaborn 0.13 bug (boxprops UnboundLocalError).
                # Workaround: omit legend=False and remove the auto-created legend after.
                bp = sns.boxplot(
                    data=clean_data,
                    x=_hue_col,
                    y=dv,
                    hue=_hue_col,
                    palette=_pal_dict,
                    order=groups,
                    ax=ax,
                    showfliers=show_outliers,
                    showcaps=showcaps,
                    width=box_width,
                    whiskerprops={"linewidth": line_width},
                    capprops={"linewidth": line_width},
                    medianprops={"linewidth": line_width},
                )
                if ax.get_legend() is not None:
                    ax.get_legend().remove()

                # Set box alpha and remove edge (replaces the old boxprops linewidth=0)
                for patch in ax.patches:
                    patch.set_alpha(box_alpha)
                    patch.set_edgecolor("none")

                # Add strip plot with per-group colors
                for i, group in enumerate(groups):
                    group_df = clean_data[clean_data[_hue_col] == group]
                    sns.stripplot(
                        data=group_df,
                        x=_hue_col,
                        y=dv,
                        order=groups,
                        color=palette[i] if isinstance(palette, (list, tuple)) else list(_pal_dict.values())[i],
                        alpha=strip_alpha,
                        size=strip_size,
                        jitter=0.2,
                        ax=ax,
                    )
                    # Overplot highlighted points in red
                    if highlight_col and highlight_col in clean_data.columns:
                        hi_df = group_df[group_df[highlight_col] > 0]
                        if len(hi_df) > 0:
                            sns.stripplot(
                                data=hi_df,
                                x=_hue_col,
                                y=dv,
                                order=groups,
                                color=highlight_color,
                                alpha=0.9,
                                size=strip_size + 1.5,
                                jitter=0.2,
                                ax=ax,
                                zorder=5,
                            )

                # Add mean markers and quartile lines (matching box color)
                for i, group in enumerate(groups):
                    group_data = clean_data[clean_data[_hue_col] == group][dv]
                    mean_val = group_data.mean()
                    q1 = group_data.quantile(0.25)
                    q3 = group_data.quantile(0.75)
                    # Mean marker
                    ax.plot(
                        i,
                        mean_val,
                        marker=mean_marker,
                        markersize=mean_marker_size,
                        color=palette[i],
                        markeredgecolor=palette[i],
                        markeredgewidth=1.5,
                        zorder=3,
                        alpha=1,
                    )
                    # Quartile lines (25% and 75%)
                    half_width = box_width / 2
                    ax.hlines(
                        [q1, q3],
                        i - half_width,
                        i + half_width,
                        colors=palette[i],
                        linewidth=line_width,
                        alpha=1,
                        zorder=3,
                    )

                # Recolor whiskers/caps/median/fliers to match box color
                for line in ax.lines:
                    xdata = line.get_xdata()
                    if xdata is None or len(xdata) == 0:
                        continue
                    # Use mean x position to infer group index
                    x_mean = float(np.mean(xdata))
                    idx = int(np.clip(round(x_mean), 0, n_groups - 1))
                    line.set_color(palette[idx])
                    line.set_alpha(1)
                
                # Perform Kruskal-Wallis test
                group_data = [clean_data[clean_data[_hue_col] == g][dv].values for g in groups]
                h_stat, kw_p = kruskal(*group_data)
                
                results[f"{dv}_{group_var}"] = {
                    'test': 'kruskal_wallis',
                    'dv': dv,
                    'group_var': group_var,
                    'H': h_stat,
                    'p_value': kw_p,
                    'n_groups': n_groups,
                    'sample_size': len(clean_data),
                    'groups': [str(g) for g in groups],
                    'pairwise': {}
                }
                
                # Pairwise comparisons
                if kw_p < 0.05 and show_pairwise:
                    pairs = list(itertools.combinations(range(n_groups), 2))
                    n_comparisons = len(pairs)
                    alpha_corrected = 0.05 / n_comparisons if bonferroni_correction else 0.05
                    
                    sig_pairs = []
                    for i, j in pairs:
                        g1, g2 = groups[i], groups[j]
                        vals1 = clean_data[clean_data[_hue_col] == g1][dv].values
                        vals2 = clean_data[clean_data[_hue_col] == g2][dv].values
                        
                        u_stat, p_pairwise = mannwhitneyu(vals1, vals2, alternative='two-sided')
                        
                        results[f"{dv}_{group_var}"]['pairwise'][f"{g1} vs {g2}"] = {
                            'group_1': str(g1),
                            'group_2': str(g2),
                            'statistic': u_stat,
                            'p_value': p_pairwise
                        }
                        
                        if p_pairwise < 0.1:
                            sig_pairs.append((i, j, p_pairwise))
                    
                    # Draw significance lines
                    if not show_stats and sig_pairs:
                        y_max = ax.get_ylim()[1]
                        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
                        line_y_start = y_max - (pairwise_line_height_increment * y_range)
                        
                        for level_idx, (i, j, p_val_pair) in enumerate(sig_pairs):
                            y_line = line_y_start - (level_idx * pairwise_line_height_increment * y_range)
                            
                            ax.plot([i, j], [y_line, y_line],
                                color='black', linewidth=1.5, linestyle=pairwise_line_style,
                                zorder=10, alpha=0.6)
                            
                            if p_val_pair < 0.001:
                                asterisk = "***"
                            elif p_val_pair < 0.01:
                                asterisk = "**"
                            elif p_val_pair < 0.05:
                                asterisk = "*"
                            else:
                                asterisk = "~"
                            
                            ax.text((i + j) / 2, y_line, asterisk,
                                ha='center', va='bottom',
                                fontsize=int(sig_marker_fontsize * pairwise_sig_modifier),
                                fontweight='bold')
                
                # Show overall significance
                if not show_stats:
                    if kw_p < 0.001:
                        sig_text = "***"
                    elif kw_p < 0.01:
                        sig_text = "**"
                    elif kw_p < 0.05:
                        sig_text = "*"
                    elif kw_p < 0.1:
                        sig_text = "~"
                    else:
                        sig_text = ""
                    
                    ax.text(sig_marker_x, sig_marker_y, sig_text, transform=ax.transAxes,
                        verticalalignment='top', horizontalalignment='center',
                        fontsize=sig_marker_fontsize, fontweight='bold')
            
            # Set ylabel only if in the first column:
            if col == 0:
                # Wrap long y-axis labels so they don't overflow the figure
                import textwrap as _tw
                _ylabel_wrapped = "\n".join(_tw.wrap(str(ylabel), width=20))
                ax.set_ylabel(_ylabel_wrapped, fontsize=ylabel_fontsize, labelpad=10)
            
            # Otherwise, delete the label AND the axis itself! 
            else:
                ax.set_ylabel('')
                ax.yaxis.set_visible(False)
                ax.yaxis.set_ticks([])

            # Set xlabel only on the bottom row AND only when xlabel=True.
            # Pass xlabel=False for panels that sit above another panel in a composite
            # figure so the axis title doesn't repeat on every row.
            if row == nrows - 1 and xlabel:
                ax.set_xlabel(dv_to_lab_short.get(group_var, group_var), fontsize=xlabel_fontsize, labelpad=10)
            else:
                ax.set_xlabel('')
                # ax.xaxis.set_visible(False)
                # ax.xaxis.set_ticks([])

            # # X-axis label only if middle column
            # if col == ncols // 2:
            #     ax.set_xlabel(dv_to_lab_short.get(group_var, group_var), fontsize=title_fontsize)
            # else:
            #     ax.set_xlabel('')

            
            # Set xtick labels
            if xtick_labels:
                ax.set_xticklabels([_wrap_label(x) for x in xtick_labels], rotation=0, ha='center')
            else:
                ax.set_xticklabels([_wrap_label(x) for x in groups], rotation=0, ha='center')
            
            # Tick font sizes: explicit overrides take precedence over multipliers
            _xtick_fs = xtick_fontsize if xtick_fontsize is not None else (
                title_fontsize * xtick_font_multiplier if xtick_font_multiplier else None
            )
            if _xtick_fs is not None:
                ax.tick_params(axis='x', labelsize=_xtick_fs)

            _ytick_fs = ytick_fontsize if ytick_fontsize is not None else (
                title_fontsize * ytick_font_multiplier if ytick_font_multiplier else None
            )
            if _ytick_fs is not None:
                ax.tick_params(axis='y', labelsize=_ytick_fs)
            
            # Apply ylim if provided
            if ylim is not None:
                ax.set_ylim(ylim)
            
            # Remove borders
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.tick_params(axis="both", which="both", length=0)


                
            print(f"Panel {panel_idx} completed successfully!")
            
        except Exception as e:
            print(f"\n{'!'*60}")
            print(f"ERROR in panel {panel_idx} (Row {panel_idx // ncols}, Col {panel_idx % ncols})")
            print(f"DV: {spec.get('dv', 'N/A')}, Group: {spec.get('group_var', 'N/A')}")
            print(f"Error type: {type(e).__name__}")
            print(f"Error message: {str(e)}")
            print(f"{'!'*60}\n")
            
            # Mark this panel as errored
            ax.text(0.5, 0.5, f"Error:\n{dv}\n{str(e)[:50]}", 
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=10, color='red')
            continue
    
    # Hide unused subplots
    for panel_idx in range(n_panels, nrows * ncols):
        row = panel_idx // ncols
        col = panel_idx % ncols
        axes[row, col].axis('off')
    
    plt.tight_layout()
    
    # Save figure and associated test summaries
    if savepath:
        base_savepath = savepath
        savepath_dir = os.path.dirname(savepath)
        if savepath_dir:
            os.makedirs(savepath_dir, exist_ok=True)
        if not savepath.endswith('.png'):
            savepath = f"{savepath}.png"
        plt.savefig(savepath, format='png', dpi=dpi, bbox_inches='tight')
        # Also save vector SVG for true-vector figure assembly
        plt.savefig(savepath.replace('.png', '.svg'), format='svg', bbox_inches='tight')

        summary_rows = []
        for panel, panel_result in results.items():
            base_row = {
                "panel": panel,
                "dv": panel_result.get("dv"),
                "group_var": panel_result.get("group_var"),
                "test": panel_result.get("test", ""),
                "n_groups": panel_result.get("n_groups"),
                "sample_size": panel_result.get("sample_size"),
                "groups": "|".join(panel_result.get("groups", [])),
            }

            if panel_result.get("test") == "mann_whitney":
                row = base_row.copy()
                row.update(
                    {
                        "row_type": "overall",
                        "comparison": "",
                        "statistic": panel_result.get("statistic"),
                        "p_value": panel_result.get("p_value"),
                        "rank_biserial_correlation": panel_result.get("rank_biserial_correlation"),
                        "group_1_label": panel_result.get("group_1_label"),
                        "group_1_n": panel_result.get("group_1_n"),
                        "group_1_median": panel_result.get("group_1_median"),
                        "group_1_iqr_q1": panel_result.get("group_1_iqr_q1"),
                        "group_1_iqr_q3": panel_result.get("group_1_iqr_q3"),
                        "group_2_label": panel_result.get("group_2_label"),
                        "group_2_n": panel_result.get("group_2_n"),
                        "group_2_median": panel_result.get("group_2_median"),
                        "group_2_iqr_q1": panel_result.get("group_2_iqr_q1"),
                        "group_2_iqr_q3": panel_result.get("group_2_iqr_q3"),
                    }
                )
                summary_rows.append(row)
            elif panel_result.get("test") == "kruskal_wallis":
                row = base_row.copy()
                row.update(
                    {
                        "row_type": "overall",
                        "comparison": "",
                        "statistic": panel_result.get("H"),
                        "p_value": panel_result.get("p_value"),
                    }
                )
                summary_rows.append(row)

                for comparison, pairwise_result in panel_result.get("pairwise", {}).items():
                    pairwise_row = base_row.copy()
                    pairwise_row.update(
                        {
                            "row_type": "pairwise",
                            "comparison": comparison,
                            "group_1": pairwise_result.get("group_1"),
                            "group_2": pairwise_result.get("group_2"),
                            "statistic": pairwise_result.get("statistic"),
                            "p_value": pairwise_result.get("p_value"),
                        }
                    )
                    summary_rows.append(pairwise_row)

        summary_df = pd.DataFrame(summary_rows)
        summary_base = os.path.splitext(base_savepath)[0]
        summary_dir = os.path.join(os.path.dirname(summary_base), "summary_results")
        os.makedirs(summary_dir, exist_ok=True)
        summary_csv_path = os.path.join(summary_dir, f"{os.path.basename(summary_base)}.csv")
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"Saved summary results to: {summary_csv_path}")

    return results

__all__ = ["multipanel_boxplot_grid"]
