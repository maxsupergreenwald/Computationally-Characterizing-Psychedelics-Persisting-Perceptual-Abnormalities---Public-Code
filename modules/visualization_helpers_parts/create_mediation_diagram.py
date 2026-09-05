"""Helper function `create_mediation_diagram`."""

from master_config import dv_to_lab_short

def create_mediation_diagram(
    med_results,
    predictor_label=None,
    mediator_label=None,
    dv_label=None,

    # Box positioning
    predictor_pos=(0, 0.5),
    mediator_pos=(0.5, 1),
    dv_pos=(1, 0.5),
    box_padding=30,  
    box_text_width=20,

    # Path names in df
    a_path_name = "a_path",
    b_path_name = "b_path",
    c_prime_path_name = "c_prime_direct",
    prop_mediated_name = "prop_mediated",
    indirect_name = "indirect_ab",

    # Arrow styling
    a_path_color='black',
    b_path_color='black',
    c_prime_path_color='black',
    color_paths_by_sign=True,
    positive_path_color=None,
    negative_path_color=None,
    sig_alpha=0.9,
    nonsig_alpha=0.5,
    avoid_label_overlap=True,
    label_box_padding=0.03,
    arrow_length_scale=0.8,
    arrow_width=6,
    arrow_style='simple',
    arrow_head_width=6,
    arrow_head_length=6,
    arrow_mutation_scale=70,

    # Text styling
    box_fontsize=28,
    box_fontweight='bold',
    stat_fontsize=18,
    stat_fontweight='normal',
    stat_color='black',

    # Stat box positioning (relative to arrow midpoint)
    a_stat_offset=(-0.12, 0.06),
    b_stat_offset=(0.12, 0.06),
    c_prime_stat_offset=(0.0, -0.1),
    indirect_text_y_offset=-0.0,

    # Figure settings
    figsize=(10, 6),
    dpi=300,

    # Stat formatting
    show_estimate=True,
    show_ci=False,
    ci_label="95% CI",       # label string shown before the interval, e.g. "94% HDI"
    show_prob=False,
    decimal_places=3,
    show_significance=True,
    coef_label="Β",          # prefix for estimate (e.g. "Β" for standardised, "Δ" for counterfactual)
    sig_thresholds=[0.001, 0.01, 0.05,0.1],
    sig_symbols=['***', '**', '*', '~'],

    # Additional customization
    use_dv_to_lab=True,
    dv_to_lab_dict=None,
    savepath = None,

    # Indirect effect display mode
    # When True, replaces "Prop. Mediated: ..." with "Mediated Effect Prob: {pd:.1f}%".
    # Used for COUNTERFACTUAL diagrams where proportion mediated is less meaningful.
    show_indirect_pd=False,
    # When show_indirect_pd=True, use this float (0–1) directly as the pd value
    # instead of computing from the combined DataFrame.  Should be sourced from
    # mc_mediation_summary.csv → p_direction for the NIE row, which uses MC draws
    # and is the authoritative pd for the indirect effect.
    indirect_pd_override=None,
    # When show_indirect_pd=True, display this value (response-scale Δ) in the
    # centre label instead of the generic "Δ = change" legend.  Should be sourced
    # from mc_mediation_summary.csv → median for the NIE row.
    nie_delta=None,
    # When True, multiply nie_delta by 100 and append "%" (use for binary DVs
    # where the counterfactual is on the probability scale, e.g. hppd_binary).
    nie_delta_as_pct=False,
):
    """
    Create a customizable mediation diagram.
    
    Parameters:
    -----------
    med_results : pd.DataFrame
        Mediation results dataframe with columns: predictor, mediator, dv, effect, 
        estimate, lower_95, upper_95, prob_above_0, prob_below_0
    predictor_label : str, optional
        Custom label for predictor box (if None, uses predictor from med_results)
    mediator_label : str, optional
        Custom label for mediator box (if None, uses mediator from med_results)
    dv_label : str, optional
        Custom label for DV box (if None, uses dv from med_results)
    predictor_pos : tuple
        (x, y) position for predictor box
    mediator_pos : tuple
        (x, y) position for mediator box
    dv_pos : tuple
        (x, y) position for DV box
    a_path_color, b_path_color, c_prime_path_color : str
        Colors for each path arrow
    color_paths_by_sign : bool
        If True, overrides path colors by effect sign
    positive_path_color, negative_path_color : str
        Colors used when `color_paths_by_sign=True`
    sig_alpha, nonsig_alpha : float
        Arrow alpha values for significant vs non-significant paths
    avoid_label_overlap : bool
        If True, clip arrow endpoints to text-label bounding boxes
    label_box_padding : float
        Extra padding (data units) around label boxes for arrow clipping
    arrow_length_scale : float
        Scale factor for arrow length in data coordinates (0.8 = 20% shorter)
    arrow_width, arrow_width, arrow_width : float
        Line widths for each path arrow
    arrow_style : str
        Matplotlib arrow style (e.g., '-|>', '->', '<->', 'fancy')
    arrow_head_width, arrow_head_length : float
        Arrow head dimensions
    box_fontsize : int
        Font size for variable names in boxes
    box_fontweight : str
        Font weight for variable names
    stat_fontsize : int
        Font size for statistics
    stat_fontweight : str
        Font weight for statistics
    stat_color : str
        Color for statistics text
    a_stat_offset, b_stat_offset, c_prime_stat_offset : tuple
        (x, y) offsets for stat text relative to arrow midpoint
    figsize : tuple
        Figure size (width, height)
    dpi : int
        Figure resolution
    show_estimate : bool
        Whether to show point estimate
    show_ci : bool
        Whether to show 95% CI
    show_prob : bool
        Whether to show probability above/below 0
    decimal_places : int
        Number of decimal places for statistics
    use_dv_to_lab : bool
        Whether to use dv_to_lab_short dictionary for labels
    dv_to_lab_dict : dict, optional
        Custom dictionary for variable name translation (if None, uses dv_to_lab_short)
    show_indirect_pd : bool
        When True, replaces "Prop. Mediated: ..." with "Mediated Effect Prob: {pd:.1f}%".
        Use for COUNTERFACTUAL diagrams where proportion mediated is less meaningful.
    indirect_pd_override : float or None
        When show_indirect_pd=True and this is not None, use this value (0–1 scale)
        directly as the pd.  Should come from mc_mediation_summary.csv → p_direction
        for the NIE row (MC-integration-based pd, more accurate than path posteriors).

    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyArrowPatch

    # Default path colors from current binary_palette (if not explicitly provided)
    if positive_path_color is None or negative_path_color is None:
        try:
            from master_config import binary_palette as _binary_palette
        except Exception:
            try:
                from modules.variable_labels import binary_palette as _binary_palette
            except Exception:
                _binary_palette = ['#faeee2', '#bc511c']
        negative_path_color = _binary_palette[0]
        positive_path_color = _binary_palette[1]
    
    # Get the dv_to_lab dictionary
    if use_dv_to_lab:
        if dv_to_lab_dict is None:
            # Assume dv_to_lab_short exists in global scope
            try:
                lab_dict = dv_to_lab_short
            except NameError:
                print("Warning: dv_to_lab_short not found, using raw variable names")
                lab_dict = {}
        else:
            lab_dict = dv_to_lab_dict
    else:
        lab_dict = {}
    
    lab_dict_nounderscore = {k.replace('_', ''): v for k, v in lab_dict.items()}
    
    ##############################################################################################################
    # Get variable names from results
    predictor_var = med_results['predictor'].iloc[0]
    mediator_var = med_results['mediator'].iloc[0]
    dv_var = med_results['dv'].iloc[0]
    
    # Set labels (use custom labels if provided, otherwise use dv_to_lab or raw names)
    def get_label(var_name, custom_label):
        if custom_label is not None:
            return custom_label
        # First try with original name
        if var_name in lab_dict:
            return lab_dict[var_name]
        # If not found, try dict that does not have underscores
        elif var_name in lab_dict_nounderscore:
            return lab_dict_nounderscore[var_name]
        # Fall back to raw variable name
        else:
            return var_name
    
    predictor_label = get_label(predictor_var, predictor_label)
    mediator_label = get_label(mediator_var, mediator_label)
    dv_label = get_label(dv_var, dv_label)
    ##############################################################################################################

    # Extract path statistics
    a_path = med_results[med_results['effect'] == a_path_name].iloc[0]
    b_path = med_results[med_results['effect'] == b_path_name].iloc[0]
    c_prime = med_results[med_results['effect'] == c_prime_path_name].iloc[0]
    
    ##############################################################################################################
    # Format statistics strings
    def _get_significance_symbol(row):
        if row['estimate'] > 0:
            p_val = 1 - row['prob_above_0']
        else:
            p_val = 1 - row['prob_below_0']

        for threshold, symbol in zip(sig_thresholds, sig_symbols):
            if p_val < threshold:
                return symbol
        return ''

    def _is_significant(row):
        # Treat only starred symbols as significant; "~" remains non-significant.
        return '*' in _get_significance_symbol(row)

    def format_stat(row):
        parts = []
        if show_estimate:
            estimate_text = f"{coef_label} = {row['estimate']:.{decimal_places}f}"
            
            # Add significance symbol if requested
            if show_significance:
                sig_symbol = _get_significance_symbol(row)
                
                if sig_symbol:
                    estimate_text += f"$^{{{sig_symbol}}}$"
            
            parts.append(estimate_text)
        
        if show_ci:
            parts.append(f"{ci_label} [{row['lower_95']:.{decimal_places}f}, {row['upper_95']:.{decimal_places}f}]")
        if show_prob:
            if row['estimate'] > 0:
                parts.append(f"P(>0) = {row['prob_above_0']:.{decimal_places}f}")
            else:
                parts.append(f"P(<0) = {row['prob_below_0']:.{decimal_places}f}")
        return '\n'.join(parts)
    
    a_stat_text = format_stat(a_path)
    b_stat_text = format_stat(b_path)
    c_prime_stat_text = format_stat(c_prime)



    ##############################################################################################################    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.3)
    ax.axis('off')

    # Draw boxes (just text, no borders)
    from textwrap import fill
    
    # Function to wrap text based on character limit
    def wrap_label(text, width=box_text_width):
        """Wrap text to specified character width"""
        return fill(text, width=width)
    
    predictor_text_artist = ax.text(
        predictor_pos[0],
        predictor_pos[1],
        wrap_label(predictor_label),
        ha='center',
        va='center',
        fontsize=box_fontsize,
        fontweight=box_fontweight,
    )
    mediator_text_artist = ax.text(
        mediator_pos[0],
        mediator_pos[1],
        wrap_label(mediator_label),
        ha='center',
        va='center',
        fontsize=box_fontsize,
        fontweight=box_fontweight,
    )
    dv_text_artist = ax.text(
        dv_pos[0],
        dv_pos[1],
        wrap_label(dv_label),
        ha='center',
        va='center',
        fontsize=box_fontsize,
        fontweight=box_fontweight,
    )

        # Add indirect effect / proportion mediated text in center
    try:
        if show_indirect_pd:
            # COUNTERFACTUAL mode: show pd (probability of direction) for indirect effect.
            # Prefer indirect_pd_override (from mc_mediation_summary.csv p_direction —
            # MC-integration-based, authoritative) over path-level posterior computation.
            if indirect_pd_override is not None:
                pd_val = indirect_pd_override * 100
            else:
                indirect_effect = med_results[med_results['effect'] == indirect_name].iloc[0]
                pd_val = max(indirect_effect['prob_above_0'], indirect_effect['prob_below_0']) * 100
            indirect_text = f"P(Δ≠0) = {pd_val:.1f}%"
        else:
            prop_mediated = med_results[med_results['effect'] == prop_mediated_name].iloc[0]

            # Format estimate with significance
            estimate_text = f"{prop_mediated['estimate']:.{decimal_places}f}"
            if show_significance:
                if prop_mediated['estimate'] > 0:
                    p_val = 1 - prop_mediated['prob_above_0']
                else:
                    p_val = 1 - prop_mediated['prob_below_0']

                sig_symbol = ''
                for threshold, symbol in zip(sig_thresholds, sig_symbols):
                    if p_val < threshold:
                        sig_symbol = symbol
                        break

                if sig_symbol:
                    estimate_text += f"$^{{{sig_symbol}}}$"

            # Build the text
            indirect_text = f"Prop. Mediated: {estimate_text} \n[{prop_mediated['lower_95']:.{decimal_places}f}, {prop_mediated['upper_95']:.{decimal_places}f}]"

            # Add probability if requested
            if show_prob:
                indirect_effect = med_results[med_results['effect'] == indirect_name].iloc[0]
                if indirect_effect['estimate'] > 0:
                    prob_text = f"\nP(>0) = {indirect_effect['prob_above_0']:.{decimal_places}f}"
                else:
                    prob_text = f"\nP(<0) = {indirect_effect['prob_below_0']:.{decimal_places}f}"
                indirect_text += prob_text

        # Calculate center position (between all three boxes)
        center_x = (predictor_pos[0] + mediator_pos[0] + dv_pos[0]) / 3
        center_y = (predictor_pos[1] + mediator_pos[1] + dv_pos[1]) / 3

        # Draw the text (bold)
        ax.text(center_x, (center_y + indirect_text_y_offset), indirect_text,
                ha='center', va='center',
                fontsize=stat_fontsize,
                fontweight='bold',
                color=stat_color)

        # In counterfactual diagrams, show the NIE Δ above the P(Δ≠0) text.
        # If nie_delta is provided (from mc_mediation_summary.csv), display the
        # actual median response-scale indirect effect; otherwise fall back to the
        # generic legend "Δ = change".
        if show_indirect_pd:
            if nie_delta is not None:
                if nie_delta_as_pct:
                    nie_label = f"\u0394 = {nie_delta * 100:+.1f}%"
                else:
                    nie_label = f"\u0394 = {nie_delta:+.3f}"
            else:
                nie_label = "\u0394 = change"
            ax.text(
                center_x,
                center_y + indirect_text_y_offset + 0.05,
                nie_label,
                ha='center',
                va='bottom',
                fontsize=stat_fontsize,
                fontstyle='italic',
                fontweight='bold',
                color=stat_color,
            )
    except (IndexError, KeyError):
        # If required row not found, skip
        pass
    
    # Resolve arrow colors from path sign if requested
    if color_paths_by_sign:
        def _sign_color(row):
            return positive_path_color if row['estimate'] >= 0 else negative_path_color

        a_path_color = _sign_color(a_path)
        b_path_color = _sign_color(b_path)
        c_prime_path_color = _sign_color(c_prime)

    a_alpha = sig_alpha if _is_significant(a_path) else nonsig_alpha
    b_alpha = sig_alpha if _is_significant(b_path) else nonsig_alpha
    c_alpha = sig_alpha if _is_significant(c_prime) else nonsig_alpha

    def _expanded_bbox_data(text_artist, pad_data):
        renderer = fig.canvas.get_renderer()
        bbox_disp = text_artist.get_window_extent(renderer=renderer)
        bbox_data = bbox_disp.transformed(ax.transData.inverted())
        return (
            bbox_data.x0 - pad_data,
            bbox_data.x1 + pad_data,
            bbox_data.y0 - pad_data,
            bbox_data.y1 + pad_data,
        )

    def _ray_exit_rect(center, direction, rect):
        x, y = center
        dx, dy = direction
        xmin, xmax, ymin, ymax = rect

        t_candidates = []
        if abs(dx) > 1e-12:
            t_candidates.append((xmax - x) / dx if dx > 0 else (xmin - x) / dx)
        if abs(dy) > 1e-12:
            t_candidates.append((ymax - y) / dy if dy > 0 else (ymin - y) / dy)

        t_positive = [t for t in t_candidates if t > 0]
        if not t_positive:
            return center

        t_exit = min(t_positive)
        return (x + t_exit * dx, y + t_exit * dy)

    def _clip_segment(start, end, start_rect=None, end_rect=None):
        sx, sy = start
        ex, ey = end
        dx, dy = ex - sx, ey - sy
        seg_len = (dx ** 2 + dy ** 2) ** 0.5
        if seg_len < 1e-12:
            return start, end

        ux, uy = dx / seg_len, dy / seg_len

        clipped_start = _ray_exit_rect(start, (dx, dy), start_rect) if start_rect else start
        clipped_end = _ray_exit_rect(end, (-dx, -dy), end_rect) if end_rect else end

        clipped_start = (clipped_start[0] + ux * label_box_padding, clipped_start[1] + uy * label_box_padding)
        clipped_end = (clipped_end[0] - ux * label_box_padding, clipped_end[1] - uy * label_box_padding)

        return clipped_start, clipped_end

    def _centered_segment(start, end, target_len):
        sx, sy = start
        ex, ey = end
        dx, dy = ex - sx, ey - sy
        seg_len = (dx ** 2 + dy ** 2) ** 0.5
        if seg_len < 1e-12:
            return start, end
        ux, uy = dx / seg_len, dy / seg_len
        mx, my = (sx + ex) / 2, (sy + ey) / 2
        hx, hy = (target_len * ux) / 2, (target_len * uy) / 2
        return (mx - hx, my - hy), (mx + hx, my + hy)

    predictor_rect = mediator_rect = dv_rect = None
    if avoid_label_overlap:
        fig.canvas.draw()
        predictor_rect = _expanded_bbox_data(predictor_text_artist, label_box_padding)
        mediator_rect = _expanded_bbox_data(mediator_text_artist, label_box_padding)
        dv_rect = _expanded_bbox_data(dv_text_artist, label_box_padding)

    predictor_to_mediator_clip = _clip_segment(
        predictor_pos,
        mediator_pos,
        predictor_rect,
        mediator_rect,
    )
    mediator_to_dv_clip = _clip_segment(
        mediator_pos,
        dv_pos,
        mediator_rect,
        dv_rect,
    )
    predictor_to_dv_clip = _clip_segment(
        predictor_pos,
        dv_pos,
        predictor_rect,
        dv_rect,
    )

    def _seg_len(seg):
        (sx, sy), (ex, ey) = seg
        return ((ex - sx) ** 2 + (ey - sy) ** 2) ** 0.5

    lengths = [
        _seg_len(predictor_to_mediator_clip),
        _seg_len(mediator_to_dv_clip),
        _seg_len(predictor_to_dv_clip),
    ]
    min_len = min([l for l in lengths if l > 1e-12]) if lengths else 0
    target_len = min_len * arrow_length_scale

    predictor_to_mediator = _centered_segment(
        predictor_to_mediator_clip[0],
        predictor_to_mediator_clip[1],
        target_len,
    )
    mediator_to_dv = _centered_segment(
        mediator_to_dv_clip[0],
        mediator_to_dv_clip[1],
        target_len,
    )
    predictor_to_dv = _centered_segment(
        predictor_to_dv_clip[0],
        predictor_to_dv_clip[1],
        target_len,
    )
    # Re-center the c-path arrow horizontally at mediator_pos[0]
    _c_mid_x = (predictor_to_dv[0][0] + predictor_to_dv[1][0]) / 2
    _c_x_shift = mediator_pos[0] - _c_mid_x
    predictor_to_dv = (
        (predictor_to_dv[0][0] + _c_x_shift, predictor_to_dv[0][1]),
        (predictor_to_dv[1][0] + _c_x_shift, predictor_to_dv[1][1]),
    )

    arrow_shrink = 0 if avoid_label_overlap else box_padding

    # Draw arrows
    # a path: predictor → mediator
    arrow_a = FancyArrowPatch(
        predictor_to_mediator[0], predictor_to_mediator[1],
        arrowstyle=arrow_style,
        mutation_scale=arrow_mutation_scale,
        linewidth=0,
        facecolor=a_path_color,
        edgecolor='none',
        alpha=a_alpha,
        shrinkA=arrow_shrink, shrinkB=arrow_shrink
    )
    ax.add_patch(arrow_a)
    
    # b path: mediator → dv
    arrow_b = FancyArrowPatch(
        mediator_to_dv[0], mediator_to_dv[1],
        arrowstyle=arrow_style,
        mutation_scale=arrow_mutation_scale,
        linewidth=0,
        facecolor=b_path_color,
        edgecolor='none',
        alpha=b_alpha,
        shrinkA=arrow_shrink, shrinkB=arrow_shrink
    )
    ax.add_patch(arrow_b)
    
    # c' path: predictor → dv (direct effect)
    arrow_c = FancyArrowPatch(
        predictor_to_dv[0], predictor_to_dv[1],
        arrowstyle=arrow_style,
        mutation_scale=arrow_mutation_scale,
        linewidth=0,
        facecolor=c_prime_path_color,
        edgecolor='none',
        alpha=c_alpha,
        shrinkA=arrow_shrink, shrinkB=arrow_shrink
    )
    ax.add_patch(arrow_c)
    
    ##############################################################################################################
    # Add statistics text
    # a path stats (between predictor and mediator)
    a_mid_x = (predictor_pos[0] + mediator_pos[0]) / 2 + a_stat_offset[0]
    a_mid_y = (predictor_pos[1] + mediator_pos[1]) / 2 + a_stat_offset[1]
    ax.text(a_mid_x, a_mid_y, a_stat_text,
            ha='center', va='bottom', fontsize=stat_fontsize, 
            fontweight=stat_fontweight, color=stat_color)
    
    # b path stats (between mediator and dv)
    b_mid_x = (mediator_pos[0] + dv_pos[0]) / 2 + b_stat_offset[0]
    b_mid_y = (mediator_pos[1] + dv_pos[1]) / 2 + b_stat_offset[1]
    ax.text(b_mid_x, b_mid_y, b_stat_text,
            ha='center', va='bottom', fontsize=stat_fontsize,
            fontweight=stat_fontweight, color=stat_color)
    
    # c' path stats (between predictor and dv) — nudge slightly upward to sit between arrows
    c_mid_x = (predictor_pos[0] + dv_pos[0]) / 2 + c_prime_stat_offset[0]
    c_mid_y = (predictor_pos[1] + dv_pos[1]) / 2 + c_prime_stat_offset[1] + 0.02
    ax.text(c_mid_x, c_mid_y, c_prime_stat_text,
            ha='center', va='top', fontsize=stat_fontsize,
            fontweight=stat_fontweight, color=stat_color)

    ##############################################################################################################

    
    plt.tight_layout()

    # Save if specified
    if savepath:
        if not savepath.endswith('.png'):
            savepath = f"{savepath}.png"
        plt.savefig(savepath, format='png', dpi=600, bbox_inches='tight')
        # Also save vector SVG for true-vector figure assembly
        plt.savefig(savepath.replace('.png', '.svg'), format='svg', bbox_inches='tight')

    return fig, ax

__all__ = ["create_mediation_diagram"]
