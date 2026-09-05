### figure_assembly.py
### Tight, publication-quality manuscript figure assembly using matplotlib.
### Panel labels are rendered in Arial Bold (fontfamily='Arial' in _add_panel_label).
###
### Usage:
###   from figure_assembly import assemble_manuscript_figure
###   path = assemble_manuscript_figure(config, base_dir="/abs/path/to/project/root")
###
### Config schema:
###   {
###     "title": "Figure 2 — ...",           # optional header above figure
###     "output_path": "figures/figure_2.tiff",
###     "fig_width_inches": 7.205,            # target width (183 mm = Sci Reports double column)
###     "dpi": 300,
###     "gap_inches": 0.05,                  # gap between panels (inches)
###     "row_gap_inches": 0.08,              # gap between rows (inches); defaults to gap_inches
###     "label_fontsize": 10,
###     "label_color": "black",
###     "rows": [
###       {
###         "panels": [
###           {"path": "results/...", "label": "a"},   # path relative to base_dir or absolute
###           {"path": "results/...", "label": "b"},
###         ]
###       },
###       ...
###     ]
###   }

import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from PIL import Image


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_aspect_ratio(path):
    """Return width/height aspect ratio for an image file."""
    with Image.open(str(path)) as img:
        w, h = img.size
    return w / h


def _resolve_path(p, base_dir):
    """Resolve a path string relative to base_dir if not absolute."""
    p = Path(p)
    if p.is_absolute():
        return p
    return (Path(base_dir) / p).resolve()


def _compute_layout(config, base_dir):
    """
    Compute the figure size and per-row layout data.

    Returns
    -------
    figsize : (width_in, height_in)
    row_data : list of dicts, one per row:
        {
          'row_height_in': float,
          'panels': [{'path': Path, 'label': str, 'ar': float}, ...],
          'width_ratios': [float, ...],
        }
    """
    fig_width = config.get('fig_width_inches', 7.0)
    gap = config.get('gap_inches', 0.05)

    row_data = []
    total_height = 0.0

    for row_cfg in config['rows']:
        panels_out = []
        ars = []
        for panel_cfg in row_cfg['panels']:
            path = _resolve_path(panel_cfg['path'], base_dir)
            if not path.exists():
                raise FileNotFoundError(f"Panel image not found: {path}")
            ar = _get_aspect_ratio(path)
            ars.append(ar)
            panel_entry = {
                'path': path,
                'label': panel_cfg.get('label', ''),
                'ar': ar,
            }
            # Preserve per-panel label position overrides so assemble_manuscript_figure
            # can apply them (without this, panel.get('label_x') always falls back to
            # the figure-level default and per-panel overrides are silently ignored).
            for _k in ('label_x', 'label_y'):
                if _k in panel_cfg:
                    panel_entry[_k] = panel_cfg[_k]
            panels_out.append(panel_entry)

        n = len(ars)
        # H * sum(ar) + (n-1)*gap = fig_width  →  H = (fig_width - (n-1)*gap) / sum(ar)
        H = (fig_width - (n - 1) * gap) / sum(ars)

        # Optional per-row height_scale: multiply H to add vertical breathing
        # room (e.g. 1.1 = 10% taller).  Panels keep their native aspect ratio
        # and are bottom-aligned within the taller row.
        height_scale = row_cfg.get('height_scale', 1.0)
        H *= height_scale

        row_data.append({
            'row_height_in': H,
            'height_scale': height_scale,
            'panels': panels_out,
            'width_ratios': ars,
        })
        total_height += H

    # Add row gaps between rows
    row_gap = config.get('row_gap_inches', gap)
    n_rows = len(row_data)
    total_height += (n_rows - 1) * row_gap

    # Add title space if requested
    title = config.get('title', '')
    title_height = 0.25 if title else 0.0
    total_height += title_height

    return (fig_width, total_height), row_data, title_height


def _add_panel_label(ax, label, fontsize, color, x_offset=0.01, y_offset=0.99):
    """Add a bold panel label (e.g., 'a') to the upper-left of an axes.

    x_offset : float
        x position in axes coordinates.  Values in [0, 1] land inside the panel;
        negative values land to the LEFT of the panel (use with label_left_margin
        in assemble_manuscript_figure so the label is not clipped by the figure edge).
        When x_offset < 0 the text is right-aligned so it hugs the panel edge.
    y_offset : float
        y position in axes coordinates (default 0.99 = just inside top edge).
        Values > 1.0 land above the panel in the row gap or top figure margin.
    """
    if not label:
        return
    ha = 'right' if x_offset < 0 else 'left'
    ax.text(
        x_offset, y_offset, label,
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight='bold',
        fontfamily='Arial',
        va='top',
        ha=ha,
        color=color,
        zorder=10,
        clip_on=False,
    )


# ---------------------------------------------------------------------------
# Font-size standardization
# ---------------------------------------------------------------------------

def compute_source_fontsize(target_pt, source_height, fig_width, gap, aspect_ratios):
    """Compute the source font size so text appears at *target_pt* in the
    assembled figure.

    The assembled figure arranges panels in a row.  Each panel shares a
    common display height *H* determined by the figure width, gaps, and
    the panels' aspect ratios.  The scale factor for a panel is ``H /
    source_height``, and the required source font size is ``target_pt /
    scale``.

    Parameters
    ----------
    target_pt : float
        Desired apparent font size in the final assembled figure (points).
    source_height : float
        Height of the source panel in inches (the *figsize* height used
        when generating the panel).
    fig_width : float
        Width of the assembled figure in inches (e.g. 7.0).
    gap : float
        Horizontal gap between adjacent panels in inches (e.g. 0.05).
    aspect_ratios : list of float
        Width-over-height aspect ratios of **all** panels in the same row,
        in display order (e.g. ``[6/15, 8/9]`` for a boxplot + correlation
        pair).

    Returns
    -------
    float
        Font size in points to use when generating the source panel so
        that it appears at *target_pt* after assembly scaling.
    """
    n = len(aspect_ratios)
    H = (fig_width - (n - 1) * gap) / sum(aspect_ratios)
    scale = H / source_height
    return target_pt / scale


# ---------------------------------------------------------------------------
# Raster → SVG conversion
# ---------------------------------------------------------------------------

def raster_to_svg(png_path, svg_path, dpi=600):
    """Wrap a raster image in a minimal SVG container.

    The SVG's physical dimensions (width/height in mm) are derived from
    the image's pixel dimensions and *dpi*, so the output prints at the
    correct physical size.  The raster data is base64-encoded inline —
    no external file references.

    Parameters
    ----------
    png_path : str or Path
        Source raster image (any format PIL can open).
    svg_path : str or Path
        Destination SVG file.
    dpi : int
        Resolution used to compute physical dimensions (pixels → mm).
    """
    import base64, io

    img = Image.open(str(png_path))
    w_px, h_px = img.size

    buf = io.BytesIO()
    img.convert('RGB').save(buf, format='PNG')
    b64_data = base64.b64encode(buf.getvalue()).decode('ascii')

    w_mm = w_px / dpi * 25.4
    h_mm = h_px / dpi * 25.4

    svg_content = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg"\n'
        '     xmlns:xlink="http://www.w3.org/1999/xlink"\n'
        f'     width="{w_mm:.2f}mm" height="{h_mm:.2f}mm"\n'
        f'     viewBox="0 0 {w_px} {h_px}">\n'
        f'  <image width="{w_px}" height="{h_px}"\n'
        f'         xlink:href="data:image/png;base64,{b64_data}"/>\n'
        '</svg>\n'
    )

    Path(svg_path).parent.mkdir(parents=True, exist_ok=True)
    with open(str(svg_path), 'w') as f:
        f.write(svg_content)


def scale_svg_to_dpi(svg_path, dpi=600, target_width_mm=None,
                     height_pad_frac=0.005, top_margin_mm=0.0):
    """Scale an SVG's display dimensions.

    When *target_width_mm* is given, the SVG artboard is set to that exact
    width (in mm) with height computed from the aspect ratio plus a small
    padding (*height_pad_frac*, default 0.5%) so the artboard extends just
    beyond the figure content.  Otherwise dimensions are scaled from points
    to *dpi*-equivalent pixels.

    *top_margin_mm* adds whitespace above the figure content by shifting the
    viewBox origin upward and increasing the artboard height accordingly.

    Modifies the file in-place.
    """
    from lxml import etree
    tree = etree.parse(str(svg_path))
    root = tree.getroot()

    # Read existing viewBox if present, otherwise use width/height attributes
    vb = root.get('viewBox')
    if vb:
        parts = vb.split()
        w_pt = float(parts[2]) - float(parts[0])
        h_pt = float(parts[3]) - float(parts[1])
    else:
        w_pt = _parse_svg_dimension(root.get('width', '0'))
        h_pt = _parse_svg_dimension(root.get('height', '0'))
    if w_pt <= 0 or h_pt <= 0:
        return

    # Convert top margin from mm to points (1mm = 2.8346pt)
    margin_pt = top_margin_mm * 2.8346
    # Shift viewBox origin upward so content renders lower on the artboard
    root.set('viewBox', f'0 {-margin_pt:.4f} {w_pt:.4f} {h_pt + margin_pt:.4f}')

    if target_width_mm is not None:
        aspect = (h_pt + margin_pt) / w_pt
        height_mm = target_width_mm * aspect * (1.0 + height_pad_frac)
        root.set('width', f'{target_width_mm}mm')
        root.set('height', f'{height_mm:.4f}mm')
    else:
        scale = dpi / 72.0
        root.set('width', f'{w_pt * scale:.2f}')
        root.set('height', f'{(h_pt + margin_pt) * scale * (1.0 + height_pad_frac):.2f}')

    tree.write(str(svg_path), xml_declaration=True, standalone='yes', encoding='ASCII')


# ---------------------------------------------------------------------------
# True-vector SVG assembly via svgutils
# ---------------------------------------------------------------------------

def _parse_svg_dimension(s):
    """Parse an SVG dimension string (e.g. '195.5pt', '200px') to float points."""
    s = str(s)
    for suffix in ('pt', 'px', 'in', 'mm', 'cm', 'em'):
        if s.endswith(suffix):
            val = float(s[:-len(suffix)])
            if suffix == 'in':
                return val * 72.0
            elif suffix == 'mm':
                return val * 72.0 / 25.4
            elif suffix == 'cm':
                return val * 72.0 / 2.54
            return val  # pt and px treated as 1:1
    return float(s)


def _get_svg_aspect_ratio(svg_path):
    """Return width/height aspect ratio from an SVG file."""
    import svgutils.transform as st
    svg = st.fromfile(str(svg_path))
    w = _parse_svg_dimension(svg.width)
    h = _parse_svg_dimension(svg.height)
    return w / h, w, h


def assemble_manuscript_figure_svg(config, output_path=None, base_dir=None):
    """Assemble a true-vector SVG figure from individual panel SVGs.

    Mirrors the layout of ``assemble_manuscript_figure()`` but keeps all
    elements as vector paths by compositing SVGs via *svgutils* instead of
    rasterizing through PIL.

    Panel paths in *config* should reference ``.png`` files — this function
    automatically looks for the ``.svg`` sibling.  If an SVG panel is missing,
    the function falls back to embedding the PNG as a raster element.

    Returns
    -------
    str
        Absolute path of the saved SVG file.
    """
    import svgutils.transform as st

    if base_dir is None:
        base_dir = os.getcwd()
    base_dir = str(Path(base_dir).resolve())

    out = output_path or config.get('output_path', 'figures/assembled/figure.svg')
    out = _resolve_path(out, base_dir)
    # Force .svg extension
    out = out.with_suffix('.svg')
    out.parent.mkdir(parents=True, exist_ok=True)

    fig_width_in = config.get('fig_width_inches', 7.0)
    gap_in = config.get('gap_inches', 0.05)
    row_gap_in = config.get('row_gap_inches', gap_in)
    label_fontsize = config.get('label_fontsize', 10)
    label_color = config.get('label_color', 'black')
    label_x_frac = config.get('label_x', 0.01)
    label_y_frac = config.get('label_y', 0.99)
    label_left_margin = config.get('label_left_margin', 0.0)

    # Work in points (1 inch = 72 pt)
    fig_width_pt = fig_width_in * 72.0
    gap_pt = gap_in * 72.0
    row_gap_pt = row_gap_in * 72.0
    content_left_pt = label_left_margin * fig_width_pt

    # Compute layout: for each row, determine panel aspect ratios and row height
    rows_layout = []
    total_height_pt = 0.0

    for row_cfg in config['rows']:
        panels_info = []
        ars = []
        for panel_cfg in row_cfg['panels']:
            png_path = _resolve_path(panel_cfg['path'], base_dir)
            svg_path = png_path.with_suffix('.svg')

            if svg_path.exists():
                ar, native_w, native_h = _get_svg_aspect_ratio(svg_path)
                panel_entry = {
                    'svg_path': svg_path,
                    'is_vector': True,
                    'native_w': native_w,
                    'native_h': native_h,
                    'ar': ar,
                    'label': panel_cfg.get('label', ''),
                    'label_x': panel_cfg.get('label_x', label_x_frac),
                    'label_y': panel_cfg.get('label_y', label_y_frac),
                }
            else:
                # Fallback: embed raster PNG
                ar = _get_aspect_ratio(png_path)
                with Image.open(str(png_path)) as img:
                    pw, ph = img.size
                panel_entry = {
                    'png_path': png_path,
                    'is_vector': False,
                    'native_w': pw,
                    'native_h': ph,
                    'ar': ar,
                    'label': panel_cfg.get('label', ''),
                    'label_x': panel_cfg.get('label_x', label_x_frac),
                    'label_y': panel_cfg.get('label_y', label_y_frac),
                }
                print(f"  [SVG assembly] no SVG for {png_path.name}, embedding raster")

            ars.append(ar)
            panels_info.append(panel_entry)

        n = len(ars)
        content_width = fig_width_pt - content_left_pt
        H_pt = (content_width - (n - 1) * gap_pt) / sum(ars)
        height_scale = row_cfg.get('height_scale', 1.0)
        H_pt *= height_scale

        rows_layout.append({
            'row_height_pt': H_pt,
            'height_scale': height_scale,
            'panels': panels_info,
            'aspect_ratios': ars,
        })
        total_height_pt += H_pt

    n_rows = len(rows_layout)
    total_height_pt += (n_rows - 1) * row_gap_pt

    # Create composed SVG
    composed = st.SVGFigure(f'{fig_width_pt}pt', f'{total_height_pt}pt')
    elements = []

    y_cursor = 0.0

    for row_idx, rd in enumerate(rows_layout):
        x_cursor = content_left_pt
        n_panels = len(rd['panels'])
        content_width = fig_width_pt - content_left_pt
        H_pt = rd['row_height_pt']

        for j, panel in enumerate(rd['panels']):
            panel_w_pt = panel['ar'] * (H_pt / rd['height_scale'])

            if panel['is_vector']:
                svg_panel = st.fromfile(str(panel['svg_path']))
                root = svg_panel.getroot()
                native_w = panel['native_w']
                native_h = panel['native_h']
                # Scale to fit: target height = H_pt / height_scale (actual content height)
                actual_h = H_pt / rd['height_scale']
                scale = actual_h / native_h
                root.scale(scale)
                # Bottom-align when height_scale > 1
                y_offset = H_pt - actual_h if rd['height_scale'] > 1.0 else 0.0
                root.moveto(x_cursor, y_cursor + y_offset)
                elements.append(root)
            else:
                # Embed raster as SVG image element
                import base64, io
                img = Image.open(str(panel['png_path']))
                buf = io.BytesIO()
                img.convert('RGB').save(buf, format='PNG')
                b64 = base64.b64encode(buf.getvalue()).decode('ascii')
                actual_h = H_pt / rd['height_scale']
                y_offset = H_pt - actual_h if rd['height_scale'] > 1.0 else 0.0
                img_elem = st.ImageElement(
                    str(panel['png_path']),
                    panel_w_pt, actual_h,
                )
                img_elem.moveto(x_cursor, y_cursor + y_offset)
                elements.append(img_elem)

            # Panel label
            if panel['label']:
                lx = panel['label_x']
                ly = panel['label_y']
                # Convert fractional label position to absolute coordinates
                label_x_pt = x_cursor + lx * panel_w_pt
                # label_y is from top (0.99 = near top)
                label_y_pt = y_cursor + (1.0 - ly) * (H_pt / rd['height_scale'])
                if rd['height_scale'] > 1.0:
                    label_y_pt += H_pt - H_pt / rd['height_scale']
                label_elem = st.TextElement(
                    label_x_pt, label_y_pt,
                    panel['label'],
                    size=label_fontsize,
                    weight='bold',
                    font='Arial',
                    color=label_color,
                )
                elements.append(label_elem)

            x_cursor += panel_w_pt + gap_pt

        y_cursor += H_pt + row_gap_pt

    composed.append(elements)
    composed.save(str(out))

    # svgutils.SVGFigure.save() does not write width/height/viewBox
    # attributes.  Inject them so downstream tools (scale_svg_to_dpi,
    # svgutils.fromfile for combining) can read the SVG dimensions.
    from lxml import etree as _etree_inject
    _inj_tree = _etree_inject.parse(str(out))
    _inj_root = _inj_tree.getroot()
    if not _inj_root.get('viewBox'):
        _inj_root.set('width', f'{fig_width_pt:.4f}pt')
        _inj_root.set('height', f'{y_cursor:.4f}pt')
        _inj_root.set('viewBox', f'0 0 {fig_width_pt:.4f} {y_cursor:.4f}')
        _inj_tree.write(str(out), xml_declaration=True,
                        standalone='yes', encoding='ASCII')

    # Set viewBox and scale display dimensions (DPI-based default; the final
    # figures section may override with target_width_mm later).
    scale_svg_to_dpi(str(out), dpi=config.get('dpi', 300))

    print(f"Saved vector SVG: {out}")
    return str(out)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def assemble_manuscript_figure(config, output_path=None, base_dir=None):
    """
    Assemble a tight, publication-quality manuscript figure.

    Parameters
    ----------
    config : dict
        Figure configuration dict. See module docstring for schema.
    output_path : str, optional
        Override config['output_path'].
    base_dir : str or Path, optional
        Directory against which relative paths in config are resolved.
        Defaults to the current working directory.

    Returns
    -------
    str
        Absolute path of the saved figure file.
    """
    if base_dir is None:
        base_dir = os.getcwd()
    base_dir = str(Path(base_dir).resolve())

    # Resolve output path
    out = output_path or config.get('output_path', 'figures/assembled/figure.png')
    out = _resolve_path(out, base_dir)
    out.parent.mkdir(parents=True, exist_ok=True)

    dpi = config.get('dpi', 300)
    gap = config.get('gap_inches', 0.05)
    row_gap = config.get('row_gap_inches', gap)
    label_fontsize = config.get('label_fontsize', 10)
    label_color = config.get('label_color', 'black')
    label_x = config.get('label_x', 0.01)
    label_y = config.get('label_y', 0.99)
    label_left_margin = config.get('label_left_margin', 0.0)  # fraction of figure width
    title = config.get('title', '')

    # Compute layout
    figsize, row_data, title_height = _compute_layout(config, base_dir)
    fig_width_in, fig_height_in = figsize

    fig = plt.figure(figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor('white')

    n_rows = len(row_data)
    row_heights_in = [rd['row_height_in'] for rd in row_data]

    # Build a top-level GridSpec with one row per figure row.
    # Use absolute heights so each row is exactly the computed height.
    # matplotlib height_ratios are proportional, so we use row_heights directly.
    # We account for the title by reserving space at the top via top margin.

    title_frac = title_height / fig_height_in  # fraction of figure height for title
    content_height = fig_height_in - title_height
    # row_gap between rows: expressed as a fraction of content height
    gap_frac_of_content = row_gap / content_height if content_height > 0 else 0

    gs_outer = gridspec.GridSpec(
        nrows=n_rows,
        ncols=1,
        figure=fig,
        top=1.0 - title_frac,
        bottom=0.0,
        left=label_left_margin,  # reserve left margin for outside-panel labels
        right=1.0,
        hspace=gap_frac_of_content * n_rows,  # hspace is relative to avg subplot height
        height_ratios=row_heights_in,
    )

    # Correct hspace: matplotlib hspace is fraction of avg subplot height.
    # avg subplot height (in axes coords) = (1 - title_frac) / n_rows (approx)
    # hspace in inches = row_gap
    # hspace (matplotlib) = row_gap / avg_subplot_height_in
    avg_subplot_h = sum(row_heights_in) / n_rows if n_rows > 0 else 1
    gs_outer.update(hspace=row_gap / avg_subplot_h)

    # Render each row
    for i, rd in enumerate(row_data):
        panels = rd['panels']
        n_panels = len(panels)
        width_ratios = rd['width_ratios']

        # Panel gap: fraction of the row width
        # wspace in matplotlib = gap / avg_subplot_width_in
        avg_panel_w = fig_width_in / n_panels
        wspace = gap / avg_panel_w if avg_panel_w > 0 else 0

        gs_row = gridspec.GridSpecFromSubplotSpec(
            1, n_panels,
            subplot_spec=gs_outer[i],
            width_ratios=width_ratios,
            wspace=wspace,
        )

        for j, panel in enumerate(panels):
            ax = fig.add_subplot(gs_row[j])
            img_arr = np.array(Image.open(str(panel['path'])).convert('RGBA'))
            ax.imshow(img_arr)
            ax.axis('off')
            # When height_scale > 1, bottom-align image so extra space is at top
            if rd.get('height_scale', 1.0) > 1.0:
                ax.set_anchor('S')
            # Per-panel label_x / label_y overrides fall back to figure-level defaults
            panel_lx = panel.get('label_x', label_x)
            panel_ly = panel.get('label_y', label_y)
            _add_panel_label(ax, panel['label'], label_fontsize, label_color, x_offset=panel_lx, y_offset=panel_ly)

    # Add title
    if title:
        fig.text(
            0.5, 1.0 - title_frac / 2,
            title,
            ha='center', va='center',
            fontsize=label_fontsize + 2,
            fontweight='bold',
        )

    # Determine format from extension
    ext = out.suffix.lower().lstrip('.')
    if ext in ('tif', 'tiff'):
        fmt = 'tiff'
    elif ext == 'pdf':
        fmt = 'pdf'
    elif ext == 'svg':
        fmt = 'svg'
    else:
        fmt = 'png'

    fig.savefig(str(out), dpi=dpi, format=fmt, bbox_inches='tight',
                facecolor='white', pad_inches=0.0)
    plt.close(fig)

    # matplotlib + bbox_inches='tight' often drops the PNG pHYs DPI chunk.
    # Re-save with PIL to embed DPI metadata so viewers/journals render at
    # the correct physical size.
    if fmt == 'png':
        _img = Image.open(str(out))
        _img.save(str(out), dpi=(dpi, dpi))

    print(f"Saved: {out}")
    return str(out)
