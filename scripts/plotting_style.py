#!/usr/bin/env python3
"""
Publication-Quality Figure Styling Configuration
Centralized style settings for all figures in the mitochondrial cachexia analysis

This module provides consistent, publication-ready styling that meets journal standards:
- Colorblind-friendly palettes
- Appropriate font sizes for print (legible at 50% scale)
- Vector-compatible rendering
- High-resolution export settings
- Consistent visual hierarchy

Usage:
    from plotting_style import setup_publication_style, COLORS, FONTS

    setup_publication_style()
    # Your plotting code here

Author: Analysis pipeline
Date: 2025-10-15
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
import seaborn as sns
import numpy as np

# ============================================================================
# COLOR SCHEMES - Colorblind-friendly and print-compatible
# ============================================================================

# Primary colors for up/downregulation
COLORS = {
    # Regulation colors (Paul Tol's colorblind-safe palette)
    'down': '#4477AA',  # Blue - downregulated
    'up': '#EE6677',    # Red/pink - upregulated
    'neutral': '#CCCCCC',  # Gray - not significant

    # Theme-specific colors
    'proteolysis': '#CC6677',  # Rose - for proteolysis theme
    'myogenesis': '#6699CC',   # Sky blue - for myogenesis theme
    'mitochondrial': '#AA3377', # Purple - for mitochondrial theme

    # Accent colors for additional categories
    'highlight': '#FFAA00',  # Orange - for highlighting
    'positive': '#228833',   # Green - for positive effects
    'negative': '#AA3377',   # Purple - for negative effects

    # Neutral/background
    'background': '#FFFFFF',
    'grid': '#E5E5E5',
    'text': '#000000',
    'axes': '#333333',
}

# Colorblind-friendly palettes for multi-category plots
PALETTES = {
    # Qualitative palette (for categories)
    'qualitative': ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'],

    # Sequential palette (for gradients)
    'sequential_blue': ['#F0F9FF', '#DEEBF7', '#C6DBEF', '#9ECAE1', '#6BAED6', '#4292C6', '#2171B5', '#084594'],
    'sequential_red': ['#FFF5F0', '#FEE0D2', '#FCBBA1', '#FC9272', '#FB6A4A', '#EF3B2C', '#CB181D', '#99000D'],

    # Diverging palette (for up/down comparisons)
    'diverging': ['#2166AC', '#4393C3', '#92C5DE', '#D1E5F0', '#F7F7F7', '#FDDBC7', '#F4A582', '#D6604D', '#B2182B'],
}

# ============================================================================
# FONT SETTINGS - Optimized for print and readability
# ============================================================================

FONTS = {
    # Base sizes (in points)
    'title': 14,           # Main figure title
    'panel_label': 16,     # Panel labels (A, B, C, etc.)
    'axis_label': 11,      # X and Y axis labels
    'tick_label': 9,       # Tick labels on axes
    'legend': 9,           # Legend text
    'annotation': 8,       # In-plot annotations
    'small': 7,            # Very small text (use sparingly)

    # Font family (standard sans-serif for science journals)
    'family': 'sans-serif',
    'sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],

    # Font weights
    'weight_normal': 'normal',
    'weight_bold': 'bold',
}

# ============================================================================
# LAYOUT SETTINGS - Dimensions and spacing
# ============================================================================

LAYOUT = {
    # Figure dimensions (inches) - standard journal sizes
    'single_column': 3.5,      # Single column width
    'double_column': 7.0,      # Double column width (full page width)
    'full_page': 10.0,         # Full page width (for large multi-panels)

    # Aspect ratios
    'aspect_square': 1.0,
    'aspect_wide': 1.5,        # Width 1.5x height
    'aspect_tall': 0.667,      # Height 1.5x width

    # Spacing and margins
    'margin': 0.5,             # Margin around figure (inches)
    'panel_spacing': 0.35,     # Space between panels
    'tight_spacing': 0.25,     # Tighter spacing for dense figures

    # Line and marker properties
    'linewidth': 1.5,
    'linewidth_thin': 1.0,
    'linewidth_thick': 2.0,
    'markersize': 6,
    'markersize_small': 4,
    'markersize_large': 8,
}

# ============================================================================
# EXPORT SETTINGS - High-resolution publication output
# ============================================================================

EXPORT = {
    # DPI settings
    'dpi_print': 600,          # For print submission
    'dpi_screen': 300,         # For screen viewing/presentations
    'dpi_draft': 150,          # For quick drafts

    # File formats
    'formats': {
        'vector': ['pdf', 'svg', 'eps'],  # Vector formats (scalable)
        'raster': ['png', 'tiff'],        # Raster formats (fixed resolution)
    },

    # PDF settings
    'pdf_fonttype': 42,        # TrueType fonts (editable in Illustrator)
    'ps_fonttype': 42,

    # Transparency
    'transparent': False,      # White background by default
    'savefig_pad_inches': 0.1, # Padding around saved figure
}

# ============================================================================
# MATPLOTLIB RCPARAMS - Global settings
# ============================================================================

def setup_publication_style(context='paper'):
    """
    Configure matplotlib with publication-quality settings

    Parameters
    ----------
    context : str
        'paper' - for manuscript figures (high quality)
        'poster' - for conference posters (larger fonts)
        'talk' - for presentations (even larger fonts)

    Returns
    -------
    None
        Modifies matplotlib rcParams globally
    """

    # Font scaling based on context
    scale_factors = {
        'paper': 1.0,
        'poster': 1.5,
        'talk': 1.8,
    }
    scale = scale_factors.get(context, 1.0)

    # Base rcParams dictionary
    rc_params = {
        # Figure settings
        'figure.facecolor': COLORS['background'],
        'figure.edgecolor': COLORS['background'],
        'figure.dpi': EXPORT['dpi_screen'],
        'savefig.dpi': EXPORT['dpi_print'],
        'savefig.format': 'pdf',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': EXPORT['savefig_pad_inches'],
        'savefig.transparent': EXPORT['transparent'],

        # Font settings
        'font.size': FONTS['tick_label'] * scale,
        'font.family': FONTS['family'],
        'font.sans-serif': FONTS['sans-serif'],

        # Axes settings
        'axes.labelsize': FONTS['axis_label'] * scale,
        'axes.titlesize': FONTS['title'] * scale,
        'axes.labelweight': FONTS['weight_normal'],
        'axes.titleweight': FONTS['weight_bold'],
        'axes.linewidth': LAYOUT['linewidth_thin'],
        'axes.edgecolor': COLORS['axes'],
        'axes.labelcolor': COLORS['text'],
        'axes.facecolor': COLORS['background'],
        'axes.grid': False,  # Turn off grid by default (enable selectively)
        'axes.spines.top': False,     # Remove top spine (cleaner look)
        'axes.spines.right': False,   # Remove right spine
        'axes.axisbelow': True,       # Put grid/ticks below data

        # Tick settings
        'xtick.labelsize': FONTS['tick_label'] * scale,
        'ytick.labelsize': FONTS['tick_label'] * scale,
        'xtick.major.width': LAYOUT['linewidth_thin'],
        'ytick.major.width': LAYOUT['linewidth_thin'],
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'xtick.minor.size': 2,
        'ytick.minor.size': 2,
        'xtick.color': COLORS['axes'],
        'ytick.color': COLORS['axes'],
        'xtick.direction': 'out',
        'ytick.direction': 'out',

        # Legend settings
        'legend.fontsize': FONTS['legend'] * scale,
        'legend.frameon': True,
        'legend.framealpha': 0.9,
        'legend.fancybox': False,
        'legend.edgecolor': COLORS['axes'],
        'legend.facecolor': COLORS['background'],

        # Line settings
        'lines.linewidth': LAYOUT['linewidth'],
        'lines.markersize': LAYOUT['markersize'] * scale,
        'lines.markeredgewidth': 0.5,

        # Grid settings (when enabled)
        'grid.color': COLORS['grid'],
        'grid.linestyle': '--',
        'grid.linewidth': 0.5,
        'grid.alpha': 0.3,

        # PDF/PS export settings
        'pdf.fonttype': EXPORT['pdf_fonttype'],
        'ps.fonttype': EXPORT['ps_fonttype'],

        # Math text
        'mathtext.default': 'regular',  # Don't italicize all math

        # Error bars
        'errorbar.capsize': 3,
    }

    # Update matplotlib rcParams
    plt.rcParams.update(rc_params)

    # Set seaborn style (combines well with above settings)
    sns.set_style('ticks', {
        'axes.grid': False,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

    print(f"✓ Publication style configured for '{context}' context")
    print(f"  - Font scaling: {scale}x")
    print(f"  - Export DPI: {EXPORT['dpi_print']}")
    print(f"  - Colors: Colorblind-friendly palette")


def reset_style():
    """Reset matplotlib to default settings"""
    mpl.rcdefaults()
    print("✓ Matplotlib style reset to defaults")


# ============================================================================
# HELPER FUNCTIONS - Common plotting utilities
# ============================================================================

def add_panel_label(ax, label, x=-0.1, y=1.05, fontsize=None, **kwargs):
    """
    Add bold panel label (A, B, C, etc.) to a subplot

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object to add label to
    label : str
        Label text (e.g., 'A', 'B', 'C')
    x, y : float
        Position in axes coordinates
    fontsize : int, optional
        Font size (default: FONTS['panel_label'])
    **kwargs : dict
        Additional text properties
    """
    if fontsize is None:
        fontsize = FONTS['panel_label']

    defaults = {
        'transform': ax.transAxes,
        'fontsize': fontsize,
        'fontweight': 'bold',
        'va': 'top',
        'ha': 'right',
    }
    defaults.update(kwargs)

    ax.text(x, y, label, **defaults)


def format_pvalue(p, threshold=0.001):
    """
    Format p-value for display

    Parameters
    ----------
    p : float
        P-value
    threshold : float
        Threshold for scientific notation

    Returns
    -------
    str
        Formatted p-value string
    """
    if p < threshold:
        return f'P < {threshold}'
    else:
        return f'P = {p:.3f}'


def add_significance_bracket(ax, x1, x2, y, pval, height=0.05, **kwargs):
    """
    Add significance bracket with p-value annotation

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to add bracket to
    x1, x2 : float
        X positions for bracket ends
    y : float
        Y position for bracket
    pval : float
        P-value to display
    height : float
        Height of bracket vertical lines
    **kwargs : dict
        Additional line properties
    """
    # Draw bracket
    ax.plot([x1, x1, x2, x2], [y, y+height, y+height, y],
            'k-', linewidth=1, **kwargs)

    # Add p-value or stars
    if pval < 0.001:
        text = '***'
    elif pval < 0.01:
        text = '**'
    elif pval < 0.05:
        text = '*'
    else:
        text = 'ns'

    ax.text((x1+x2)/2, y+height, text, ha='center', va='bottom',
            fontsize=FONTS['annotation'])


def set_axis_labels(ax, xlabel=None, ylabel=None, title=None,
                    xlabel_italic=False, ylabel_italic=False):
    """
    Set axis labels with optional gene name formatting

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object
    xlabel, ylabel : str, optional
        Axis labels
    title : str, optional
        Subplot title
    xlabel_italic, ylabel_italic : bool
        Whether to italicize labels (for gene names)
    """
    if xlabel:
        style = 'italic' if xlabel_italic else 'normal'
        ax.set_xlabel(xlabel, fontsize=FONTS['axis_label'],
                     style=style, fontweight='bold')

    if ylabel:
        style = 'italic' if ylabel_italic else 'normal'
        ax.set_ylabel(ylabel, fontsize=FONTS['axis_label'],
                     style=style, fontweight='bold')

    if title:
        ax.set_title(title, fontsize=FONTS['title'],
                    fontweight='bold', pad=10)


def save_publication_figure(fig, filepath, dpi=None, formats=None):
    """
    Save figure in multiple publication formats

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save
    filepath : str or Path
        Base filepath (without extension)
    dpi : int, optional
        Resolution (default: EXPORT['dpi_print'])
    formats : list, optional
        List of formats to save (default: ['pdf', 'png'])

    Returns
    -------
    list
        List of saved file paths
    """
    from pathlib import Path

    if dpi is None:
        dpi = EXPORT['dpi_print']

    if formats is None:
        formats = ['pdf', 'png']

    filepath = Path(filepath)
    saved_files = []

    for fmt in formats:
        output_path = filepath.with_suffix(f'.{fmt}')
        fig.savefig(output_path, dpi=dpi, bbox_inches='tight',
                   pad_inches=EXPORT['savefig_pad_inches'],
                   transparent=EXPORT['transparent'])
        saved_files.append(output_path)
        print(f"  ✓ Saved: {output_path}")

    return saved_files


def get_diverging_colors(n_colors, palette='diverging'):
    """
    Get n colors from a diverging colormap

    Parameters
    ----------
    n_colors : int
        Number of colors needed
    palette : str
        Palette name from PALETTES

    Returns
    -------
    list
        List of hex color codes
    """
    if palette in PALETTES:
        colors = PALETTES[palette]
        # Sample evenly from palette
        indices = np.linspace(0, len(colors)-1, n_colors, dtype=int)
        return [colors[i] for i in indices]
    else:
        # Fallback to matplotlib colormap
        cmap = plt.cm.get_cmap(palette)
        return [mpl.colors.rgb2hex(cmap(i/(n_colors-1))) for i in range(n_colors)]


def format_gene_name(gene_name, italics=True):
    """
    Format gene name for display (italics for genes)

    Parameters
    ----------
    gene_name : str
        Gene symbol
    italics : bool
        Whether to apply italic formatting

    Returns
    -------
    str
        Formatted gene name for matplotlib
    """
    if italics:
        return f'$\\it{{{gene_name}}}$'  # Use math mode for italics
    else:
        return gene_name


# ============================================================================
# MAIN - Test/demo when run directly
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("PUBLICATION STYLE MODULE - Configuration Test")
    print("="*70)

    # Setup style
    setup_publication_style('paper')

    # Create test figure
    fig, axes = plt.subplots(2, 2, figsize=(LAYOUT['double_column'],
                                             LAYOUT['double_column']))
    fig.suptitle('Publication Style Test Figure', fontweight='bold')

    # Test different plot types
    # Panel A: Bar plot
    ax = axes[0, 0]
    x = ['Category 1', 'Category 2', 'Category 3']
    y_down = [15, 25, 10]
    y_up = [5, 10, 3]
    ax.barh(x, y_down, color=COLORS['down'], label='Down', alpha=0.8)
    ax.barh(x, y_up, left=y_down, color=COLORS['up'], label='Up', alpha=0.8)
    set_axis_labels(ax, 'Number of genes', title='Panel A: Stacked bars')
    ax.legend()
    add_panel_label(ax, 'A')

    # Panel B: Scatter plot
    ax = axes[0, 1]
    np.random.seed(42)
    x = np.random.randn(100)
    y = np.random.randn(100)
    colors = [COLORS['down'] if x[i] < 0 else COLORS['up'] for i in range(100)]
    ax.scatter(x, y, c=colors, s=30, alpha=0.6, edgecolors='k', linewidth=0.5)
    set_axis_labels(ax, 'log2 fold change', 'Statistical significance',
                   title='Panel B: Scatter plot')
    add_panel_label(ax, 'B')

    # Panel C: Line plot
    ax = axes[1, 0]
    x = np.linspace(0, 10, 100)
    y1 = np.sin(x)
    y2 = np.cos(x)
    ax.plot(x, y1, color=COLORS['down'], label='Condition 1', linewidth=2)
    ax.plot(x, y2, color=COLORS['up'], label='Condition 2', linewidth=2)
    set_axis_labels(ax, 'Time (hours)', 'Expression level',
                   title='Panel C: Time series')
    ax.legend()
    ax.grid(True, alpha=0.3)
    add_panel_label(ax, 'C')

    # Panel D: Heatmap-like
    ax = axes[1, 1]
    data = np.random.randn(10, 5)
    im = ax.imshow(data, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
    ax.set_xticks(range(5))
    ax.set_xticklabels([f'Gene {i+1}' for i in range(5)], rotation=45)
    ax.set_yticks(range(10))
    ax.set_yticklabels([f'Sample {i+1}' for i in range(10)])
    set_axis_labels(ax, title='Panel D: Expression heatmap')
    plt.colorbar(im, ax=ax, label='Z-score')
    add_panel_label(ax, 'D')

    plt.tight_layout()

    # Save test figure
    from pathlib import Path
    test_dir = Path(__file__).parent.parent / "results" / "test_figures"
    test_dir.mkdir(exist_ok=True, parents=True)

    save_publication_figure(fig, test_dir / "publication_style_test",
                           formats=['pdf', 'png'])

    plt.close()

    print("\n" + "="*70)
    print("TEST COMPLETE")
    print("="*70)
    print(f"\nTest figure saved to: {test_dir}")
    print("\nStyle configuration:")
    print(f"  • Colors: {len(COLORS)} defined")
    print(f"  • Fonts: {len(FONTS)} settings")
    print(f"  • Export DPI: {EXPORT['dpi_print']}")
    print(f"  • Panel label size: {FONTS['panel_label']}pt")
    print("\nReady for use in figure generation scripts!")
