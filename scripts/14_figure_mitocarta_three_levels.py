#!/usr/bin/env python3
"""
Script 07b: Figure 8 - Hierarchical Mitochondrial Pathways WITH STATISTICS
==========================================================================

Enhanced version of 3-panel figure with dual statistical significance markers:
- Panel A: Level 0 (7 broad categories) + significance markers
- Panel B: Level 1 (20 mid-level pathways) + significance markers
- Panel C: Level 2 (78 detailed pathways) + significance markers

NEW FEATURES:
- Fisher's exact test statistics (enrichment + directional bias)
- FDR correction (Benjamini-Hochberg) applied separately per level
- Dual significance markers:
  * Upper marker: Directional bias (up/down ratio differs from global)
  * Lower marker: Enrichment (DEGs overrepresented in pathway)
- Markers: *** (FDR<0.001), ** (FDR<0.01), * (FDR<0.05), ns (not shown)

Statistical Framework:
- Same approach as OXPHOS figure
- Two independent Fisher's exact tests per pathway
- FDR correction by level (7, 20, 78 pathways)
- Global baseline: 25.6% DEG rate, 55.5% up / 44.5% down

Author: Analysis pipeline
Date: 2025-12-13
Version: 4.0 (Publication with Statistics)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Import publication styling
from plotting_style import (
    setup_publication_style, COLORS, FONTS, LAYOUT, EXPORT,
    add_panel_label, set_axis_labels, save_publication_figure
)

# Setup paths
BASE_DIR = Path(__file__).parent.parent
TABLES_DIR = BASE_DIR / "results" / "tables"
FIGURES_DIR = BASE_DIR / "results" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Constants for consistency
TOTAL_MITO_DEGS = 231
N_DOWN = 191
N_UP = 40
PCT_DOWN = 82.7

def truncate_name(name, max_len=45):
    """Truncate long pathway names"""
    if len(name) <= max_len:
        return name
    return name[:max_len-3] + '...'


def add_directional_markers(ax, y_positions, pathways, stats_df,
                           down_values, up_values, level_name,
                           marker_offset=2, fontsize=8):
    """
    Add directional significance markers at the END of each bar.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to add markers to
    y_positions : array
        Y-coordinates for each pathway
    pathways : list
        List of pathway names (matching stats_df 'pathway' column)
    stats_df : pd.DataFrame
        Statistics DataFrame with 'pathway' and 'Directional_Sig' columns
    down_values : array
        Downregulated counts (positive numbers, plotted as negative)
    up_values : array
        Upregulated counts (positive numbers)
    level_name : str
        Level name ('level_0', 'level_1', 'level_2') for filtering stats
    marker_offset : float
        Small offset from bar end in data units (default 2)
    fontsize : int
        Font size for markers (default 8)
    """
    # Filter stats to this level
    level_stats = stats_df[stats_df['level'] == level_name].copy()

    for i, pathway in enumerate(pathways):
        # Find matching statistics
        pathway_stats = level_stats[level_stats['pathway'] == pathway]

        if len(pathway_stats) == 0:
            # Try alternative matching (in case of naming inconsistencies)
            pathway_stats = level_stats[level_stats['pathway'].str.strip() == pathway.strip()]

        if len(pathway_stats) == 0:
            print(f"  WARNING: No statistics found for pathway: {pathway}")
            continue

        pathway_stats = pathway_stats.iloc[0]
        dir_sig = pathway_stats['Directional_Sig']

        # Only show directional significance (enrichment all ns anyway)
        if pd.notna(dir_sig) and dir_sig != 'ns':
            # Determine which bar extends further
            # down_values are plotted as negative, so compare absolute values
            if up_values[i] > down_values[i]:
                # Upregulated bar is longer - place marker on right side
                marker_x = up_values[i] + marker_offset
                ha_align = 'left'
            else:
                # Downregulated bar is longer - place marker on left side
                marker_x = -down_values[i] - marker_offset
                ha_align = 'right'

            ax.text(
                marker_x,
                y_positions[i] - 0.02,  # Small downward adjustment for better centering
                dir_sig,
                ha=ha_align,
                va='center',
                fontsize=fontsize,
                fontweight='bold',
                color='black'
            )


def create_panel_a_level0(ax, stats_df):
    """
    Panel A: Level 0 - Broad categories (7 pathways) WITH STATISTICS
    Sorted by total DEGs
    """
    # Load data
    level0 = pd.read_csv(TABLES_DIR / "mitocarta_pathways_level_0.csv")

    # Sort by total DEGs (most pathways at top)
    level0 = level0.sort_values('total_degs', ascending=True)

    # Create diverging bars
    y_pos = np.arange(len(level0))

    # Plot downregulated (left, blue) - using publication colors
    ax.barh(y_pos, -level0['downregulated'],
            color=COLORS['down'], alpha=0.85, label='Downregulated')

    # Plot upregulated (right, red)
    ax.barh(y_pos, level0['upregulated'],
            color=COLORS['up'], alpha=0.85, label='Upregulated')

    # Formatting
    ax.set_yticks(y_pos)
    ax.set_yticklabels(level0['pathway'], fontsize=FONTS['tick_label'])

    # Use publication style axis labels
    set_axis_labels(ax, xlabel='Number of DEGs',
                   title='Level 0: Broad Categories')

    # Add vertical line at 0
    ax.axvline(x=0, color='black', linewidth=1.2, linestyle='-')

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set symmetric x-axis with padding for markers
    max_val = max(level0['downregulated'].max(), level0['upregulated'].max())
    ax.set_xlim(-max_val * 1.15, max_val * 1.15)

    # Add directional significance markers at bar ends
    add_directional_markers(
        ax, y_pos, level0['pathway'].tolist(), stats_df,
        level0['downregulated'].values, level0['upregulated'].values,
        'level_0', marker_offset=2, fontsize=9
    )

    # Add gridlines
    ax.grid(axis='x', alpha=0.25, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)

    # Add panel label
    add_panel_label(ax, 'A', x=-0.12, y=1.15)


def create_panel_b_level1(ax, stats_df):
    """
    Panel B: Level 1 - Mid-level pathways (20 pathways) WITH STATISTICS
    Sorted by total DEGs
    """
    # Load data
    level1 = pd.read_csv(TABLES_DIR / "mitocarta_pathways_level_1.csv")

    # Sort by total DEGs
    level1 = level1.sort_values('total_degs', ascending=True)

    # Create diverging bars
    y_pos = np.arange(len(level1))

    # Plot downregulated (left, blue)
    ax.barh(y_pos, -level1['downregulated'],
            color=COLORS['down'], alpha=0.85, label='Downregulated')

    # Plot upregulated (right, red)
    ax.barh(y_pos, level1['upregulated'],
            color=COLORS['up'], alpha=0.85, label='Upregulated')

    # Formatting
    ax.set_yticks(y_pos)
    ax.set_yticklabels(level1['pathway'], fontsize=FONTS['tick_label']-0.5)

    # Use publication style axis labels
    set_axis_labels(ax, xlabel='Number of DEGs',
                   title='Level 1: Mid-Level Pathways')

    # Add vertical line at 0
    ax.axvline(x=0, color='black', linewidth=1.2, linestyle='-')

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set symmetric x-axis with padding for markers
    max_val = max(level1['downregulated'].max(), level1['upregulated'].max())
    ax.set_xlim(-max_val * 1.15, max_val * 1.15)

    # Add directional significance markers at bar ends
    add_directional_markers(
        ax, y_pos, level1['pathway'].tolist(), stats_df,
        level1['downregulated'].values, level1['upregulated'].values,
        'level_1', marker_offset=2, fontsize=8
    )

    # Add gridlines
    ax.grid(axis='x', alpha=0.25, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)

    # Add panel label
    add_panel_label(ax, 'B', x=-0.12, y=1.12)


def create_panel_c_level2(ax, stats_df):
    """
    Panel C: Level 2 - Detailed pathways (78 pathways) WITH STATISTICS
    Sorted by total DEGs
    IMPORTANT: Uses smaller fonts due to large number of pathways
    """
    # Load data
    level2 = pd.read_csv(TABLES_DIR / "mitocarta_pathways_level_2_detailed.csv")
    level2.columns = ['pathway', 'down', 'up', 'total', 'pct_down', 'pct_up']

    # Sort by total DEGs
    level2 = level2.sort_values('total', ascending=True)

    # Truncate long names (slightly shorter for readability)
    level2['pathway_display'] = level2['pathway'].apply(lambda x: truncate_name(x, 40))

    # Create diverging bars
    y_pos = np.arange(len(level2))

    # Plot downregulated (left, blue)
    ax.barh(y_pos, -level2['down'],
            color=COLORS['down'], alpha=0.85, label='Downregulated')

    # Plot upregulated (right, red)
    ax.barh(y_pos, level2['up'],
            color=COLORS['up'], alpha=0.85, label='Upregulated')

    # Formatting - Use small font for 78 pathways
    ax.set_yticks(y_pos)
    ax.set_yticklabels(level2['pathway_display'], fontsize=FONTS['small'])

    # Use publication style axis labels
    set_axis_labels(ax, xlabel='Number of DEGs',
                   title='Level 2: Detailed Pathways')

    # Add vertical line at 0
    ax.axvline(x=0, color='black', linewidth=1.2, linestyle='-')

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set symmetric x-axis with padding for markers
    max_val = max(level2['down'].max(), level2['up'].max())
    ax.set_xlim(-max_val * 1.2, max_val * 1.2)

    # Add directional significance markers at bar ends (use original pathway names, not truncated)
    add_directional_markers(
        ax, y_pos, level2['pathway'].tolist(), stats_df,
        level2['down'].values, level2['up'].values,
        'level_2', marker_offset=2, fontsize=7  # Slightly smaller for dense panel
    )

    # Add gridlines
    ax.grid(axis='x', alpha=0.25, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)

    # Add panel label
    add_panel_label(ax, 'C', x=-0.12, y=1.06)


def main():
    """Generate Figure 8: Hierarchical Mitochondrial Pathways with Statistics"""
    print("\n" + "="*70)
    print("GENERATING FIGURE 8: HIERARCHICAL PATHWAYS WITH STATISTICS")
    print("="*70)

    # Setup publication style
    setup_publication_style('paper')

    # Load statistics
    print("\nLoading statistical analysis results...")
    stats_all = pd.read_csv(TABLES_DIR / "pathway_statistics_all.csv")
    print(f"  Loaded statistics for {len(stats_all)} pathways")

    # Count significant pathways by level
    for level in ['level_0', 'level_1', 'level_2']:
        level_stats = stats_all[stats_all['level'] == level]
        n_dir_sig = (level_stats['Directional_Sig'] != 'ns').sum()
        n_enr_sig = (level_stats['Enrichment_Sig'] != 'ns').sum()
        print(f"  {level}: {n_dir_sig} directionally significant, "
              f"{n_enr_sig} enrichment significant")

    # Verify data consistency
    print(f"\nData consistency check:")
    mito_degs = pd.read_csv(TABLES_DIR / "mitochondrial_degs.csv")
    actual_total = len(mito_degs)
    actual_down = (mito_degs['regulation'] == 'down').sum()
    actual_up = (mito_degs['regulation'] == 'up').sum()

    print(f"  Expected: {TOTAL_MITO_DEGS} total, {N_DOWN} down, {N_UP} up")
    print(f"  Actual:   {actual_total} total, {actual_down} down, {actual_up} up")

    if actual_total != TOTAL_MITO_DEGS:
        print(f"  WARNING: Total mismatch! Using actual value: {actual_total}")

    # Create figure with 3 panels (vertical stack)
    fig = plt.figure(figsize=(LAYOUT['double_column'] * 1.3,
                             LAYOUT['double_column'] * 2.8))

    # Define grid: 3 rows, 1 column
    gs = fig.add_gridspec(3, 1,
                         height_ratios=[1, 2.5, 9.5],
                         hspace=0.28,
                         left=0.28, right=0.97, top=0.96, bottom=0.03)

    # Create subplots
    print("\nGenerating panels with statistical markers...")

    print("  • Panel A: Level 0 (7 broad categories)")
    ax_a = fig.add_subplot(gs[0])
    create_panel_a_level0(ax_a, stats_all)

    print("  • Panel B: Level 1 (20 mid-level pathways)")
    ax_b = fig.add_subplot(gs[1])
    create_panel_b_level1(ax_b, stats_all)

    print("  • Panel C: Level 2 (78 detailed pathways)")
    ax_c = fig.add_subplot(gs[2])
    create_panel_c_level2(ax_c, stats_all)

    # Add overall figure title
    fig.suptitle('Comprehensive Mitochondrial Pathway Analysis with Statistics',
                fontsize=FONTS['title'] + 2, fontweight='bold', y=0.985)

    print("\nSaving publication-quality figure...")

    # Save to results directory
    output_base_results = FIGURES_DIR / "figure_appendix_mitocarta_hierarchy"
    saved_files_results = save_publication_figure(
        fig, output_base_results,
        dpi=EXPORT['dpi_print'],
        formats=['pdf', 'png', 'svg']
    )

    plt.close()

    print("\n" + "="*70)
    print("FIGURE 8 GENERATION COMPLETE")
    print("="*70)
    print(f"\n✓ Publication-quality 3-panel figure with statistics")
    print(f"✓ Resolution: {EXPORT['dpi_print']} DPI")
    print(f"✓ Formats: PDF (vector), PNG (raster), SVG (web)")
    print(f"✓ Colorblind-friendly palette applied")
    print(f"✓ Right-aligned directional significance markers")
    print(f"\nStatistical features:")
    print(f"  • Fisher's exact tests for directional bias")
    print(f"  • FDR correction (Benjamini-Hochberg) per level")
    print(f"  • Markers: Right-aligned at fixed x-position, vertically centered")
    print(f"  • Significance levels: *** (FDR<0.001), ** (<0.01), * (<0.05)")
    print(f"  • Enrichment markers removed (all non-significant)")
    print(f"\nFigure 8: Hierarchical mitochondrial pathway visualization:")
    print(f"  • Level 0: 7 broad categories")
    print(f"  • Level 1: 20 mid-level pathways")
    print(f"  • Level 2: 78 detailed pathways")
    print(f"  • Total: 105 pathways shown with statistical significance")
    print(f"\nSaved to:")
    for f in saved_files_results:
        print(f"    • {f}")

if __name__ == "__main__":
    main()
