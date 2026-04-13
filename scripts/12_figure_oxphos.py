#!/usr/bin/env python3
"""
OXPHOS Hierarchical Figure Generation
======================================

Creates a publication-quality hierarchical bar plot showing MitoCarta OXPHOS
pathway regulation across three hierarchical levels:

- Level 0: OXPHOS parent (all OXPHOS genes)
- Level 1: Major categories (OXPHOS subunits/assembly, Complexes I-V)
- Level 2: Subcategories (CI/CIV subunits, assembly factors, etc.)

The figure shows:
- Stacked horizontal bars (downregulated/upregulated)
- Visual hierarchy with different bar heights and transparency
- Dual significance markers (directional and enrichment)
- Percentage labels within bars
- DEG counts outside bars

Author: Bioinformatics Team
Date: December 2025
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import numpy as np

# Import publication styling
from plotting_style import (
    setup_publication_style, COLORS, FONTS,
    save_publication_figure, add_panel_label
)


def filter_core_oxphos_categories(stats_df):
    """
    Filter statistics to core OXPHOS categories only.

    Excludes auxiliary pathways like apoptosis, lipid metabolism, etc.
    Keeps only electron transport chain complexes and their components.

    Parameters
    ----------
    stats_df : pd.DataFrame
        Full OXPHOS statistics

    Returns
    -------
    pd.DataFrame
        Filtered statistics with core OXPHOS categories only
    """
    # Define core OXPHOS categories to include
    core_categories = {
        # Level 0 - always include parent
        'level_0': ['OXPHOS'],

        # Level 1 - core ETC complexes and assembly
        'level_1': [
            'OXPHOS subunits',
            'OXPHOS assembly factors',
            'Complex I',
            'Complex II',
            'Complex III',
            'Complex IV',
            'Complex V',
        ],

        # Level 2 - subunits and assembly factors for each complex
        'level_2': [
            'CI subunits',
            'CI assembly factors',
            'CII subunits',
            'CII assembly factors',
            'CIII subunits',
            'CIII assembly factors',
            'CIV subunits',
            'CIV assembly factors',
            'CV subunits',
            'CV assembly factors',
        ]
    }

    # Filter to core categories
    filtered_rows = []
    for level_name, categories in core_categories.items():
        level_data = stats_df[
            (stats_df['level'] == level_name) &
            (stats_df['category'].isin(categories))
        ]
        filtered_rows.append(level_data)

    filtered_df = pd.concat(filtered_rows, ignore_index=True)

    return filtered_df


def build_hierarchical_data_structure(stats_df):
    """
    Build hierarchical data structure for plotting.

    Organizes categories into parent > L1 branches > L2 children
    with appropriate separators and visual hierarchy.

    Parameters
    ----------
    stats_df : pd.DataFrame
        Filtered OXPHOS statistics

    Returns
    -------
    list
        List of dictionaries with plotting data and metadata
    """
    panel_data = []

    # Level 0: OXPHOS parent
    parent = stats_df[stats_df['level'] == 'level_0'].iloc[0]
    panel_data.append({
        'label': 'OXPHOS (all)',
        'type': 'parent',
        'level': 0,
        'total_genes': parent['Total_Genes'],
        'degs': parent['DEGs'],
        'up': parent['Upregulated'],
        'down': parent['Downregulated'],
        'pct_up': parent['Pct_Up'],
        'pct_down': parent['Pct_Down'],
        'dir_sig': parent['Directional_Sig'],
        'enr_sig': parent['Enrichment_Sig']
    })

    panel_data.append({'type': 'separator'})

    # Get L1 and L2 data
    l1_data = stats_df[stats_df['level'] == 'level_1'].copy()
    l2_data = stats_df[stats_df['level'] == 'level_2'].copy()

    # Define L1 category order and their L2 children
    l1_order = [
        ('OXPHOS subunits', []),
        ('OXPHOS assembly factors', []),
        ('Complex I', ['CI subunits', 'CI assembly factors']),
        ('Complex IV', ['CIV subunits', 'CIV assembly factors']),
        ('Complex V', ['CV subunits', 'CV assembly factors']),
        ('Complex III', ['CIII subunits', 'CIII assembly factors']),
        ('Complex II', ['CII subunits', 'CII assembly factors']),
    ]

    # Build structure for each L1 category
    for l1_name, l2_children_names in l1_order:
        # Get L1 data
        l1_row = l1_data[l1_data['category'] == l1_name]

        if len(l1_row) == 0:
            continue  # Skip if category not in filtered data

        l1_row = l1_row.iloc[0]

        # Abbreviate long names
        display_name = l1_name
        if l1_name == 'OXPHOS subunits':
            display_name = 'OXPHOS subunits'
        elif l1_name == 'OXPHOS assembly factors':
            display_name = 'OXPHOS assembly'

        # Add L1 category
        panel_data.append({
            'label': display_name,
            'type': 'l1_parent',
            'level': 1,
            'total_genes': l1_row['Total_Genes'],
            'degs': l1_row['DEGs'],
            'up': l1_row['Upregulated'],
            'down': l1_row['Downregulated'],
            'pct_up': l1_row['Pct_Up'],
            'pct_down': l1_row['Pct_Down'],
            'dir_sig': l1_row['Directional_Sig'],
            'enr_sig': l1_row['Enrichment_Sig']
        })

        # Add L2 children if any
        for l2_name in l2_children_names:
            l2_row = l2_data[l2_data['category'] == l2_name]

            if len(l2_row) == 0:
                continue  # Skip if not in filtered data

            l2_row = l2_row.iloc[0]

            # Abbreviate L2 names (keep on single line)
            display_name = l2_name
            if 'subunits' in l2_name.lower():
                display_name = l2_name  # Keep as-is, e.g. "CI subunits"
            elif 'assembly factors' in l2_name.lower():
                display_name = l2_name.replace(' assembly factors', ' assembly')  # e.g. "CI assembly"

            # Indent L2 labels
            display_name = '  ' + display_name

            panel_data.append({
                'label': display_name,
                'type': 'l2_child',
                'level': 2,
                'total_genes': l2_row['Total_Genes'],
                'degs': l2_row['DEGs'],
                'up': l2_row['Upregulated'],
                'down': l2_row['Downregulated'],
                'pct_up': l2_row['Pct_Up'],
                'pct_down': l2_row['Pct_Down'],
                'dir_sig': l2_row['Directional_Sig'],
                'enr_sig': l2_row['Enrichment_Sig']
            })

        # Add separator after each L1 branch
        panel_data.append({'type': 'separator'})

    return panel_data


def create_oxphos_hierarchical_figure(stats_df):
    """
    Create OXPHOS hierarchical bar plot figure.

    Parameters
    ----------
    stats_df : pd.DataFrame
        OXPHOS statistics with all levels

    Returns
    -------
    matplotlib.figure.Figure
        Publication-quality figure
    """
    # Setup publication style
    setup_publication_style('paper')

    # Filter to core OXPHOS categories
    filtered_stats = filter_core_oxphos_categories(stats_df)

    print(f"\nFiltered to {len(filtered_stats)} core OXPHOS categories:")
    for level in ['level_0', 'level_1', 'level_2']:
        n_cats = (filtered_stats['level'] == level).sum()
        print(f"  {level}: {n_cats} categories")

    # Build hierarchical data structure
    panel_data = build_hierarchical_data_structure(filtered_stats)

    # Count actual bars (exclude separators)
    n_bars = sum(1 for item in panel_data if item['type'] != 'separator')
    print(f"\nTotal bars to plot: {n_bars}")

    # Create figure
    fig_height = 4.5  # Fixed height per manuscript feedback
    fig = plt.figure(figsize=(8, fig_height))
    ax = fig.add_subplot(111)

    # Visual hierarchy parameters - emphasize bar thickness differences
    visual_params = {
        'parent': {'height': 0.85, 'alpha': 0.9},
        'l1_parent': {'height': 0.6, 'alpha': 0.9},
        'l2_child': {'height': 0.35, 'alpha': 0.9}
    }

    # Single blue color for all bars (downregulated emphasis)
    bar_color = COLORS['down']

    # Plot bars
    y_current = 0
    y_positions = []
    y_labels = []

    for item in panel_data:
        if item['type'] == 'separator':
            y_current -= 0.4  # Gap between branches
            continue

        y_positions.append(y_current)
        y_labels.append(item['label'])

        # Get visual parameters
        params = visual_params.get(item['type'], visual_params['l1_parent'])
        bar_height = params['height']
        alpha = params['alpha']

        # Plot single bar (total DEGs) in uniform blue color
        total_degs = item['degs']
        ax.barh(
            y_current, total_degs,
            height=bar_height,
            color=bar_color,
            alpha=alpha,
            edgecolor='white',
            linewidth=0.5
        )

        # Add DEG count and significance marker together outside bar
        total_degs = item['degs']
        deg_text = f"{item['degs']}/{item['total_genes']}"
        dir_sig = item['dir_sig']

        # Combine deg count and significance marker
        if dir_sig != 'ns':
            combined_text = f"{deg_text}{dir_sig}"
        else:
            combined_text = deg_text

        ax.text(
            total_degs + 1.5, y_current, combined_text,
            ha='left', va='center',
            fontsize=8, color='black'
        )

        y_current -= 1.0  # Move to next bar position

    # Configure axes with bold labels for parent and L1 terms (not L2)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=9)
    # Apply bold to parent and L1 term labels only
    label_idx = 0
    for item in panel_data:
        if item['type'] == 'separator':
            continue
        tick = ax.get_yticklabels()[label_idx]
        if item['type'] in ['parent', 'l1_parent']:
            tick.set_fontweight('bold')
        label_idx += 1
    ax.set_xlabel('Number of DEGs', fontsize=12, fontweight='bold')
    ax.set_title('OXPHOS Pathway Hierarchical Analysis',
                 fontsize=14, fontweight='bold', pad=15)

    # Set x-axis limits
    max_degs = max(item['degs'] for item in panel_data if item['type'] != 'separator')
    ax.set_xlim(0, max_degs + 10)

    # Add grid
    ax.grid(True, axis='x', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    return fig


def main():
    """Main execution function."""

    # Define paths
    base_dir = Path(__file__).parent.parent
    stats_file = base_dir / "results" / "tables" / "oxphos_all_stats.csv"
    output_base = base_dir / "results" / "figures" / "figure_oxphos_hierarchy"

    print("="*70)
    print("OXPHOS Hierarchical Figure Generation")
    print("="*70)
    print(f"\nInput file: {stats_file}")

    # Load statistics
    print("\nLoading OXPHOS statistics...")
    stats_df = pd.read_csv(stats_file)
    print(f"Loaded {len(stats_df)} OXPHOS categories")

    # Create figure
    print("\nGenerating hierarchical figure...")
    fig = create_oxphos_hierarchical_figure(stats_df)

    # Save figure in multiple formats
    output_base.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nSaving figure to {output_base}...")
    saved_files = save_publication_figure(fig, output_base, formats=['pdf', 'png', 'svg'])

    plt.close(fig)

    print("\n" + "="*70)
    print("OXPHOS Figure Generation Complete!")
    print("="*70)
    print(f"\nFigure saved in 3 formats:")
    for f in saved_files:
        print(f"  • {f}")


if __name__ == "__main__":
    main()
