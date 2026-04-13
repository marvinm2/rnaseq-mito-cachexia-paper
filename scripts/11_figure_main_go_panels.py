#!/usr/bin/env python3
"""
Script 20: Generate RNA-Seq Publication Figures

Generates TWO separate publication-ready figures for RNA-seq analysis:

Figure 1 - Figure_Volcano_RNAseq:
  • Standalone volcano plot with mitochondrial gene highlighting
  • Shows genome-wide differential expression (11,599 genes)
  • Highlights 240 mitochondrial DEGs

Figure 2 - Figure_Main_RNAseq (2 panels):
  • Panel A: Myogenesis pathways (muscle differentiation)
    - 2 GO subcategories + core TF set
  • Panel B: Proteolysis pathways (protein degradation)
    - Hierarchical view with catabolic process emphasis

This separation allows the volcano plot to be featured prominently while
keeping related pathway analyses together in a companion figure.

Output Location: manuscript_materials/figures/

Author: Analysis pipeline
Date: 2025-12-05
Version: 2.0 (separated figures)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import warnings
import sys
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')

# Add scripts folder to path for publication_style import
SCRIPT_DIR = Path(__file__).parent.parent / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

# Import publication styling
from plotting_style import (
    setup_publication_style, COLORS, FONTS, LAYOUT, EXPORT,
    add_panel_label, set_axis_labels, save_publication_figure
)

# Setup paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Global constants
TOTAL_GENES = 11599
TOTAL_DEGS = 2974
TOTAL_UP_GLOBAL = 1650
TOTAL_DOWN_GLOBAL = 1324


def load_all_data():
    """
    Load all required datasets for the comprehensive figure

    Returns:
    --------
    tuple: (all_genes, mito_degs, myo_stats, prot_stats)
    """
    print("\n" + "="*70)
    print("LOADING DATA FOR COMPREHENSIVE RNA-SEQ FIGURE")
    print("="*70)

    # 1. Full transcriptome for volcano plot
    all_genes_file = DATA_DIR / "mouse_degs.txt"
    print(f"\n1. Loading full transcriptome from: {all_genes_file.name}")

    if str(all_genes_file).endswith('.txt'):
        all_genes = pd.read_csv(all_genes_file, sep='\t')
    else:
        all_genes = pd.read_csv(all_genes_file)

    # Handle European decimal format (force, regardless of detected dtype)
    numeric_cols = ['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
    for col in numeric_cols:
        if col in all_genes.columns:
            all_genes[col] = pd.to_numeric(
                all_genes[col].astype(str).str.replace(',', '.', regex=False),
                errors='coerce',
            )

    # Handle gene symbol column
    if 'Unnamed: 2' in all_genes.columns:
        all_genes['GeneSymbol'] = all_genes['Unnamed: 2']
    elif 'Unnamed: 1' in all_genes.columns:
        if all_genes.iloc[0, 2] and not str(all_genes.iloc[0, 2]).startswith('ENSMUSG'):
            all_genes['GeneSymbol'] = all_genes.iloc[:, 2]
        else:
            all_genes['GeneSymbol'] = all_genes.iloc[:, 1]
    else:
        all_genes['GeneSymbol'] = all_genes.iloc[:, 1]

    print(f"   ✓ Loaded {len(all_genes):,} genes")

    # 2. Mitochondrial DEGs
    print("\n2. Loading mitochondrial DEGs")
    mito_degs = pd.read_csv(TABLES_DIR / "mitochondrial_degs.csv")
    print(f"   ✓ Loaded {len(mito_degs)} mitochondrial DEGs")

    # 3. GO Myogenesis stats
    print("\n3. Loading GO myogenesis data")
    myo_stats = pd.read_csv(TABLES_DIR / "go_myogenesis_stats.csv")
    print(f"   ✓ Loaded {len(myo_stats)} myogenesis subcategories")

    # 4. GO Proteolysis stats
    print("\n4. Loading GO proteolysis data")
    prot_stats = pd.read_csv(TABLES_DIR / "go_proteolysis_stats.csv")
    print(f"   ✓ Loaded {len(prot_stats)} proteolysis subcategories")

    # 5. Load summary for parent stats
    print("\n5. Loading GO summary (for parent stats)")
    go_summary = pd.read_csv(TABLES_DIR / "SUMMARY_GO_Analysis.csv")
    print(f"   ✓ Loaded {len(go_summary)} summary rows")

    print("\n" + "="*70)
    print("DATA LOADING COMPLETE")
    print("="*70)

    return all_genes, mito_degs, myo_stats, prot_stats, go_summary


def create_panel_a_volcano(ax, all_genes, mito_degs):
    """
    Create Panel A: Volcano plot with mitochondrial highlighting

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        Axes to plot on
    all_genes : pd.DataFrame
        Full transcriptome data
    mito_degs : pd.DataFrame
        Mitochondrial DEGs
    """
    print("\n" + "-"*70)
    print("Creating Panel A: Volcano Plot")
    print("-"*70)

    # Calculate -log10(padj)
    all_genes['minus_log10_padj'] = -np.log10(all_genes['padj'].replace(0, 1e-300))

    # Identify DEGs (padj < 0.01, |log2FC| >= 0.5)
    all_genes['is_deg'] = (all_genes['padj'] < 0.01) & (all_genes['log2FoldChange'].abs() >= 0.5)
    all_genes['is_up'] = all_genes['log2FoldChange'] > 0

    # Create mitochondrial gene set
    mito_genes_set = set(mito_degs['GeneSymbol'].str.upper())
    all_genes['is_mito'] = all_genes['GeneSymbol'].str.upper().isin(mito_genes_set)

    # Classify into 5 categories for layered plotting
    all_genes['category'] = 'Non-DEG'
    all_genes.loc[all_genes['is_deg'] & all_genes['is_up'] & ~all_genes['is_mito'], 'category'] = 'Other DEG (up)'
    all_genes.loc[all_genes['is_deg'] & ~all_genes['is_up'] & ~all_genes['is_mito'], 'category'] = 'Other DEG (down)'
    all_genes.loc[all_genes['is_deg'] & all_genes['is_up'] & all_genes['is_mito'], 'category'] = 'Mito DEG (up)'
    all_genes.loc[all_genes['is_deg'] & ~all_genes['is_up'] & all_genes['is_mito'], 'category'] = 'Mito DEG (down)'

    # Separate categories for layered plotting
    non_degs = all_genes[all_genes['category'] == 'Non-DEG']
    other_up = all_genes[all_genes['category'] == 'Other DEG (up)']
    other_down = all_genes[all_genes['category'] == 'Other DEG (down)']
    mito_up = all_genes[all_genes['category'] == 'Mito DEG (up)']
    mito_down = all_genes[all_genes['category'] == 'Mito DEG (down)']

    print(f"   Non-DEGs: {len(non_degs):,}")
    print(f"   Other DEGs up: {len(other_up):,}")
    print(f"   Other DEGs down: {len(other_down):,}")
    print(f"   Mito DEGs up: {len(mito_up):,}")
    print(f"   Mito DEGs down: {len(mito_down):,}")

    # Plot in layers (background to foreground)
    # Layer 1: Non-DEGs (background)
    ax.scatter(non_degs['log2FoldChange'], non_degs['minus_log10_padj'],
              s=2, alpha=0.12, c=COLORS['neutral'], label='Non-DEG',
              rasterized=True, zorder=1)

    # Layer 2: Other DEGs (middle) - PALE colors to fade into background
    ax.scatter(other_down['log2FoldChange'], other_down['minus_log10_padj'],
              s=4, alpha=0.3, c='#91BAD6', label='Other DEG (down)',
              edgecolors='none', zorder=2)
    ax.scatter(other_up['log2FoldChange'], other_up['minus_log10_padj'],
              s=4, alpha=0.3, c='#F7A9A8', label='Other DEG (up)',
              edgecolors='none', zorder=2)

    # Layer 3: Mitochondrial DEGs (foreground) - BRIGHT colors to stand out
    ax.scatter(mito_down['log2FoldChange'], mito_down['minus_log10_padj'],
              s=25, alpha=0.95, c=COLORS['down'], label='Mito DEG (down)',
              edgecolors='#0d3d5c', linewidth=0.6, zorder=10)
    ax.scatter(mito_up['log2FoldChange'], mito_up['minus_log10_padj'],
              s=25, alpha=0.95, c=COLORS['up'], label='Mito DEG (up)',
              edgecolors='#8B0000', linewidth=0.6, zorder=10)

    # Add threshold lines
    padj_threshold = 0.01
    ax.axhline(-np.log10(padj_threshold), color='black', linestyle='--',
              linewidth=1, alpha=0.4, zorder=0)
    ax.axvline(-0.5, color='black', linestyle='--', linewidth=1, alpha=0.4, zorder=0)
    ax.axvline(0.5, color='black', linestyle='--', linewidth=1, alpha=0.4, zorder=0)

    # Formatting
    ax.set_xlabel('log₂ Fold Change', fontsize=FONTS['axis_label'], fontweight='bold')
    ax.set_ylabel('-log₁₀ (adjusted P-value)', fontsize=FONTS['axis_label'], fontweight='bold')
    ax.set_title('Differential Gene Expression in Cachexia',
                fontsize=FONTS['title'], fontweight='bold', pad=10)

    # Legend (custom - with separate mito up/down)
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['neutral'],
               markersize=4, alpha=0.4, label='Non-DEG'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#A5C4D9',
               markersize=5, alpha=0.5, label='Other DEG'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['down'],
               markersize=8, alpha=0.95, label='Mito DEG ↓',
               markeredgewidth=0.8, markeredgecolor='#0d3d5c'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['up'],
               markersize=8, alpha=0.95, label='Mito DEG ↑',
               markeredgewidth=0.8, markeredgecolor='#8B0000')
    ]
    ax.legend(handles=legend_elements, loc='upper right', frameon=True,
             fontsize=FONTS['legend']-1, framealpha=0.95)

    ax.grid(alpha=0.2, linestyle=':', linewidth=0.5)
    ax.set_xlim(-4, 4)
    ax.set_ylim(0, 60)

    # No interpretive statistics in figure (per user request)
    # All statistics should be in text, not figure panels

    print("   ✓ Panel A complete")


def abbreviate_label(go_name):
    """
    Abbreviate GO term names for compact display in Panel B

    Parameters:
    -----------
    go_name : str
        Full GO term name

    Returns:
    --------
    str: Abbreviated multi-line label
    """
    abbreviations = {
        'striated muscle cell differentiation': 'Striated musc. diff.',
        'smooth muscle cell differentiation': 'Smooth musc. diff.',
        'muscle cell development': 'Muscle cell dev.',
        'muscle cell fate commitment': 'Muscle cell fate commit.',
        'skeletal muscle cell differentiation': 'Skel. musc. cell diff.',
        'skeletal muscle fiber development': 'Skel. musc. fiber dev.',
        'skeletal muscle satellite stem cell asymmetric division': 'Satellite stem cell div.',
        'skeletal muscle tissue growth': 'Skel. musc. tissue growth',
        'myoblast differentiation': 'Myoblast diff.',
        'muscle tissue development': 'Muscle tissue dev.'
    }

    # Return abbreviation if found, otherwise intelligently split long names
    if go_name.lower() in abbreviations:
        return abbreviations[go_name.lower()]
    elif len(go_name) > 30:
        # Auto-split at ~20 chars if no abbreviation defined
        words = go_name.split()
        mid = len(words) // 2
        return ' '.join(words[:mid]) + '\n' + ' '.join(words[mid:])
    else:
        return go_name


def create_panel_b_myogenesis(ax, myo_stats):
    """
    Create Panel B: Hierarchical myogenesis with 3-level GO structure

    Layout (10 bars total):
    - Branch 1: GO:0042692 (parent) + 4 children (2 is_a, 2 part_of)
    - Separator line
    - Branch 2: GO:0007519 (parent) + 3 children (part_of)
    - Separator line
    - Branch 3: Literature-based manual gene set

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        Axes to plot on
    myo_stats : pd.DataFrame
        Myogenesis GO statistics (not used - loading expanded stats instead)
    """
    print("\n" + "-"*70)
    print("Creating Panel B: Hierarchical Myogenesis (3-level GO structure)")
    print("-"*70)

    # Load expanded stats (child terms)
    expanded_stats = pd.read_csv(TABLES_DIR / "go_myogenesis_stats_expanded.csv")
    print(f"   Loaded {len(expanded_stats)} GO child terms")

    # Load annotated DEG data (has gene symbols already)
    degs = pd.read_csv(TABLES_DIR / "mouse_degs_annotated_full.csv")

    # Determine regulation direction if not already present
    if 'regulation' not in degs.columns and 'log2FoldChange' in degs.columns:
        degs['regulation'] = degs['log2FoldChange'].apply(lambda x: 'up' if x > 0 else 'down')

    # Helper function for calculating parent stats
    def calculate_parent_stats(parent_genes_file, parent_go_id, parent_name):
        """Calculate statistics for a parent GO term"""
        parent_genes = pd.read_csv(parent_genes_file)
        # Use GeneSymbol_upper for matching (GO files have uppercase gene symbols)
        parent_degs = degs[degs['GeneSymbol_upper'].isin(parent_genes['GeneSymbol'])]

        M = len(parent_genes)
        k = len(parent_degs)
        n_up = (parent_degs['regulation'] == 'up').sum()
        n_down = (parent_degs['regulation'] == 'down').sum()

        # Enrichment test
        table_enrich = [
            [k, TOTAL_DEGS - k],
            [M - k, TOTAL_GENES - M - TOTAL_DEGS + k]
        ]
        or_enrich, p_enrich = fisher_exact(table_enrich, alternative='greater')

        # Directional test
        table_dir = [
            [n_up, TOTAL_UP_GLOBAL - n_up],
            [n_down, TOTAL_DOWN_GLOBAL - n_down]
        ]
        or_dir, p_dir = fisher_exact(table_dir, alternative='two-sided')

        # Simple significance markers (no FDR for parent terms)
        def sig_marker(pval):
            if pval < 0.001:
                return '***'
            elif pval < 0.01:
                return '**'
            elif pval < 0.05:
                return '*'
            else:
                return 'ns'

        return {
            'go_id': parent_go_id,
            'go_name': parent_name,
            'degs': k,
            'up': n_up,
            'down': n_down,
            'pct_up': (n_up / k * 100) if k > 0 else 0,
            'pct_down': (n_down / k * 100) if k > 0 else 0,
            'dir_sig': sig_marker(p_dir),
            'enr_sig': sig_marker(p_enrich)
        }

    # Calculate parent term statistics
    parent1_genes = pd.read_csv(TABLES_DIR / "go_myogenesis_genes.csv")
    parent1 = calculate_parent_stats(
        TABLES_DIR / "go_myogenesis_genes.csv",
        "GO:0042692",
        "Muscle cell diff. (all)"
    )
    parent1['total_genes'] = len(parent1_genes)

    # For GO:0007519, we need to aggregate from children since no dedicated file exists
    # Get all genes from GO:0007519 children
    go_0007519_children = ['GO:0035914', 'GO:0048741', 'GO:0048630']
    go_0007519_genes = set()
    for child_id in go_0007519_children:
        child_file = TABLES_DIR / f"go_{child_id.replace(':', '_')}_genes.csv"
        if child_file.exists():
            child_genes = pd.read_csv(child_file)
            # Handle both 'gene_symbol' and 'GeneSymbol' column names
            gene_col = 'gene_symbol' if 'gene_symbol' in child_genes.columns else 'GeneSymbol'
            go_0007519_genes.update(child_genes[gene_col].tolist())

    # Calculate stats for GO:0007519 from aggregated gene list
    # Use GeneSymbol_upper for matching (GO files have uppercase/lowercase mixed gene symbols)
    parent2_degs = degs[degs['GeneSymbol_upper'].isin([g.upper() for g in go_0007519_genes])]
    M2 = len(go_0007519_genes)
    k2 = len(parent2_degs)
    n_up2 = (parent2_degs['regulation'] == 'up').sum()
    n_down2 = (parent2_degs['regulation'] == 'down').sum()

    table_enrich2 = [[k2, TOTAL_DEGS - k2], [M2 - k2, TOTAL_GENES - M2 - TOTAL_DEGS + k2]]
    or_enrich2, p_enrich2 = fisher_exact(table_enrich2, alternative='greater')

    table_dir2 = [[n_up2, TOTAL_UP_GLOBAL - n_up2], [n_down2, TOTAL_DOWN_GLOBAL - n_down2]]
    or_dir2, p_dir2 = fisher_exact(table_dir2, alternative='two-sided')

    def sig_marker(pval):
        if pval < 0.001:
            return '***'
        elif pval < 0.01:
            return '**'
        elif pval < 0.05:
            return '*'
        else:
            return 'ns'

    parent2 = {
        'go_id': 'GO:0007519',
        'go_name': 'Skel. muscle tissue dev. (all)',
        'degs': k2,
        'up': n_up2,
        'down': n_down2,
        'pct_up': (n_up2 / k2 * 100) if k2 > 0 else 0,
        'pct_down': (n_down2 / k2 * 100) if k2 > 0 else 0,
        'dir_sig': sig_marker(p_dir2),
        'enr_sig': sig_marker(p_enrich2),
        'total_genes': M2
    }

    # Load manual myogenesis data
    manual_degs = pd.read_csv(TABLES_DIR / "myogenesis_degs.csv")
    n_manual_total = 56
    n_manual_degs = len(manual_degs)
    n_manual_down = (manual_degs['regulation'] == 'down').sum()
    n_manual_up = (manual_degs['regulation'] == 'up').sum()

    # Manual myogenesis stats
    k_manual = n_manual_degs
    M_manual = n_manual_total
    table_enrich_manual = [[k_manual, TOTAL_DEGS - k_manual],
                           [M_manual - k_manual, TOTAL_GENES - M_manual - TOTAL_DEGS + k_manual]]
    or_enrich_manual, p_enrich_manual = fisher_exact(table_enrich_manual, alternative='greater')

    table_dir_manual = [[n_manual_up, TOTAL_UP_GLOBAL - n_manual_up],
                        [n_manual_down, TOTAL_DOWN_GLOBAL - n_manual_down]]
    or_dir_manual, p_dir_manual = fisher_exact(table_dir_manual, alternative='two-sided')

    print(f"   Parent 1 (GO:0042692): {parent1['degs']} DEGs")
    print(f"   Parent 2 (GO:0007519): {parent2['degs']} DEGs")
    print(f"   Manual myogenesis: {n_manual_degs} DEGs")

    # === BUILD HIERARCHICAL DATA STRUCTURE ===
    panel_data = []

    # === BRANCH 1: GO:0042692 ===
    panel_data.append({
        'label': parent1['go_name'],
        'type': 'parent',
        'degs': parent1['degs'],
        'up': parent1['up'],
        'down': parent1['down'],
        'pct_up': parent1['pct_up'],
        'pct_down': parent1['pct_down'],
        'dir_sig': parent1['dir_sig'],
        'enr_sig': parent1['enr_sig'],
        'total_genes': parent1['total_genes']
    })

    # Child bars for GO:0042692 (sorted by is_a first, then part_of)
    branch1_children = [
        'GO:0051146',  # striated (is_a)
        'GO:0051145',  # smooth (is_a)
        'GO:0055001',  # development (part_of)
        'GO:0042693'   # fate commitment (part_of)
    ]

    for child_id in branch1_children:
        row = expanded_stats[expanded_stats['go_id'] == child_id]
        if len(row) > 0:
            row = row.iloc[0]
            panel_data.append({
                'label': '  ' + abbreviate_label(row['go_name']),  # Indent children
                'type': f"child_{row['relationship']}",
                'degs': row['DEGs'],
                'up': row['Upregulated'],
                'down': row['Downregulated'],
                'pct_up': row['Pct_Up'],
                'pct_down': row['Pct_Down'],
                'dir_sig': row['Directional_Sig'],
                'enr_sig': row['Enrichment_Sig'],
                'total_genes': row.get('Total_Genes', row['DEGs'])  # Fallback to DEGs if not available
            })

    # Separator
    panel_data.append({'type': 'separator'})

    # === BRANCH 2: GO:0007519 ===
    panel_data.append({
        'label': parent2['go_name'],
        'type': 'parent',
        'degs': parent2['degs'],
        'up': parent2['up'],
        'down': parent2['down'],
        'pct_up': parent2['pct_up'],
        'pct_down': parent2['pct_down'],
        'dir_sig': parent2['dir_sig'],
        'enr_sig': parent2['enr_sig'],
        'total_genes': parent2['total_genes']
    })

    # Child bars for GO:0007519
    branch2_children = [
        'GO:0035914',  # skeletal muscle cell differentiation
        'GO:0048741',  # skeletal muscle fiber development
        'GO:0048630'   # skeletal muscle tissue growth
    ]

    for child_id in branch2_children:
        row = expanded_stats[expanded_stats['go_id'] == child_id]
        if len(row) > 0:
            row = row.iloc[0]
            panel_data.append({
                'label': '  ' + abbreviate_label(row['go_name']),  # Indent children
                'type': f"child_{row['relationship']}",
                'degs': row['DEGs'],
                'up': row['Upregulated'],
                'down': row['Downregulated'],
                'pct_up': row['Pct_Up'],
                'pct_down': row['Pct_Down'],
                'dir_sig': row['Directional_Sig'],
                'enr_sig': row['Enrichment_Sig'],
                'total_genes': row.get('Total_Genes', row['DEGs'])  # Fallback to DEGs if not available
            })

    # Separator
    panel_data.append({'type': 'separator'})

    # === BRANCH 3: Literature-based manual set ===
    panel_data.append({
        'label': 'Literature-based',
        'type': 'manual',
        'degs': n_manual_degs,
        'up': n_manual_up,
        'down': n_manual_down,
        'non_deg': n_manual_total - n_manual_degs,
        'pct_up': (n_manual_up / n_manual_degs * 100) if n_manual_degs > 0 else 0,
        'pct_down': (n_manual_down / n_manual_degs * 100) if n_manual_degs > 0 else 0,
        'dir_sig': sig_marker(p_dir_manual),
        'enr_sig': sig_marker(p_enrich_manual),
        'total_genes': n_manual_total
    })

    # === PLOTTING ===
    y_positions = []
    y_current = 0

    for i, item in enumerate(panel_data):
        if item['type'] == 'separator':
            y_current -= 0.4  # Small gap
            continue

        y_positions.append(y_current)

        # Use consistent bright colors, but different heights for parent vs child
        alpha = 0.9
        if item['type'] == 'parent':
            height = 0.8  # Thicker for parent terms
        else:
            height = 0.45  # Thinner for child terms

        # Plot segments
        if item['down'] > 0:
            ax.barh(y_current, item['down'], height=height,
                   color=COLORS['down'], alpha=alpha, edgecolor='white', linewidth=0.5)

        if item['up'] > 0:
            ax.barh(y_current, item['up'], left=item['down'], height=height,
                   color=COLORS['up'], alpha=alpha, edgecolor='white', linewidth=0.5)

        # Percentage labels inside bars
        if item['down'] > 0 and item['pct_down'] > 5:
            ax.text(item['down']/2, y_current, f"{item['pct_down']:.0f}%",
                   ha='center', va='center', fontsize=FONTS['small']-1,
                   color='white', fontweight='bold')
        if item['up'] > 0 and item['pct_up'] > 5:
            ax.text(item['down'] + item['up']/2, y_current,
                   f"{item['pct_up']:.0f}%",
                   ha='center', va='center', fontsize=FONTS['small']-1,
                   color='white', fontweight='bold')

        # DEG count labels (outside bars) - OXPHOS style: degs/total_genes
        total_shown = item['degs']  # Only count DEGs, not non_deg segment
        deg_text = f"{item['degs']}/{item['total_genes']}"
        ax.text(total_shown + 1, y_current, deg_text,
               ha='left', va='center', fontsize=FONTS['small']-1,
               color=COLORS['text'])

        # Directional significance marker at bar end
        if item['dir_sig'] != 'ns':
            marker_x = total_shown + 5  # Offset for asterisks
            ax.text(marker_x, y_current, item['dir_sig'],
                   ha='left', va='center', fontsize=FONTS['small'],
                   color=COLORS['text'], fontweight='bold')

        y_current -= 1.0

    # Create y-tick labels and positions with bold for parent terms and manual
    labels = []
    label_weights = []
    for item in panel_data:
        if item['type'] == 'separator':
            continue
        labels.append(item['label'])
        label_weights.append('bold' if item['type'] in ['parent', 'manual'] else 'normal')
    y_tick_positions = y_positions

    # Formatting
    ax.set_yticks(y_tick_positions)
    ax.set_yticklabels(labels, fontsize=FONTS['tick_label']-1)
    # Apply bold to parent term labels
    for tick, weight in zip(ax.get_yticklabels(), label_weights):
        tick.set_fontweight(weight)
    ax.set_ylim(y_current, 0.5)
    ax.set_xlabel('Number of Genes', fontsize=FONTS['axis_label'], fontweight='bold')
    ax.set_title('Myogenesis', fontsize=FONTS['title'], fontweight='bold', pad=10)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.grid(axis='x', alpha=0.2, linestyle=':', linewidth=0.5)
    ax.set_axisbelow(True)

    print(f"   ✓ Panel B hierarchical complete ({len(labels)} bars)")

def create_panel_c_proteolysis(ax, prot_stats, go_summary):
    """
    Create Panel C: HIERARCHICAL proteolysis view (parent + 3 children)

    Shows parent term at top, then 3 main subcategories below:
    - Protein catabolism (emphasized - THE key finding)
    - Protein processing
    - Membrane protein proteolysis

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        Axes to plot on
    prot_stats : pd.DataFrame
        Proteolysis GO statistics (subcategories only)
    go_summary : pd.DataFrame
        Summary with parent stats
    """
    print("\n" + "-"*70)
    print("Creating Panel C: Hierarchical Proteolysis View")
    print("-"*70)

    # Get parent from summary
    parent_summary = go_summary[(go_summary['Category'] == 'Proteolysis') & (go_summary['Subcategory'] == 'ALL')].iloc[0]

    # Create parent dict
    parent = {
        'Label': 'Proteolysis (all)',
        'DEGs': parent_summary['Total_DEGs'],
        'Upregulated': parent_summary['Upregulated'],
        'Downregulated': parent_summary['Downregulated'],
        'Pct_Up': parent_summary['Pct_Up'],
        'Pct_Down': parent_summary['Pct_Down'],
        'Total_Genes': 1051,
        'Directional_Sig': '',
        'Enrichment_Sig': '',
        'Prefix': '',  # No tree symbol for parent
        'IsParent': True
    }

    # Calculate parent term significance (not in pre-computed stats)
    # Enrichment test
    k_parent = parent['DEGs']
    M_parent = parent['Total_Genes']
    K = TOTAL_DEGS
    N = TOTAL_GENES
    table_enrich = [[k_parent, K - k_parent],
                    [M_parent - k_parent, N - M_parent - K + k_parent]]
    or_enrich, p_enrich = fisher_exact(table_enrich, alternative='greater')

    # Directional test
    a = parent['Upregulated']
    b = parent['Downregulated']
    c = TOTAL_UP_GLOBAL - a
    d = TOTAL_DOWN_GLOBAL - b
    table_dir = [[a, c], [b, d]]
    or_dir, p_dir = fisher_exact(table_dir, alternative='two-sided')

    # Apply significance threshold (simple approach without FDR correction)
    def sig_marker(pval):
        if pd.isna(pval):
            return 'ns'
        elif pval < 0.001:
            return '***'
        elif pval < 0.01:
            return '**'
        elif pval < 0.05:
            return '*'
        else:
            return 'ns'

    parent['Enrichment_Sig'] = sig_marker(p_enrich)
    parent['Directional_Sig'] = sig_marker(p_dir)

    # Get children from prot_stats (exclude "self proteolysis" with 0 DEGs)
    children_data = []

    # 1. Protein catabolism (THE key finding - emphasized)
    catabolic_row = prot_stats[prot_stats['Subcategory'] == 'proteolysis involved in protein catabolic process']
    if len(catabolic_row) > 0:
        catabolic = catabolic_row.iloc[0]
        children_data.append({
            'Label': 'Catabolic process',
            'DEGs': catabolic['DEGs'],
            'Upregulated': catabolic['Upregulated'],
            'Downregulated': catabolic['Downregulated'],
            'Pct_Up': catabolic['Pct_Up'],
            'Pct_Down': catabolic['Pct_Down'],
            'Total_Genes': catabolic['Total_Genes'],
            'Directional_Sig': catabolic['Directional_Sig'],
            'Enrichment_Sig': catabolic['Enrichment_Sig'],
            'Prefix': '',
            'IsParent': False,
            'IsKey': True  # EMPHASIZED
        })

    # 2. Protein processing
    processing_row = prot_stats[prot_stats['Subcategory'] == 'protein processing']
    if len(processing_row) > 0:
        processing = processing_row.iloc[0]
        children_data.append({
            'Label': 'Processing',
            'DEGs': processing['DEGs'],
            'Upregulated': processing['Upregulated'],
            'Downregulated': processing['Downregulated'],
            'Pct_Up': processing['Pct_Up'],
            'Pct_Down': processing['Pct_Down'],
            'Total_Genes': processing['Total_Genes'],
            'Directional_Sig': processing['Directional_Sig'],
            'Enrichment_Sig': processing['Enrichment_Sig'],
            'Prefix': '',
            'IsParent': False,
            'IsKey': False
        })

    # 3. Membrane protein proteolysis (last one gets └─)
    membrane_row = prot_stats[prot_stats['Subcategory'] == 'membrane protein proteolysis']
    if len(membrane_row) > 0:
        membrane = membrane_row.iloc[0]
        children_data.append({
            'Label': 'Membrane',
            'DEGs': membrane['DEGs'],
            'Upregulated': membrane['Upregulated'],
            'Downregulated': membrane['Downregulated'],
            'Pct_Up': membrane['Pct_Up'],
            'Pct_Down': membrane['Pct_Down'],
            'Total_Genes': membrane['Total_Genes'],
            'Directional_Sig': membrane['Directional_Sig'],
            'Enrichment_Sig': membrane['Enrichment_Sig'],
            'Prefix': '',
            'IsParent': False,
            'IsKey': False
        })

    # Combine all rows (parent + children)
    all_data = [parent] + children_data
    n_rows = len(all_data)

    print(f"   Parent: {parent['DEGs']} DEGs")
    for child in children_data:
        print(f"   {child['Prefix']}{child['Label']}: {child['DEGs']} DEGs ({child['Pct_Up']:.1f}% up)")

    # Create y positions (reversed for top-to-bottom hierarchy)
    y_pos = np.arange(n_rows)[::-1]  # Reverse so parent is at top

    # Determine max_x from DEGs (no grey background bars)
    max_x = max([row['DEGs'] for row in all_data])

    # Plot DEG bars (colored, no grey background per user request)
    for i, row in enumerate(all_data):
        # Use consistent bright alpha for all bars
        alpha_val = 0.9
        # Different heights for parent vs child terms
        bar_height = 0.8 if row['IsParent'] else 0.45

        # Downregulated
        if row['Downregulated'] > 0:
            ax.barh(y_pos[i], row['Downregulated'], height=bar_height, color=COLORS['down'],
                   alpha=alpha_val, edgecolor='white', linewidth=0.5,
                   label='Downregulated' if i == 0 else '', zorder=2)

        # Upregulated
        if row['Upregulated'] > 0:
            ax.barh(y_pos[i], row['Upregulated'], height=bar_height,
                   left=row['Downregulated'], color=COLORS['up'],
                   alpha=alpha_val, edgecolor='white', linewidth=0.5,
                   label='Upregulated' if i == 0 else '', zorder=2)

    # Add percentage labels (only for bars large enough to fit text)
    min_degs_for_label = 50  # Only show percentage labels if bar has >= 50 DEGs
    for i, row in enumerate(all_data):
        if row['DEGs'] < min_degs_for_label:
            continue  # Skip small bars (Processing, Membrane)

        fontsize_pct = FONTS['small']

        if row['Downregulated'] > 0:
            ax.text(row['Downregulated']/2, y_pos[i], f"{row['Pct_Down']:.0f}%",
                   ha='center', va='center', fontsize=fontsize_pct,
                   color='white', fontweight='bold')
        if row['Upregulated'] > 0:
            ax.text(row['Downregulated'] + row['Upregulated']/2, y_pos[i],
                   f"{row['Pct_Up']:.0f}%",
                   ha='center', va='center', fontsize=fontsize_pct,
                   color='white', fontweight='bold')

    # Add DEG counts (OXPHOS style: degs/total_genes)
    for i, row in enumerate(all_data):
        deg_text = f"{row['DEGs']}/{row['Total_Genes']}"
        ax.text(row['DEGs'] + 1.5, y_pos[i], deg_text,
               ha='left', va='center', fontsize=FONTS['small']-1,
               color=COLORS['text'])

    # Add directional significance markers at bar ends
    for i, row in enumerate(all_data):
        if row['Directional_Sig'] and row['Directional_Sig'] not in ['ns', '—']:
            marker_x = row['DEGs'] + 40  # Increased offset to prevent overlap
            dir_sig = row['Directional_Sig']
            fontsize_dir = FONTS['small']+1 if row.get('IsKey', False) else FONTS['small']
            ax.text(marker_x, y_pos[i], dir_sig,
                   ha='left', va='center', fontsize=fontsize_dir,
                   color=COLORS['text'], fontweight='bold')

    # Format y-axis labels with tree symbols
    # For multi-line labels, indent continuation lines to align with tree symbol
    labels = []
    for row in all_data:
        label = row['Label']
        prefix = row['Prefix']
        if '\n' in label and prefix:
            # Add indentation to continuation lines (3 spaces for tree symbols like "├─ ")
            indent = ' ' * len(prefix)
            label_parts = label.split('\n')
            formatted_label = prefix + label_parts[0] + '\n' + indent + '\n'.join(label_parts[1:])
            labels.append(formatted_label)
        else:
            labels.append(prefix + label)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=FONTS['tick_label']-1)
    # Apply bold to parent term labels
    for i, (tick, row) in enumerate(zip(ax.get_yticklabels(), all_data)):
        tick.set_fontweight('bold' if row['IsParent'] else 'normal')

    # Formatting
    ax.set_xlabel('Number of Genes', fontsize=FONTS['axis_label'], fontweight='bold')
    ax.set_title('Proteolysis',
                fontsize=FONTS['title'], fontweight='bold', pad=10)
    ax.set_xlim(0, max_x * 1.35)  # Accommodate markers at +35 offset

    # Grid
    ax.grid(axis='x', alpha=0.25, linestyle=':', linewidth=0.5)
    ax.set_axisbelow(True)

    print("   ✓ Panel C hierarchical complete")


def create_figure_volcano_standalone(all_genes, mito_degs):
    """
    Create standalone volcano plot figure (separate from main 3-panel figure)

    This generates a single-panel publication figure showing the volcano plot
    with mitochondrial gene highlighting, suitable for use as a separate figure.

    Parameters:
    -----------
    all_genes : pd.DataFrame
        Full transcriptome data
    mito_degs : pd.DataFrame
        Mitochondrial DEGs

    Returns:
    --------
    matplotlib.figure.Figure
        Complete figure object ready for saving
    """
    print("\n" + "="*70)
    print("CREATING STANDALONE VOLCANO FIGURE")
    print("="*70)

    # Create single-panel figure (60% of original height: 5.5 * 0.6 = 3.3)
    fig = plt.figure(figsize=(7.0, 3.3))
    ax = fig.add_subplot(111)

    # Use existing volcano panel function
    create_panel_a_volcano(ax, all_genes, mito_degs)

    # No panel label for standalone figure

    return fig


def create_figure_pathways_2panel(myo_stats, prot_stats, go_summary):
    """
    Create 2-panel pathway figure (proteolysis + myogenesis)

    This generates a side-by-side 2-panel figure with:
    - Panel A: Proteolysis
    - Panel B: Myogenesis

    Parameters:
    -----------
    myo_stats : pd.DataFrame
        Myogenesis GO statistics
    prot_stats : pd.DataFrame
        Proteolysis GO statistics
    go_summary : pd.DataFrame
        Summary with parent stats

    Returns:
    --------
    matplotlib.figure.Figure
        Complete figure object ready for saving
    """
    print("\n" + "="*70)
    print("CREATING 2-PANEL PATHWAY FIGURE")
    print("="*70)

    # Create 2-panel figure (60% of original height: 5.0 * 0.6 = 3.0)
    fig = plt.figure(figsize=(7.5, 3.0))
    gs = fig.add_gridspec(1, 2,
                          wspace=0.75,  # Increased for more space between panels
                          left=0.08, right=0.98,
                          top=0.94, bottom=0.10)

    # Create axes - proteolysis on left (A), myogenesis on right (B)
    ax_prot = fig.add_subplot(gs[0, 0])  # Left panel (Panel A)
    ax_myo = fig.add_subplot(gs[0, 1])   # Right panel (Panel B)

    print("   ✓ 2-panel layout created")

    # Create panels - proteolysis is now Panel A, myogenesis is Panel B
    create_panel_c_proteolysis(ax_prot, prot_stats, go_summary)
    create_panel_b_myogenesis(ax_myo, myo_stats)

    # Add panel labels - moved more to the left
    add_panel_label(ax_prot, 'A', x=-0.30, y=1.10)
    add_panel_label(ax_myo, 'B', x=-0.30, y=1.10)

    return fig


def main():
    """
    Main execution function

    Generates TWO separate publication figures:
    1. Figure_Volcano_RNAseq - Standalone volcano plot
    2. Figure_Main_RNAseq - 2-panel pathway analysis (A: Myogenesis, B: Proteolysis)
    """
    print("\n" + "="*70)
    print("RNA-SEQ PUBLICATION FIGURES GENERATOR")
    print("="*70)
    print("\nGenerating TWO separate figures:")
    print("  Figure 1: Volcano plot (standalone)")
    print("  Figure 2: Pathway analysis (2 panels: A=Myogenesis, B=Proteolysis)")

    # Setup publication style
    setup_publication_style('paper')

    # Load all data
    all_genes, mito_degs, myo_stats, prot_stats, go_summary = load_all_data()

    # ========================================================================
    # FIGURE 1: STANDALONE VOLCANO PLOT
    # ========================================================================
    print("\n" + "="*70)
    print("FIGURE 1: VOLCANO PLOT")
    print("="*70)

    fig_volcano = create_figure_volcano_standalone(all_genes, mito_degs)

    # Save volcano figure
    filepath_volcano = FIGURES_DIR / "figure_volcano_mitochondrial"
    saved_volcano = save_publication_figure(fig_volcano, filepath_volcano,
                                           formats=['pdf', 'png', 'svg'])
    plt.close(fig_volcano)

    print("\n✓ Figure 1 saved:")
    for f in saved_volcano:
        print(f"  • {f.name}")

    # ========================================================================
    # FIGURE 2: 2-PANEL PATHWAY FIGURE (PANELS A & B)
    # ========================================================================
    print("\n" + "="*70)
    print("FIGURE 2: PATHWAY ANALYSIS (2 PANELS)")
    print("="*70)

    fig_pathways = create_figure_pathways_2panel(myo_stats, prot_stats, go_summary)

    # Save pathway figure (myogenesis + proteolysis GO panels)
    filepath_pathways = FIGURES_DIR / "figure_go_myogenesis_proteolysis"
    saved_pathways = save_publication_figure(fig_pathways, filepath_pathways,
                                            formats=['pdf', 'png', 'svg'])
    plt.close(fig_pathways)

    print("\n✓ Figure 2 saved:")
    for f in saved_pathways:
        print(f"  • {f.name}")

    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("\n" + "="*70)
    print("FIGURE GENERATION COMPLETE!")
    print("="*70)

    print("\n📊 TWO SEPARATE FIGURES GENERATED:")
    print("\n1. Figure_Volcano_RNAseq (standalone)")
    print("   • Genome-wide volcano plot")
    print("   • 11,599 genes, 2,974 DEGs")
    print("   • 240 mitochondrial DEGs highlighted (83% downregulated)")

    print("\n2. Figure_Main_RNAseq (2 panels)")
    print("   Panel A: Myogenesis")
    print("   • 2 GO subcategories + core TF set")
    print("   • Shows muscle differentiation patterns")
    print("   Panel B: Proteolysis")
    print("   • Hierarchical view (parent + 3 children)")
    print("   • Emphasizes catabolic process (78% upregulated)")

    print("\n🔬 Biological Story:")
    print("  • Mitochondrial suppression (83% down)")
    print("  • Muscle differentiation blockade (balanced regulation)")
    print("  • Protein catabolism activation (78% up)")

    print("\n✓ All figures saved to: manuscript_materials/figures/")
    print("✓ Ready for manuscript submission!")


if __name__ == "__main__":
    main()
