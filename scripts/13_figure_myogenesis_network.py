#!/usr/bin/env python3
"""
Script 23: Myogenesis Literature Gene Network Figure

Visualizes 56 literature-curated myogenesis genes as a TF-target regulatory
network using TFLink data. Nodes encode expression changes (color, size, border)
and are grouped by functional category.

Inputs:
- results/tables/Supplementary_Table_Literature_56_Myogenesis_Genes.csv
- data/tflink_mice.tsv/TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0.tsv

Outputs:
- results/figures/figure_myogenesis_network.{pdf,png,svg}

Author: Analysis pipeline
Date: 2026-03-30
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import networkx as nx
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Setup paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"

import sys
sys.path.append(str(BASE_DIR / "scripts"))
from plotting_style import (
    setup_publication_style, COLORS, FONTS, LAYOUT, EXPORT,
    save_publication_figure, add_panel_label, format_gene_name
)

# ============================================================================
# CATEGORY COLORS (8 categories, colorblind-friendly)
# ============================================================================

CATEGORY_COLORS = {
    'Myogenic Regulatory Factors': '#CC6677',   # Rose
    'Satellite Cell Markers':      '#882255',   # Wine
    'Growth Factors & Signaling':  '#AA3377',   # Magenta
    'Cell Cycle Regulators':       '#EE6677',   # Salmon
    'Myosin Heavy Chains':         '#44AA99',   # Teal
    'Myosin Light Chains':         '#88CCEE',   # Light blue
    'Actin & Troponin Complex':    '#4477AA',   # Blue
    'Sarcomere & Structural':      '#999933',   # Olive
}

CATEGORY_ORDER = [
    'Myogenic Regulatory Factors',
    'Satellite Cell Markers',
    'Growth Factors & Signaling',
    'Cell Cycle Regulators',
    'Myosin Heavy Chains',
    'Myosin Light Chains',
    'Actin & Troponin Complex',
    'Sarcomere & Structural',
]

# Short labels for the figure
CATEGORY_SHORT = {
    'Myogenic Regulatory Factors': 'Myogenic\nRegulatory\nFactors',
    'Satellite Cell Markers':      'Satellite\nCell Markers',
    'Growth Factors & Signaling':  'Growth Factors\n& Signaling',
    'Cell Cycle Regulators':       'Cell Cycle\nRegulators',
    'Myosin Heavy Chains':         'Myosin\nHeavy Chains',
    'Myosin Light Chains':         'Myosin\nLight Chains',
    'Actin & Troponin Complex':    'Actin &\nTroponin',
    'Sarcomere & Structural':      'Sarcomere &\nStructural',
}


def load_gene_data():
    """Load the 56 literature-curated myogenesis genes with expression data.

    The hand-curated gene list (Gene + Category) lives in gene_sets/. Expression
    values (log2FC, padj, DEG status) are pulled fresh from the raw mouse DEG
    file so nothing depends on an intermediate pre-baked supplementary table.
    """
    ref = pd.read_csv(BASE_DIR / "gene_sets" / "literature_56_myogenesis_genes.csv")

    raw = pd.read_csv(DATA_DIR / "mouse_degs.txt", sep="\t")
    # Third unnamed column is the gene symbol (e.g. Gnai3)
    if 'Unnamed: 2' in raw.columns:
        raw['GeneSymbol'] = raw['Unnamed: 2']
    else:
        raw['GeneSymbol'] = raw.iloc[:, 2]

    # European decimal format (comma → period) on numeric columns
    for col in ['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']:
        if col in raw.columns:
            raw[col] = pd.to_numeric(
                raw[col].astype(str).str.replace(',', '.', regex=False),
                errors='coerce',
            )

    merged = ref.merge(
        raw[['GeneSymbol', 'baseMean', 'log2FoldChange', 'padj']],
        how='left', left_on='Gene', right_on='GeneSymbol',
    )

    merged.rename(columns={'Gene': 'Gene Symbol', 'log2FoldChange': 'log2FC'}, inplace=True)
    merged['padj_num'] = merged['padj']

    is_deg = (merged['padj'] < 0.01) & (merged['log2FC'].abs() >= 0.5)
    merged['DEG (padj<0.01, |log2FC|>=0.5)'] = np.where(is_deg, 'Yes', 'No')
    merged['Direction'] = np.where(
        ~is_deg, 'ns',
        np.where(merged['log2FC'] > 0, 'Up', 'Down'),
    )

    missing = merged[merged['GeneSymbol'].isna()]['Gene Symbol'].tolist()
    if missing:
        print(f"WARNING: {len(missing)} literature genes not found in DEG table: {missing}")

    return merged


def load_tflink_edges(gene_set):
    """Load TFLink edges internal to the 56-gene set."""
    tflink_path = DATA_DIR / "tflink_mice.tsv" / "TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0.tsv"
    tf = pd.read_csv(tflink_path, sep='\t')

    # Filter to internal edges (both TF and target in our set)
    internal = tf[
        (tf['Name.TF'].isin(gene_set)) &
        (tf['Name.Target'].isin(gene_set)) &
        (tf['Name.TF'] != tf['Name.Target'])  # remove self-loops
    ].copy()

    # Deduplicate TF-Target pairs, keeping track of evidence quality
    edges = internal.groupby(['Name.TF', 'Name.Target']).agg(
        has_small_scale=('Small-scale.evidence', lambda x: 'Yes' in x.values),
        n_evidence=('Detection.method', 'count')
    ).reset_index()

    return edges


def build_network(gene_df, edges_df):
    """Build networkx DiGraph with node and edge attributes."""
    G = nx.DiGraph()

    # Add nodes
    for _, row in gene_df.iterrows():
        gene = row['Gene Symbol']
        deg_status = row['DEG (padj<0.01, |log2FC|>=0.5)']
        direction = row['Direction']
        log2fc = row['log2FC'] if pd.notna(row['log2FC']) else 0
        category = row['Category']

        # Node color by regulation
        if deg_status == 'Yes' and direction == 'Down':
            node_color = COLORS['down']
        elif deg_status == 'Yes' and direction == 'Up':
            node_color = COLORS['up']
        elif deg_status == 'Not detected':
            node_color = '#FFFFFF'
        else:
            node_color = COLORS['neutral']

        # Node size by |log2FC|
        if deg_status == 'Not detected':
            node_size = 180
        else:
            node_size = 180 + abs(log2fc) * 280

        # Border style
        if deg_status == 'Yes':
            edge_width = 2.5
            edge_color = '#00CC00'  # bright green for DEG significance
        elif deg_status == 'Not detected':
            edge_width = 1.0
            edge_color = '#999999'
        else:
            edge_width = 0.8
            edge_color = '#999999'

        G.add_node(gene,
                   color=node_color,
                   size=node_size,
                   edge_width=edge_width,
                   edge_color=edge_color,
                   category=category,
                   is_deg=(deg_status == 'Yes'),
                   is_detected=(deg_status != 'Not detected'),
                   log2fc=log2fc,
                   direction=direction)

    # Add edges
    for _, row in edges_df.iterrows():
        tf_name = row['Name.TF']
        target = row['Name.Target']
        if tf_name in G and target in G:
            # Check if both endpoints are DEGs
            both_deg = G.nodes[tf_name].get('is_deg', False) and G.nodes[target].get('is_deg', False)
            G.add_edge(tf_name, target,
                       high_confidence=row['has_small_scale'],
                       both_deg=both_deg)

    return G


def compute_grouped_layout(G):
    """
    Compute node positions with category-based grouping.
    Categories arranged in a ring; genes within each category in a local cluster.
    MRFs placed in the center since they are the hub TFs.
    """
    pos = {}

    # Category angles around a ring (MRFs at center, others around)
    # Place categories in a biologically meaningful arrangement:
    # Top: upstream regulators, Bottom: structural effectors
    ring_categories = [
        'Satellite Cell Markers',       # top-left
        'Growth Factors & Signaling',   # left
        'Cell Cycle Regulators',        # bottom-left
        'Myosin Light Chains',          # bottom
        'Actin & Troponin Complex',     # bottom-right
        'Sarcomere & Structural',       # right
        'Myosin Heavy Chains',          # top-right
    ]

    # Get nodes by category
    cat_nodes = {}
    for node in G.nodes():
        cat = G.nodes[node]['category']
        if cat not in cat_nodes:
            cat_nodes[cat] = []
        cat_nodes[cat].append(node)

    # Place MRFs in center cluster - wider spacing to avoid label overlap
    mrf_genes = cat_nodes.get('Myogenic Regulatory Factors', [])
    n_mrf = len(mrf_genes)
    center_radius = 0.9
    for i, gene in enumerate(mrf_genes):
        angle = 2 * np.pi * i / n_mrf - np.pi / 2
        pos[gene] = (center_radius * np.cos(angle), center_radius * np.sin(angle))

    # Place other categories in a ring
    ring_radius = 4.2
    n_ring = len(ring_categories)

    for idx, cat in enumerate(ring_categories):
        genes = cat_nodes.get(cat, [])
        if not genes:
            continue

        # Category center angle
        angle = 2 * np.pi * idx / n_ring - np.pi / 2

        cx = ring_radius * np.cos(angle)
        cy = ring_radius * np.sin(angle)

        # Arrange genes in local cluster
        n = len(genes)
        if n == 1:
            pos[genes[0]] = (cx, cy)
        elif n <= 6:
            # Small sub-ring
            local_radius = 0.35 + 0.1 * n
            for j, gene in enumerate(genes):
                local_angle = 2 * np.pi * j / n
                lx = cx + local_radius * np.cos(local_angle)
                ly = cy + local_radius * np.sin(local_angle)
                pos[gene] = (lx, ly)
        else:
            # Larger categories (like Actin & Troponin with 14):
            # use a grid-like arrangement
            cols = int(np.ceil(np.sqrt(n)))
            rows = int(np.ceil(n / cols))
            spacing = 0.7
            for j, gene in enumerate(genes):
                row = j // cols
                col = j % cols
                lx = cx + (col - (cols - 1) / 2) * spacing
                ly = cy + (row - (rows - 1) / 2) * spacing
                pos[gene] = (lx, ly)

    return pos


def draw_network(ax, G, pos):
    """Render the network on the given axes."""

    # Draw edges first (behind nodes)
    # Separate edges by type for layered drawing
    regular_edges = [(u, v) for u, v, d in G.edges(data=True) if not d.get('both_deg', False)]
    deg_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('both_deg', False)]
    hc_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('high_confidence', False) and not d.get('both_deg', False)]

    # Regular edges: very light
    nx.draw_networkx_edges(
        G, pos, edgelist=regular_edges, ax=ax,
        edge_color='#CCCCCC', alpha=0.2, width=0.4,
        arrows=True, arrowsize=5, arrowstyle='-|>',
        connectionstyle='arc3,rad=0.08',
        min_source_margin=12, min_target_margin=12,
    )

    # High-confidence edges: slightly more visible
    nx.draw_networkx_edges(
        G, pos, edgelist=hc_edges, ax=ax,
        edge_color='#AAAAAA', alpha=0.35, width=0.8,
        arrows=True, arrowsize=6, arrowstyle='-|>',
        connectionstyle='arc3,rad=0.08',
        min_source_margin=12, min_target_margin=12,
    )

    # DEG-to-DEG edges: highlighted
    nx.draw_networkx_edges(
        G, pos, edgelist=deg_edges, ax=ax,
        edge_color='#555555', alpha=0.55, width=1.5,
        arrows=True, arrowsize=8, arrowstyle='-|>',
        connectionstyle='arc3,rad=0.1',
        min_source_margin=12, min_target_margin=12,
    )

    # Draw nodes
    node_colors = [G.nodes[n]['color'] for n in G.nodes()]
    node_sizes = [G.nodes[n]['size'] for n in G.nodes()]
    edge_widths = [G.nodes[n]['edge_width'] for n in G.nodes()]
    edge_colors = [G.nodes[n]['edge_color'] for n in G.nodes()]

    # Draw not-detected nodes with dashed border (separate pass)
    detected_nodes = [n for n in G.nodes() if G.nodes[n]['is_detected']]
    undetected_nodes = [n for n in G.nodes() if not G.nodes[n]['is_detected']]

    # Detected nodes
    nx.draw_networkx_nodes(
        G, pos, nodelist=detected_nodes, ax=ax,
        node_color=[G.nodes[n]['color'] for n in detected_nodes],
        node_size=[G.nodes[n]['size'] for n in detected_nodes],
        edgecolors=[G.nodes[n]['edge_color'] for n in detected_nodes],
        linewidths=[G.nodes[n]['edge_width'] for n in detected_nodes],
    )

    # Undetected nodes (white fill, gray border)
    if undetected_nodes:
        nx.draw_networkx_nodes(
            G, pos, nodelist=undetected_nodes, ax=ax,
            node_color='white',
            node_size=[G.nodes[n]['size'] for n in undetected_nodes],
            edgecolors='#999999',
            linewidths=1.0,
        )

    # Draw all gene labels inside nodes in black
    for node in G.nodes():
        x, y = pos[node]
        is_deg = G.nodes[node]['is_deg']

        if is_deg:
            ax.text(x, y, node,
                    ha='center', va='center',
                    fontsize=FONTS['annotation'],
                    fontweight='bold', fontstyle='italic',
                    color='black',
                    zorder=5)
        else:
            ax.text(x, y, node,
                    ha='center', va='center',
                    fontsize=FONTS['small'] - 0.5,
                    fontstyle='italic',
                    color='black',
                    zorder=5)


def draw_category_patches(ax, G, pos):
    """Draw light background patches for each category cluster."""
    from matplotlib.patches import FancyBboxPatch
    from scipy.spatial import ConvexHull

    for cat in CATEGORY_ORDER:
        nodes = [n for n in G.nodes() if G.nodes[n]['category'] == cat]
        if not nodes:
            continue

        coords = np.array([pos[n] for n in nodes])
        color = CATEGORY_COLORS[cat]

        # Compute bounding area
        cx, cy = coords.mean(axis=0)
        if len(coords) >= 3:
            # Use convex hull + padding
            try:
                hull = ConvexHull(coords)
                hull_pts = coords[hull.vertices]
                # Expand hull outward
                expanded = []
                for pt in hull_pts:
                    direction = pt - np.array([cx, cy])
                    norm = np.linalg.norm(direction)
                    if norm > 0:
                        expanded.append(pt + direction / norm * 0.45)
                    else:
                        expanded.append(pt + 0.45)
                expanded = np.array(expanded)
                from matplotlib.patches import Polygon
                patch = Polygon(expanded, closed=True,
                                facecolor=color, alpha=0.08,
                                edgecolor=color, linewidth=1.0, linestyle='--',
                                zorder=0)
                ax.add_patch(patch)
            except Exception:
                # Fallback: circle around center
                radius = np.max(np.linalg.norm(coords - [cx, cy], axis=1)) + 0.5
                circle = plt.Circle((cx, cy), radius, facecolor=color, alpha=0.08,
                                    edgecolor=color, linewidth=1.0, linestyle='--', zorder=0)
                ax.add_patch(circle)
        elif len(coords) == 2:
            radius = np.linalg.norm(coords[0] - coords[1]) / 2 + 0.5
            circle = plt.Circle((cx, cy), radius, facecolor=color, alpha=0.08,
                                edgecolor=color, linewidth=1.0, linestyle='--', zorder=0)
            ax.add_patch(circle)
        else:
            circle = plt.Circle((cx, cy), 0.5, facecolor=color, alpha=0.08,
                                edgecolor=color, linewidth=1.0, linestyle='--', zorder=0)
            ax.add_patch(circle)

        # Category label - always placed below the gene group
        label_x = cx
        label_y = coords[:, 1].min() - 0.6

        ax.text(label_x, label_y, CATEGORY_SHORT[cat],
                ha='center', va='top',
                fontsize=FONTS['small'] + 0.5,
                fontweight='bold',
                color=color, alpha=0.9,
                zorder=1)


def draw_legend(ax):
    """Draw a clean legend panel with round circles matching the network."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_aspect('equal')
    ax.axis('off')

    y = 9.5
    dy = 0.55
    x_circle = 1.0
    r_node = 0.22  # circle radius matching network nodes
    x_text = 1.8

    # --- Node Color ---
    ax.text(0.2, y, 'Node Color (Regulation)', fontsize=FONTS['annotation'] + 1,
            fontweight='bold', va='center')
    y -= dy * 1.2

    legend_items = [
        (COLORS['down'], 'black', 2.0, 'Downregulated DEG'),
        (COLORS['up'], 'black', 2.0, 'Upregulated DEG'),
        (COLORS['neutral'], '#999999', 0.8, 'Not significant'),
        ('white', '#999999', 0.8, 'Not detected'),
    ]
    for color, ec, lw, label in legend_items:
        circle = plt.Circle((x_circle, y), r_node, facecolor=color,
                            edgecolor=ec, linewidth=lw, transform=ax.transData)
        ax.add_patch(circle)
        ax.text(x_text, y, label, fontsize=FONTS['annotation'],
                va='center', transform=ax.transData)
        y -= dy

    # --- Node Size ---
    y -= dy * 0.4
    ax.text(0.2, y, 'Node Size', fontsize=FONTS['annotation'] + 1,
            fontweight='bold', va='center')
    y -= dy * 1.2

    for r, label in [(0.15, 'Small |log2FC|'), (0.32, 'Large |log2FC|')]:
        circle = plt.Circle((x_circle, y), r, facecolor=COLORS['neutral'],
                            edgecolor='#999999', linewidth=0.8, transform=ax.transData)
        ax.add_patch(circle)
        ax.text(x_text, y, label, fontsize=FONTS['annotation'],
                va='center', transform=ax.transData)
        y -= dy

    # --- Node Border ---
    y -= dy * 0.4
    ax.text(0.2, y, 'Node Border', fontsize=FONTS['annotation'] + 1,
            fontweight='bold', va='center')
    y -= dy * 1.2

    circle = plt.Circle((x_circle, y), r_node, facecolor=COLORS['neutral'],
                        edgecolor='#00CC00', linewidth=2.5, transform=ax.transData)
    ax.add_patch(circle)
    ax.text(x_text, y, 'DEG (green border)', fontsize=FONTS['annotation'],
            va='center', transform=ax.transData)
    y -= dy

    circle = plt.Circle((x_circle, y), r_node, facecolor=COLORS['neutral'],
                        edgecolor='#999999', linewidth=0.8, transform=ax.transData)
    ax.add_patch(circle)
    ax.text(x_text, y, 'Non-DEG (thin gray)', fontsize=FONTS['annotation'],
            va='center', transform=ax.transData)

    # --- Edges ---
    y -= dy * 1.4
    ax.text(0.2, y, 'Edges (TFLink)', fontsize=FONTS['annotation'] + 1,
            fontweight='bold', va='center')
    y -= dy * 1.2

    edge_items = [
        ('#555555', 1.5, 0.55, 'DEG \u2192 DEG'),
        ('#AAAAAA', 0.8, 0.35, 'High-confidence'),
        ('#CCCCCC', 0.4, 0.2, 'ChIP-seq evidence'),
    ]
    for color, lw, alpha, label in edge_items:
        ax.annotate('', xy=(1.4, y), xytext=(0.4, y),
                    arrowprops=dict(arrowstyle='-|>', color=color, lw=lw,
                                    alpha=alpha),
                    transform=ax.transData)
        ax.text(x_text, y, label, fontsize=FONTS['annotation'],
                va='center', transform=ax.transData)
        y -= dy


def main():
    setup_publication_style('paper')

    print("\n" + "=" * 70)
    print("GENERATING MYOGENESIS NETWORK FIGURE")
    print("=" * 70)

    # Load data
    print("\nLoading gene expression data...")
    gene_df = load_gene_data()
    gene_set = set(gene_df['Gene Symbol'].tolist())
    print(f"  Loaded {len(gene_df)} genes")

    print("Loading TFLink edges...")
    edges_df = load_tflink_edges(gene_set)
    print(f"  Found {len(edges_df)} internal edges")
    print(f"  High-confidence: {edges_df['has_small_scale'].sum()}")

    # Build network
    print("Building network...")
    G = build_network(gene_df, edges_df)
    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    n_deg = sum(1 for n in G.nodes() if G.nodes[n]['is_deg'])
    n_down = sum(1 for n in G.nodes() if G.nodes[n]['direction'] == 'Down')
    n_up = sum(1 for n in G.nodes() if G.nodes[n]['direction'] == 'Up')
    print(f"  DEGs: {n_deg} ({n_down} down, {n_up} up)")

    # Compute layout
    print("Computing layout...")
    pos = compute_grouped_layout(G)

    # Create figure
    print("Drawing figure...")
    fig = plt.figure(figsize=(LAYOUT['full_page'], 9.5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2.8, 1], wspace=0.05)

    # Panel A: Network
    ax_net = fig.add_subplot(gs[0])
    ax_net.set_aspect('equal')
    ax_net.axis('off')

    draw_category_patches(ax_net, G, pos)
    draw_network(ax_net, G, pos)

    # Auto-scale axes with symmetric ranges
    all_x = [p[0] for p in pos.values()]
    all_y = [p[1] for p in pos.values()]
    margin = 1.2
    x_center = (max(all_x) + min(all_x)) / 2
    y_center = (max(all_y) + min(all_y)) / 2
    half_range = max(max(all_x) - min(all_x), max(all_y) - min(all_y)) / 2 + margin
    ax_net.set_xlim(x_center - half_range, x_center + half_range)
    ax_net.set_ylim(y_center - half_range, y_center + half_range)

    # Panel B: Legend
    ax_leg = fig.add_subplot(gs[1])
    draw_legend(ax_leg)

    plt.tight_layout()

    # Title centered over the network panel's actual position
    net_bbox = ax_net.get_position()
    title_x = (net_bbox.x0 + net_bbox.x1) / 2
    fig.text(title_x, net_bbox.y1 + 0.04,
             'TF-Target Regulatory Network of 56 Literature-Curated Myogenesis Genes',
             fontsize=FONTS['title'], fontweight='bold',
             ha='center', va='bottom')

    # Save
    print("\nSaving figure...")
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    save_publication_figure(fig, FIGURES_DIR / "figure_myogenesis_network",
                            formats=['pdf', 'png', 'svg'])

    plt.close()

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == '__main__':
    main()
