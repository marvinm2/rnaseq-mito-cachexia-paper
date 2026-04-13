#!/usr/bin/env python3
"""
MitoCarta OXPHOS Statistical Analysis
======================================

This script calculates enrichment and directional bias statistics for all
OXPHOS categories using Fisher's exact test.

Two independent statistical tests:
1. Enrichment test (one-tailed): Tests if DEGs are overrepresented in the
   category compared to genome-wide DEG rate (25.6%)
2. Directional bias test (two-tailed): Tests if up/down regulation ratio
   differs from global distribution (55.5% up / 44.5% down)

FDR correction applied using Benjamini-Hochberg method.

Input:
------
- OXPHOS gene lists from mitocarta_oxphos_annotation.py
- DEG data from mitochondrial_degs.csv

Output:
-------
- oxphos_all_stats.csv: Combined statistics for all OXPHOS categories

Author: Bioinformatics Team
Date: December 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


# Global constants (from RNA-seq analysis)
TOTAL_GENES = 11599
TOTAL_DEGS = 2974
TOTAL_UP_GLOBAL = 1650
TOTAL_DOWN_GLOBAL = 1324


def calculate_oxphos_category_stats(category_genes, all_degs, category_name, level_name):
    """
    Calculate enrichment and directional statistics for an OXPHOS category.

    Parameters
    ----------
    category_genes : list
        List of gene symbols in the category
    all_degs : pd.DataFrame
        DataFrame with columns: GeneSymbol, log2FoldChange, regulation
    category_name : str
        Name of the OXPHOS category
    level_name : str
        Hierarchy level (level_0, level_1, level_2)

    Returns
    -------
    dict
        Dictionary with statistics including:
        - Total_Genes, DEGs, Upregulated, Downregulated
        - Pct_Up, Pct_Down
        - Enrichment_P, Directional_P
    """
    # Match genes to DEGs
    category_degs = all_degs[all_degs['GeneSymbol'].isin(category_genes)].copy()

    M = len(category_genes)  # Total genes in category
    k = len(category_degs)    # DEGs in category
    n_up = (category_degs['regulation'] == 'up').sum()
    n_down = (category_degs['regulation'] == 'down').sum()

    # Enrichment test (one-tailed)
    # H0: DEG rate in category = global DEG rate
    # H1: DEG rate in category > global DEG rate
    table_enrich = [
        [k, TOTAL_DEGS - k],                              # DEGs: in category, not in category
        [M - k, TOTAL_GENES - M - TOTAL_DEGS + k]        # non-DEGs: in category, not in category
    ]
    try:
        or_enrich, p_enrich = fisher_exact(table_enrich, alternative='greater')
    except:
        or_enrich, p_enrich = np.nan, 1.0

    # Directional bias test (two-tailed)
    # H0: up/down ratio in category = global ratio
    # H1: up/down ratio in category ≠ global ratio
    if k > 0:
        table_dir = [
            [n_up, TOTAL_UP_GLOBAL - n_up],              # Up: in category, not in category
            [n_down, TOTAL_DOWN_GLOBAL - n_down]         # Down: in category, not in category
        ]
        try:
            or_dir, p_dir = fisher_exact(table_dir, alternative='two-sided')
        except:
            or_dir, p_dir = np.nan, 1.0
    else:
        or_dir, p_dir = np.nan, 1.0

    return {
        'category': category_name,
        'level': level_name,
        'Total_Genes': M,
        'DEGs': k,
        'Upregulated': n_up,
        'Downregulated': n_down,
        'Pct_Up': (n_up / k * 100) if k > 0 else 0,
        'Pct_Down': (n_down / k * 100) if k > 0 else 0,
        'Enrichment_P': p_enrich,
        'Directional_P': p_dir,
        'Enrichment_OR': or_enrich,
        'Directional_OR': or_dir
    }


def add_significance_markers(df, p_col, sig_col):
    """
    Add significance markers based on FDR-corrected p-values.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with p-value column
    p_col : str
        Name of FDR-corrected p-value column
    sig_col : str
        Name of new column to store significance markers

    Returns
    -------
    pd.DataFrame
        DataFrame with added significance column
    """
    def get_marker(p):
        if pd.isna(p):
            return 'ns'
        elif p < 0.001:
            return '***'
        elif p < 0.01:
            return '**'
        elif p < 0.05:
            return '*'
        else:
            return 'ns'

    df[sig_col] = df[p_col].apply(get_marker)
    return df


def load_oxphos_gene_lists(gene_dir):
    """
    Load all OXPHOS gene lists from CSV files.

    Parameters
    ----------
    gene_dir : Path
        Directory containing OXPHOS gene list CSV files

    Returns
    -------
    dict
        Dictionary with structure:
        {
            'level_0': {category_name: [genes]},
            'level_1': {category_name: [genes]},
            'level_2': {category_name: [genes]}
        }
    """
    gene_lists = {'level_0': {}, 'level_1': {}, 'level_2': {}}

    for level in ['level_0', 'level_1', 'level_2']:
        # Find all files for this level
        pattern = f"oxphos_{level}_*.csv"
        files = list(gene_dir.glob(pattern))

        for file in files:
            df = pd.read_csv(file)
            category = df['category'].iloc[0]
            genes = df['gene_symbol'].tolist()
            gene_lists[level][category] = genes

    return gene_lists


def main():
    """Main execution function."""

    # Define paths
    base_dir = Path(__file__).parent.parent
    gene_dir = base_dir / "results" / "tables" / "oxphos_genes"
    deg_file = base_dir / "results" / "tables" / "mitochondrial_degs.csv"
    output_file = base_dir / "results" / "tables" / "oxphos_all_stats.csv"

    print("="*70)
    print("MitoCarta OXPHOS Statistical Analysis")
    print("="*70)
    print(f"\nInput directory: {gene_dir}")
    print(f"DEG file: {deg_file}")
    print(f"Output file: {output_file}")

    # Load DEG data
    print("\nLoading DEG data...")
    degs = pd.read_csv(deg_file)

    # Add regulation column if not present
    if 'regulation' not in degs.columns:
        degs['regulation'] = degs['log2FoldChange'].apply(
            lambda x: 'up' if x > 0 else 'down'
        )

    print(f"Loaded {len(degs)} mitochondrial DEGs")
    print(f"  Up-regulated: {(degs['regulation'] == 'up').sum()}")
    print(f"  Down-regulated: {(degs['regulation'] == 'down').sum()}")

    # Load OXPHOS gene lists
    print("\nLoading OXPHOS gene lists...")
    gene_lists = load_oxphos_gene_lists(gene_dir)

    total_categories = sum(len(cats) for cats in gene_lists.values())
    print(f"Loaded {total_categories} OXPHOS categories:")
    for level, cats in gene_lists.items():
        print(f"  {level}: {len(cats)} categories")

    # Calculate statistics for all categories
    print("\nCalculating statistics...")
    all_stats = []

    for level_name in ['level_0', 'level_1', 'level_2']:
        categories = gene_lists[level_name]

        for category_name, genes in categories.items():
            stats = calculate_oxphos_category_stats(
                genes, degs, category_name, level_name
            )
            all_stats.append(stats)

            # Print progress
            pct_down = stats['Pct_Down']
            print(f"  {level_name} - {category_name}: {stats['DEGs']}/{stats['Total_Genes']} DEGs ({pct_down:.1f}% down)")

    # Create DataFrame
    stats_df = pd.DataFrame(all_stats)

    # Apply FDR correction separately for each test type
    print("\nApplying FDR correction (Benjamini-Hochberg)...")

    # Enrichment FDR
    valid_enrich = stats_df['Enrichment_P'].notna()
    if valid_enrich.any():
        reject, pvals_fdr, _, _ = multipletests(
            stats_df.loc[valid_enrich, 'Enrichment_P'],
            alpha=0.05,
            method='fdr_bh'
        )
        stats_df.loc[valid_enrich, 'Enrichment_FDR'] = pvals_fdr
    else:
        stats_df['Enrichment_FDR'] = np.nan

    # Directional FDR
    valid_dir = stats_df['Directional_P'].notna()
    if valid_dir.any():
        reject, pvals_fdr, _, _ = multipletests(
            stats_df.loc[valid_dir, 'Directional_P'],
            alpha=0.05,
            method='fdr_bh'
        )
        stats_df.loc[valid_dir, 'Directional_FDR'] = pvals_fdr
    else:
        stats_df['Directional_FDR'] = np.nan

    # Add significance markers
    stats_df = add_significance_markers(stats_df, 'Enrichment_FDR', 'Enrichment_Sig')
    stats_df = add_significance_markers(stats_df, 'Directional_FDR', 'Directional_Sig')

    # Reorder columns for clarity
    column_order = [
        'level', 'category',
        'Total_Genes', 'DEGs', 'Upregulated', 'Downregulated',
        'Pct_Up', 'Pct_Down',
        'Enrichment_P', 'Enrichment_FDR', 'Enrichment_OR', 'Enrichment_Sig',
        'Directional_P', 'Directional_FDR', 'Directional_OR', 'Directional_Sig'
    ]
    stats_df = stats_df[column_order]

    # Sort by level and DEG count
    stats_df = stats_df.sort_values(['level', 'DEGs'], ascending=[True, False])

    # Save results
    print(f"\nSaving results to {output_file}...")
    stats_df.to_csv(output_file, index=False)

    # Print summary
    print("\n" + "="*70)
    print("Statistical Analysis Summary")
    print("="*70)

    for level in ['level_0', 'level_1', 'level_2']:
        level_data = stats_df[stats_df['level'] == level]
        print(f"\n{level.upper().replace('_', ' ')}:")
        print(f"  Categories: {len(level_data)}")
        print(f"  Enrichment significant (FDR<0.05): {(level_data['Enrichment_Sig'] != 'ns').sum()}")
        print(f"  Directional significant (FDR<0.05): {(level_data['Directional_Sig'] != 'ns').sum()}")

        # Show top categories by DEG count
        print(f"\n  Top categories by DEG count:")
        for idx, row in level_data.head(5).iterrows():
            print(f"    {row['category']}: {row['DEGs']} DEGs ({row['Pct_Down']:.1f}% down) "
                  f"[Enrich: {row['Enrichment_Sig']}, Dir: {row['Directional_Sig']}]")

    print("\n" + "="*70)
    print("OXPHOS statistical analysis complete!")
    print("="*70)
    print(f"\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
