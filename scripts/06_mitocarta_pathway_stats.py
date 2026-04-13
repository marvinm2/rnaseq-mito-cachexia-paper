#!/usr/bin/env python3
"""
MitoCarta Pathway Statistical Analysis
=======================================

This script calculates enrichment and directional bias statistics for all
MitoCarta pathways across 3 hierarchical levels using Fisher's exact test.

Two independent statistical tests per pathway:
1. Enrichment test (one-tailed): Tests if DEGs are overrepresented in the
   pathway compared to genome-wide DEG rate (25.6%)
2. Directional bias test (two-tailed): Tests if up/down regulation ratio
   differs from global distribution (55.5% up / 44.5% down)

FDR correction applied using Benjamini-Hochberg method separately for each
level and test type.

Input:
------
- Pathway gene lists from mitocarta_pathway_gene_extraction.py
- DEG data from mitochondrial_degs.csv

Output:
-------
- pathway_statistics_level_0.csv: Statistics for 7 Level 0 pathways
- pathway_statistics_level_1.csv: Statistics for 20 Level 1 pathways
- pathway_statistics_level_2.csv: Statistics for 78 Level 2 pathways
- pathway_statistics_all.csv: Combined statistics for all 105 pathways

Author: Bioinformatics Team
Date: December 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')


# Global constants (from RNA-seq analysis)
TOTAL_GENES = 11599
TOTAL_DEGS = 2974
TOTAL_UP_GLOBAL = 1650
TOTAL_DOWN_GLOBAL = 1324


def load_pathway_genes(gene_list_file):
    """
    Load gene list for a pathway from CSV file.

    Parameters
    ----------
    gene_list_file : Path
        Path to pathway gene list CSV

    Returns
    -------
    list
        List of gene symbols in the pathway
    """
    df = pd.read_csv(gene_list_file)
    return df['gene_symbol'].tolist()


def calculate_pathway_statistics(pathway_genes, all_degs, pathway_name, level):
    """
    Calculate enrichment and directional statistics for a pathway.

    Parameters
    ----------
    pathway_genes : list
        List of gene symbols in the pathway
    all_degs : pd.DataFrame
        DataFrame with columns: GeneSymbol, log2FoldChange, regulation
    pathway_name : str
        Name of the pathway
    level : str
        Hierarchy level (level_0, level_1, level_2)

    Returns
    -------
    dict
        Dictionary with statistics including:
        - pathway, level, Total_Genes, DEGs, Upregulated, Downregulated
        - Pct_Up, Pct_Down, Enrichment_P, Directional_P, odds ratios
    """
    # Match genes to DEGs
    pathway_degs = all_degs[all_degs['GeneSymbol'].isin(pathway_genes)].copy()

    M = len(pathway_genes)  # Total genes in pathway
    k = len(pathway_degs)    # DEGs in pathway
    n_up = (pathway_degs['regulation'] == 'up').sum()
    n_down = (pathway_degs['regulation'] == 'down').sum()

    # Enrichment test (one-tailed)
    # H0: DEG rate in pathway = global DEG rate
    # H1: DEG rate in pathway > global DEG rate
    table_enrich = [
        [k, TOTAL_DEGS - k],                              # DEGs: in pathway, not in pathway
        [M - k, TOTAL_GENES - M - TOTAL_DEGS + k]        # non-DEGs: in pathway, not in pathway
    ]
    try:
        or_enrich, p_enrich = fisher_exact(table_enrich, alternative='greater')
    except:
        or_enrich, p_enrich = np.nan, 1.0

    # Directional bias test (two-tailed)
    # H0: up/down ratio in pathway = global ratio
    # H1: up/down ratio in pathway ≠ global ratio
    if k > 0:
        table_dir = [
            [n_up, TOTAL_UP_GLOBAL - n_up],              # Up: in pathway, not in pathway
            [n_down, TOTAL_DOWN_GLOBAL - n_down]         # Down: in pathway, not in pathway
        ]
        try:
            or_dir, p_dir = fisher_exact(table_dir, alternative='two-sided')
        except:
            or_dir, p_dir = np.nan, 1.0
    else:
        or_dir, p_dir = np.nan, 1.0

    return {
        'pathway': pathway_name,
        'level': level,
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


def apply_fdr_correction_by_level(stats_df):
    """
    Apply FDR correction separately for each level and test type.

    Rationale:
    - Level 0: 7 pathways
    - Level 1: 20 pathways
    - Level 2: 78 pathways

    Separate correction accounts for different pathway counts at each level.

    Parameters
    ----------
    stats_df : pd.DataFrame
        DataFrame with statistics for all pathways

    Returns
    -------
    pd.DataFrame
        DataFrame with added FDR columns and significance markers
    """
    # Initialize FDR columns
    stats_df['Enrichment_FDR'] = np.nan
    stats_df['Directional_FDR'] = np.nan

    for level in ['level_0', 'level_1', 'level_2']:
        level_mask = stats_df['level'] == level
        level_data = stats_df[level_mask].copy()

        if len(level_data) == 0:
            continue

        # Enrichment FDR
        valid_enrich = level_data['Enrichment_P'].notna()
        if valid_enrich.any():
            reject, pvals_fdr, _, _ = multipletests(
                level_data.loc[valid_enrich, 'Enrichment_P'],
                alpha=0.05,
                method='fdr_bh'
            )
            stats_df.loc[level_mask & stats_df.index.isin(level_data.index[valid_enrich]), 'Enrichment_FDR'] = pvals_fdr

        # Directional FDR
        valid_dir = level_data['Directional_P'].notna()
        if valid_dir.any():
            reject, pvals_fdr, _, _ = multipletests(
                level_data.loc[valid_dir, 'Directional_P'],
                alpha=0.05,
                method='fdr_bh'
            )
            stats_df.loc[level_mask & stats_df.index.isin(level_data.index[valid_dir]), 'Directional_FDR'] = pvals_fdr

    # Add significance markers
    stats_df['Enrichment_Sig'] = stats_df['Enrichment_FDR'].apply(get_sig_marker)
    stats_df['Directional_Sig'] = stats_df['Directional_FDR'].apply(get_sig_marker)

    return stats_df


def get_sig_marker(fdr):
    """
    Convert FDR p-value to significance marker.

    Parameters
    ----------
    fdr : float
        FDR-corrected p-value

    Returns
    -------
    str
        Significance marker: ***, **, *, or ns
    """
    if pd.isna(fdr):
        return 'ns'
    elif fdr < 0.001:
        return '***'
    elif fdr < 0.01:
        return '**'
    elif fdr < 0.05:
        return '*'
    else:
        return 'ns'


def main():
    """Main execution function."""

    # Define paths
    base_dir = Path(__file__).parent.parent
    gene_dir = base_dir / "results" / "tables" / "pathway_genes"
    deg_file = base_dir / "results" / "tables" / "mitochondrial_degs.csv"
    output_dir = base_dir / "results" / "tables"

    print("="*70)
    print("MitoCarta Pathway Statistical Analysis")
    print("="*70)
    print(f"\nInput directory: {gene_dir}")
    print(f"DEG file: {deg_file}")
    print(f"Output directory: {output_dir}")

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

    # Process each level
    all_stats = []

    for level_num in range(3):
        level_name = f'level_{level_num}'

        print(f"\n{'='*70}")
        print(f"Processing {level_name.upper()}")
        print(f"{'='*70}")

        # Find all gene list files for this level
        pattern = f"{level_name}_*.csv"
        gene_files = list(gene_dir.glob(pattern))

        print(f"Found {len(gene_files)} pathways")

        # Calculate statistics for each pathway
        for gene_file in sorted(gene_files):
            # Extract pathway name from filename
            # Format: level_X_Pathway_Name_genes.csv
            filename = gene_file.stem  # Remove .csv
            pathway_name = filename.replace(f'{level_name}_', '').replace('_genes', '').replace('_', ' ')

            # Load genes
            pathway_genes = load_pathway_genes(gene_file)

            # Calculate statistics
            stats = calculate_pathway_statistics(
                pathway_genes, degs, pathway_name, level_name
            )
            all_stats.append(stats)

            # Print progress
            pct_down = stats['Pct_Down']
            print(f"  {pathway_name}: {stats['DEGs']}/{stats['Total_Genes']} DEGs ({pct_down:.1f}% down)")

    # Create DataFrame
    stats_df = pd.DataFrame(all_stats)

    # Apply FDR correction
    print(f"\n{'='*70}")
    print("Applying FDR correction (Benjamini-Hochberg)")
    print(f"{'='*70}")
    stats_df = apply_fdr_correction_by_level(stats_df)

    # Reorder columns for clarity
    column_order = [
        'level', 'pathway',
        'Total_Genes', 'DEGs', 'Upregulated', 'Downregulated',
        'Pct_Up', 'Pct_Down',
        'Enrichment_P', 'Enrichment_FDR', 'Enrichment_OR', 'Enrichment_Sig',
        'Directional_P', 'Directional_FDR', 'Directional_OR', 'Directional_Sig'
    ]
    stats_df = stats_df[column_order]

    # Save combined statistics
    output_file_all = output_dir / "pathway_statistics_all.csv"
    stats_df.to_csv(output_file_all, index=False)
    print(f"\nSaved combined statistics: {output_file_all}")

    # Save level-specific statistics
    for level_num in range(3):
        level_name = f'level_{level_num}'
        level_data = stats_df[stats_df['level'] == level_name]

        output_file = output_dir / f"pathway_statistics_{level_name}.csv"
        level_data.to_csv(output_file, index=False)
        print(f"Saved {level_name} statistics: {output_file}")

    # Print summary
    print(f"\n{'='*70}")
    print("Statistical Analysis Summary")
    print(f"{'='*70}")

    for level in ['level_0', 'level_1', 'level_2']:
        level_data = stats_df[stats_df['level'] == level]
        print(f"\n{level.upper().replace('_', ' ')}:")
        print(f"  Total pathways: {len(level_data)}")
        print(f"  Enrichment significant (FDR<0.05): {(level_data['Enrichment_Sig'] != 'ns').sum()}")
        print(f"  Directional significant (FDR<0.05): {(level_data['Directional_Sig'] != 'ns').sum()}")

        # Show top pathways by DEG count
        print(f"\n  Top pathways by DEG count:")
        for idx, row in level_data.nlargest(5, 'DEGs').iterrows():
            print(f"    {row['pathway']}: {row['DEGs']} DEGs ({row['Pct_Down']:.1f}% down) "
                  f"[Enrich: {row['Enrichment_Sig']}, Dir: {row['Directional_Sig']}]")

    print(f"\n{'='*70}")
    print("Pathway statistical analysis complete!")
    print(f"{'='*70}")
    print(f"\nOutput files:")
    print(f"  - pathway_statistics_all.csv ({len(stats_df)} pathways)")
    print(f"  - pathway_statistics_level_0.csv")
    print(f"  - pathway_statistics_level_1.csv")
    print(f"  - pathway_statistics_level_2.csv")


if __name__ == "__main__":
    main()
