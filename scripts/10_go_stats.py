#!/usr/bin/env python3
"""
Script 17: Statistical Analysis of GO-Annotated Gene Sets

Performs comprehensive statistical analysis of proteolysis and myogenesis
GO-annotated gene sets against the DEG dataset.

Statistical Tests:
1. Enrichment Test (Fisher's exact, one-tailed): Tests if DEGs are overrepresented
   in each GO category/subcategory
2. Directional Bias Test (Fisher's exact, two-tailed): Tests if the ratio of
   upregulated vs downregulated DEGs differs from the global distribution

FDR Correction:
- Benjamini-Hochberg method applied across ALL p-values in each analysis
- Separate markers for enrichment (E) and directional (D) significance

Outputs:
- DEG lists matched to GO categories
- Per-subcategory statistical summaries
- Validation report comparing with manual curation

Author: Analysis pipeline
Date: 2025-11-28
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Setup paths
BASE_DIR = Path(__file__).parent.parent
RESULTS_DIR = BASE_DIR / "results" / "tables"

# Global dataset constants
TOTAL_GENES = 11599   # Total expressed genes
TOTAL_DEGS = 2974     # Total DEGs (FDR < 0.01)
TOTAL_UP_GLOBAL = 1650    # Global upregulated (55.5%)
TOTAL_DOWN_GLOBAL = 1324  # Global downregulated (44.5%)


def load_deg_data():
    """Load and prepare DEG dataset"""
    print("\n" + "="*70)
    print("LOADING DEG DATA")
    print("="*70)

    # Load DEGs
    degs = pd.read_csv(RESULTS_DIR / "mouse_degs_filtered.csv")

    # Handle column naming issues
    if 'GeneSymbol' not in degs.columns:
        if 'Unnamed: 2' in degs.columns:
            degs['GeneSymbol'] = degs['Unnamed: 2']
        elif 'Unnamed: 1' in degs.columns:
            degs['GeneSymbol'] = degs['Unnamed: 1']

    # Ensure uppercase gene symbols
    degs['GeneSymbol'] = degs['GeneSymbol'].str.upper()

    # Determine regulation direction
    if 'regulation' not in degs.columns:
        if 'log2FoldChange' in degs.columns:
            degs['regulation'] = degs['log2FoldChange'].apply(lambda x: 'up' if x > 0 else 'down')
        elif 'logFC' in degs.columns:
            degs['regulation'] = degs['logFC'].apply(lambda x: 'up' if x > 0 else 'down')

    print(f"\nDEG Dataset:")
    print(f"  Total DEGs: {len(degs)}")
    print(f"  Upregulated: {(degs['regulation'] == 'up').sum()}")
    print(f"  Downregulated: {(degs['regulation'] == 'down').sum()}")

    return degs


def match_degs_to_go_category(degs, go_genes_file, go_subcats_file, category_name):
    """
    Match DEGs to a GO category and its subcategories

    Parameters:
    -----------
    degs : DataFrame
        DEG data
    go_genes_file : Path
        Path to GO gene list CSV
    go_subcats_file : Path
        Path to GO subcategory mapping CSV
    category_name : str
        Category name for output files

    Returns:
    --------
    tuple: (degs_in_category, subcategory_stats)
    """
    print(f"\n{'-'*70}")
    print(f"MATCHING DEGs TO {category_name.upper()}")
    print(f"{'-'*70}")

    # Load GO annotations
    go_genes = pd.read_csv(go_genes_file)
    go_subcats = pd.read_csv(go_subcats_file)

    # Ensure uppercase
    go_genes['GeneSymbol'] = go_genes['GeneSymbol'].str.upper()
    go_subcats['GeneSymbol'] = go_subcats['GeneSymbol'].str.upper()

    # Match DEGs to category
    degs_in_cat = degs[degs['GeneSymbol'].isin(go_genes['GeneSymbol'])].copy()

    # Add subcategory annotations
    degs_in_cat = degs_in_cat.merge(
        go_subcats[['GeneSymbol', 'Subcategory', 'GO_Terms', 'Description']],
        on='GeneSymbol',
        how='left'
    )

    print(f"\n{category_name} matching:")
    print(f"  Total GO genes in dataset: {len(go_genes)}")
    print(f"  DEGs in {category_name}: {len(degs_in_cat)} / {len(go_genes)} ({len(degs_in_cat)/len(go_genes)*100:.1f}%)")
    print(f"  Upregulated: {(degs_in_cat['regulation'] == 'up').sum()}")
    print(f"  Downregulated: {(degs_in_cat['regulation'] == 'down').sum()}")

    return degs_in_cat, go_subcats


def fisher_exact_enrichment(k, M, K, N):
    """
    Fisher's exact test for enrichment (one-tailed, alternative='greater')

    Tests if DEGs are overrepresented in a category.

    Contingency table:
                  In Category    Not in Category
    DEG                k              K - k
    Not DEG          M - k      N - M - K + k

    Parameters:
    -----------
    k : int
        DEGs in category
    M : int
        Total genes in category
    K : int
        Total DEGs
    N : int
        Total genes

    Returns:
    --------
    tuple: (odds_ratio, p_value)
    """
    table = [
        [k, K - k],
        [M - k, N - M - K + k]
    ]

    odds_ratio, p_value = fisher_exact(table, alternative='greater')
    return odds_ratio, p_value


def fisher_exact_directional(a, b, c, d):
    """
    Fisher's exact test for directional bias (two-tailed)

    Tests if up/down ratio differs from global distribution.

    Contingency table:
                    This Category    Rest of DEGs
    Upregulated          a                c
    Downregulated        b                d

    Parameters:
    -----------
    a : int
        Up in this category
    b : int
        Down in this category
    c : int
        Up in rest of DEGs (TOTAL_UP_GLOBAL - a)
    d : int
        Down in rest of DEGs (TOTAL_DOWN_GLOBAL - b)

    Returns:
    --------
    tuple: (odds_ratio, p_value)
    """
    table = [
        [a, c],
        [b, d]
    ]

    odds_ratio, p_value = fisher_exact(table, alternative='two-sided')
    return odds_ratio, p_value


def calculate_subcategory_stats(degs_in_cat, go_subcats, category_name):
    """
    Calculate statistics for each subcategory

    Parameters:
    -----------
    degs_in_cat : DataFrame
        DEGs matched to category
    go_subcats : DataFrame
        Subcategory mappings
    category_name : str
        Category name

    Returns:
    --------
    DataFrame with statistical results for each subcategory
    """
    print(f"\nCalculating subcategory statistics for {category_name}...")

    # Get unique subcategories
    subcategories = go_subcats['Subcategory'].unique()

    results = []
    p_enrich_list = []
    p_dir_list = []

    for subcat in subcategories:
        # Get genes in this subcategory
        subcat_go_genes = go_subcats[go_subcats['Subcategory'] == subcat]
        M = len(subcat_go_genes)  # Total genes in subcategory

        # Match to DEGs
        subcat_degs = degs_in_cat[degs_in_cat['Subcategory'] == subcat]
        k = len(subcat_degs)  # DEGs in subcategory

        if k == 0:
            # No DEGs in this subcategory
            results.append({
                'Subcategory': subcat,
                'GO_Terms': subcat_go_genes.iloc[0]['GO_Terms'],
                'Description': subcat_go_genes.iloc[0]['Description'],
                'Total_Genes': M,
                'DEGs': 0,
                'Upregulated': 0,
                'Downregulated': 0,
                'Pct_Up': 0.0,
                'Pct_Down': 0.0,
                'Enrichment_OR': np.nan,
                'Enrichment_P': 1.0,
                'Directional_OR': np.nan,
                'Directional_P': 1.0
            })
            p_enrich_list.append(1.0)
            p_dir_list.append(1.0)
            continue

        # Count regulation directions
        n_up = (subcat_degs['regulation'] == 'up').sum()
        n_down = (subcat_degs['regulation'] == 'down').sum()

        # Enrichment test
        or_enrich, p_enrich = fisher_exact_enrichment(k, M, TOTAL_DEGS, TOTAL_GENES)
        p_enrich_list.append(p_enrich)

        # Directional bias test
        a = n_up
        b = n_down
        c = TOTAL_UP_GLOBAL - a
        d = TOTAL_DOWN_GLOBAL - b

        or_dir, p_dir = fisher_exact_directional(a, b, c, d)
        p_dir_list.append(p_dir)

        results.append({
            'Subcategory': subcat,
            'GO_Terms': subcat_go_genes.iloc[0]['GO_Terms'],
            'Description': subcat_go_genes.iloc[0]['Description'],
            'Total_Genes': M,
            'DEGs': k,
            'Upregulated': n_up,
            'Downregulated': n_down,
            'Pct_Up': n_up / k * 100 if k > 0 else 0.0,
            'Pct_Down': n_down / k * 100 if k > 0 else 0.0,
            'Enrichment_OR': or_enrich,
            'Enrichment_P': p_enrich,
            'Directional_OR': or_dir,
            'Directional_P': p_dir
        })

    # Create DataFrame
    stats_df = pd.DataFrame(results)

    # FDR correction - combine all p-values
    all_pvals = p_enrich_list + p_dir_list
    _, all_fdr, _, _ = multipletests(all_pvals, method='fdr_bh')

    # Split FDR values back
    stats_df['Enrichment_FDR'] = all_fdr[:len(p_enrich_list)]
    stats_df['Directional_FDR'] = all_fdr[len(p_enrich_list):]

    # Add significance markers
    def sig_marker(fdr):
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

    stats_df['Enrichment_Sig'] = stats_df['Enrichment_FDR'].apply(sig_marker)
    stats_df['Directional_Sig'] = stats_df['Directional_FDR'].apply(sig_marker)

    return stats_df


def create_summary_table(prot_stats, myo_stats):
    """Create overall summary table"""

    summary_rows = []

    # Overall proteolysis
    prot_total = prot_stats[['DEGs']].sum()['DEGs']
    prot_up = prot_stats[['Upregulated']].sum()['Upregulated']
    prot_down = prot_stats[['Downregulated']].sum()['Downregulated']

    summary_rows.append({
        'Category': 'Proteolysis',
        'GO_Term': 'GO:0006508',
        'Subcategory': 'ALL',
        'Total_DEGs': prot_total,
        'Upregulated': prot_up,
        'Downregulated': prot_down,
        'Pct_Up': prot_up / prot_total * 100 if prot_total > 0 else 0,
        'Pct_Down': prot_down / prot_total * 100 if prot_total > 0 else 0
    })

    # Proteolysis subcategories
    for _, row in prot_stats.iterrows():
        summary_rows.append({
            'Category': 'Proteolysis',
            'GO_Term': row['GO_Terms'],
            'Subcategory': row['Subcategory'],
            'Total_DEGs': row['DEGs'],
            'Upregulated': row['Upregulated'],
            'Downregulated': row['Downregulated'],
            'Pct_Up': row['Pct_Up'],
            'Pct_Down': row['Pct_Down']
        })

    # Overall myogenesis
    myo_total = myo_stats[['DEGs']].sum()['DEGs']
    myo_up = myo_stats[['Upregulated']].sum()['Upregulated']
    myo_down = myo_stats[['Downregulated']].sum()['Downregulated']

    summary_rows.append({
        'Category': 'Myogenesis',
        'GO_Term': 'GO:0007519',
        'Subcategory': 'ALL',
        'Total_DEGs': myo_total,
        'Upregulated': myo_up,
        'Downregulated': myo_down,
        'Pct_Up': myo_up / myo_total * 100 if myo_total > 0 else 0,
        'Pct_Down': myo_down / myo_total * 100 if myo_total > 0 else 0
    })

    # Myogenesis subcategories
    for _, row in myo_stats.iterrows():
        summary_rows.append({
            'Category': 'Myogenesis',
            'GO_Term': row['GO_Terms'],
            'Subcategory': row['Subcategory'],
            'Total_DEGs': row['DEGs'],
            'Upregulated': row['Upregulated'],
            'Downregulated': row['Downregulated'],
            'Pct_Up': row['Pct_Up'],
            'Pct_Down': row['Pct_Down']
        })

    return pd.DataFrame(summary_rows)


def main():
    """Main statistical analysis pipeline"""
    print("\n" + "="*70)
    print("GO-BASED STATISTICAL ANALYSIS")
    print("="*70)
    print("\nGlobal Dataset Constants:")
    print(f"  Total expressed genes: {TOTAL_GENES:,}")
    print(f"  Total DEGs: {TOTAL_DEGS:,}")
    print(f"  Upregulated: {TOTAL_UP_GLOBAL:,} ({TOTAL_UP_GLOBAL/TOTAL_DEGS*100:.1f}%)")
    print(f"  Downregulated: {TOTAL_DOWN_GLOBAL:,} ({TOTAL_DOWN_GLOBAL/TOTAL_DEGS*100:.1f}%)")

    # Load DEG data
    degs = load_deg_data()

    # Proteolysis analysis
    prot_degs, prot_subcats = match_degs_to_go_category(
        degs,
        RESULTS_DIR / "go_proteolysis_genes.csv",
        RESULTS_DIR / "go_proteolysis_subcats.csv",
        "Proteolysis"
    )

    prot_stats = calculate_subcategory_stats(prot_degs, prot_subcats, "Proteolysis")

    # Myogenesis analysis
    myo_degs, myo_subcats = match_degs_to_go_category(
        degs,
        RESULTS_DIR / "go_myogenesis_genes.csv",
        RESULTS_DIR / "go_myogenesis_subcats.csv",
        "Myogenesis"
    )

    myo_stats = calculate_subcategory_stats(myo_degs, myo_subcats, "Myogenesis")

    # ==== NEW: Expanded myogenesis analysis (multi-relationship GO terms) ====
    print("\n" + "="*70)
    print("EXPANDED MYOGENESIS ANALYSIS (Multi-relationship)")
    print("="*70)

    # List of GO term files generated by script 16
    go_term_files = [
        ('GO:0051145', 'GO_0051145'),  # smooth muscle cell differentiation
        ('GO:0051146', 'GO_0051146'),  # striated muscle cell differentiation
        ('GO:0055001', 'GO_0055001'),  # muscle cell development
        ('GO:0042693', 'GO_0042693'),  # muscle cell fate commitment
        ('GO:0014707', 'GO_0014707'),  # skeletal muscle tissue development (is_a)
        ('GO:0002074', 'GO_0002074'),  # skeletal muscle tissue development (is_a)
        ('GO:0035914', 'GO_0035914'),  # skeletal muscle satellite stem cell asymmetric division
        ('GO:0048741', 'GO_0048741'),  # skeletal muscle fiber development
        ('GO:0048630', 'GO_0048630'),  # skeletal muscle tissue growth
    ]

    expanded_stats = []

    for go_id, go_file_id in go_term_files:
        gene_file = RESULTS_DIR / f"go_{go_file_id}_genes.csv"

        if not gene_file.exists():
            print(f"  Warning: {gene_file.name} not found, skipping...")
            continue

        # Load GO term gene list
        go_genes = pd.read_csv(gene_file)
        go_name = go_genes['go_name'].iloc[0] if 'go_name' in go_genes.columns else go_id
        relationship = go_genes['relationship'].iloc[0] if 'relationship' in go_genes.columns else 'is_a'

        # Match to DEGs
        go_degs = degs[degs['GeneSymbol'].isin(go_genes['gene_symbol'])]

        # Calculate statistics
        M = len(go_genes)  # Total genes in this GO term
        k = len(go_degs)   # DEGs in this GO term

        if k == 0:
            print(f"  {go_id} - {go_name}: No DEGs found")
            continue

        n_up = (go_degs['regulation'] == 'up').sum()
        n_down = (go_degs['regulation'] == 'down').sum()

        # Enrichment test (one-tailed)
        # H0: DEGs not enriched in this GO term
        # Table: [[k, M-k], [TOTAL_DEGS-k, TOTAL_GENES-M-TOTAL_DEGS+k]]
        table_enrich = [
            [k, TOTAL_DEGS - k],
            [M - k, TOTAL_GENES - M - TOTAL_DEGS + k]
        ]
        or_enrich, p_enrich = fisher_exact(table_enrich, alternative='greater')

        # Directional bias test (two-tailed)
        # H0: up/down ratio same as global
        # Table: [[up_in_go, up_global-up_in_go], [down_in_go, down_global-down_in_go]]
        table_dir = [
            [n_up, TOTAL_UP_GLOBAL - n_up],
            [n_down, TOTAL_DOWN_GLOBAL - n_down]
        ]
        or_dir, p_dir = fisher_exact(table_dir, alternative='two-sided')

        expanded_stats.append({
            'go_id': go_id,
            'go_name': go_name,
            'relationship': relationship,
            'Total_Genes': M,
            'DEGs': k,
            'Upregulated': n_up,
            'Downregulated': n_down,
            'Pct_Up': (n_up / k * 100) if k > 0 else 0,
            'Pct_Down': (n_down / k * 100) if k > 0 else 0,
            'Enrichment_P': p_enrich,
            'Directional_P': p_dir
        })

        print(f"  {go_id} - {go_name}: {k} DEGs ({n_up}↑/{n_down}↓) [{relationship}]")

    # Create DataFrame
    expanded_df = pd.DataFrame(expanded_stats)

    # Apply FDR correction
    if len(expanded_df) > 0:
        all_p_values = list(expanded_df['Enrichment_P']) + list(expanded_df['Directional_P'])
        _, all_fdr, _, _ = multipletests(all_p_values, method='fdr_bh')

        n_terms = len(expanded_df)
        expanded_df['Enrichment_FDR'] = all_fdr[:n_terms]
        expanded_df['Directional_FDR'] = all_fdr[n_terms:]

        # Add significance markers
        def sig_marker(fdr):
            if fdr < 0.001:
                return '***'
            elif fdr < 0.01:
                return '**'
            elif fdr < 0.05:
                return '*'
            else:
                return 'ns'

        expanded_df['Enrichment_Sig'] = expanded_df['Enrichment_FDR'].apply(sig_marker)
        expanded_df['Directional_Sig'] = expanded_df['Directional_FDR'].apply(sig_marker)

        # Save expanded stats
        expanded_file = RESULTS_DIR / "go_myogenesis_stats_expanded.csv"
        expanded_df.to_csv(expanded_file, index=False)
        print(f"\nExpanded myogenesis statistics: {expanded_file}")
        print(f"  Total GO terms: {len(expanded_df)}")

    # Save results
    print("\n" + "="*70)
    print("SAVING RESULTS")
    print("="*70)

    # Save DEG lists
    prot_degs_file = RESULTS_DIR / "go_proteolysis_degs.csv"
    prot_degs.to_csv(prot_degs_file, index=False)
    print(f"\nProteolysis DEGs: {prot_degs_file}")
    print(f"  Total: {len(prot_degs)}")

    myo_degs_file = RESULTS_DIR / "go_myogenesis_degs.csv"
    myo_degs.to_csv(myo_degs_file, index=False)
    print(f"\nMyogenesis DEGs: {myo_degs_file}")
    print(f"  Total: {len(myo_degs)}")

    # Save statistical summaries
    prot_stats_file = RESULTS_DIR / "go_proteolysis_stats.csv"
    prot_stats.to_csv(prot_stats_file, index=False)
    print(f"\nProteolysis statistics: {prot_stats_file}")

    myo_stats_file = RESULTS_DIR / "go_myogenesis_stats.csv"
    myo_stats.to_csv(myo_stats_file, index=False)
    print(f"\nMyogenesis statistics: {myo_stats_file}")

    # Create and save summary table
    summary = create_summary_table(prot_stats, myo_stats)
    summary_file = RESULTS_DIR / "SUMMARY_GO_Analysis.csv"
    summary.to_csv(summary_file, index=False)
    print(f"\nSummary table: {summary_file}")

    # Display statistics
    print("\n" + "="*70)
    print("PROTEOLYSIS SUBCATEGORY STATISTICS")
    print("="*70)
    print(prot_stats[['Subcategory', 'DEGs', 'Upregulated', 'Downregulated',
                       'Enrichment_Sig', 'Directional_Sig']].to_string(index=False))

    print("\n" + "="*70)
    print("MYOGENESIS SUBCATEGORY STATISTICS")
    print("="*70)
    print(myo_stats[['Subcategory', 'DEGs', 'Upregulated', 'Downregulated',
                      'Enrichment_Sig', 'Directional_Sig']].to_string(index=False))

    # Summary
    print("\n" + "="*70)
    print("STATISTICAL ANALYSIS COMPLETE")
    print("="*70)
    print("\nOutputs:")
    print(f"  1. {prot_degs_file}")
    print(f"  2. {myo_degs_file}")
    print(f"  3. {prot_stats_file}")
    print(f"  4. {myo_stats_file}")
    print(f"  5. {summary_file}")
    print("\nNext steps:")
    print("  Run Script 18 to generate figures")


if __name__ == "__main__":
    main()
