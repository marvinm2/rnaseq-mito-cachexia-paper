#!/usr/bin/env python3
"""
Script 03: Comprehensive Pathway Analysis (3 approaches)
1. MitoCarta pathway hierarchy (unbiased)
2. Statistical enrichment on MitoCarta pathways
3. g:Profiler enrichment (fully unbiased)

Author: Analysis pipeline
Date: 2025-10-06
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# Setup paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results" / "tables"

def parse_mitocarta_pathway_hierarchy(pathway_str):
    """Parse MitoCarta pathway string into hierarchical levels"""
    if pd.isna(pathway_str) or str(pathway_str) == '0':
        return []

    # Split by pipe (multiple pathways)
    pathways = str(pathway_str).split('|')

    all_levels = []
    for pathway in pathways:
        pathway = pathway.strip()
        if pathway and pathway != '0':
            # Split by > to get hierarchy
            levels = [level.strip() for level in pathway.split('>')]
            all_levels.append(levels)

    return all_levels

def extract_pathway_at_level(pathway_hierarchies, level=0):
    """Extract pathway names at a specific hierarchical level"""
    pathways = []
    if not pathway_hierarchies:
        return pathways

    # pathway_hierarchies is a list of lists of hierarchical levels
    # e.g. [[['OXPHOS', 'Complex I', 'CI subunits'], ['OXPHOS', 'OXPHOS subunits']]]
    for hierarchy in pathway_hierarchies:
        # hierarchy is a list like ['OXPHOS', 'Complex I', 'CI subunits']
        if isinstance(hierarchy, list) and len(hierarchy) > level:
            pathways.append(hierarchy[level])

    return pathways

def analyze_mitocarta_pathways_unbiased(mito_degs, mitocarta_all):
    """
    APPROACH 1: Analyze MitoCarta pathways WITHOUT bias
    Count DEGs at each level of the pathway hierarchy
    """
    print("\n" + "="*70)
    print("APPROACH 1: MitoCarta Pathway Hierarchy (UNBIASED)")
    print("="*70)

    # Parse all pathways for DEGs
    mito_degs['pathway_hierarchy'] = mito_degs['MitoCarta3.0_MitoPathways'].apply(
        parse_mitocarta_pathway_hierarchy
    )

    # Count pathways at each level
    results_by_level = {}

    for level in range(3):  # Level 0, 1, 2
        print(f"\n--- Hierarchy Level {level} ---")

        # Extract pathways at this level
        all_pathways = []
        up_pathways = []
        down_pathways = []

        for idx, row in mito_degs.iterrows():
            pathways_at_level = extract_pathway_at_level(row['pathway_hierarchy'], level)
            all_pathways.extend(pathways_at_level)

            if row['regulation'] == 'up':
                up_pathways.extend(pathways_at_level)
            else:
                down_pathways.extend(pathways_at_level)

        # Count occurrences
        pathway_counts = Counter(all_pathways)
        up_counts = Counter(up_pathways)
        down_counts = Counter(down_pathways)

        # Create summary (ALL pathways, not top 20 — downstream figure scripts
        # need the full set)
        summary = []
        for pathway, total in pathway_counts.most_common():
            up = up_counts.get(pathway, 0)
            down = down_counts.get(pathway, 0)
            pct_down = (down / total * 100) if total > 0 else 0

            summary.append({
                'level': level,
                'pathway': pathway,
                'total_degs': total,
                'upregulated': up,
                'downregulated': down,
                'pct_downregulated': pct_down
            })

        results_by_level[f'level_{level}'] = pd.DataFrame(summary)

        # Print top 10
        print(pd.DataFrame(summary).head(10).to_string(index=False))

    return results_by_level

def calculate_enrichment_statistics(mito_degs, mitocarta_all):
    """
    APPROACH 2: Statistical enrichment analysis
    Test if specific pathways are over-represented in DEGs
    """
    print("\n" + "="*70)
    print("APPROACH 2: Statistical Enrichment (Fisher's Exact Test)")
    print("="*70)

    # Parse pathways for all mitochondrial genes
    mitocarta_all['pathway_hierarchy'] = mitocarta_all['MitoCarta3.0_MitoPathways'].apply(
        parse_mitocarta_pathway_hierarchy
    )

    # Get DEG gene symbols
    deg_symbols = set(mito_degs['GeneSymbol'])

    # Extract level 1 pathways (most informative level)
    all_mito_pathways = {}
    for idx, row in mitocarta_all.iterrows():
        gene = str(row['Symbol'])
        pathways = extract_pathway_at_level(row['pathway_hierarchy'], level=1)
        if pathways:
            all_mito_pathways[gene] = pathways

    # Count genes per pathway
    pathway_gene_map = {}
    for gene, pathways in all_mito_pathways.items():
        for pathway in pathways:
            if pathway not in pathway_gene_map:
                pathway_gene_map[pathway] = set()
            pathway_gene_map[pathway].add(gene)

    # Perform Fisher's exact test for each pathway
    enrichment_results = []

    total_mito_genes = len(all_mito_pathways)
    total_deg_genes = len(deg_symbols)

    print(f"\nBackground: {total_mito_genes} mitochondrial genes")
    print(f"DEGs: {total_deg_genes} mitochondrial DEGs")
    print(f"\nTesting {len(pathway_gene_map)} pathways...\n")

    for pathway, pathway_genes in pathway_gene_map.items():
        # 2x2 contingency table
        deg_in_pathway = len(pathway_genes & deg_symbols)
        deg_not_in_pathway = total_deg_genes - deg_in_pathway
        not_deg_in_pathway = len(pathway_genes) - deg_in_pathway
        not_deg_not_in_pathway = total_mito_genes - len(pathway_genes) - deg_not_in_pathway

        # Fisher's exact test
        odds_ratio, p_value = fisher_exact([
            [deg_in_pathway, deg_not_in_pathway],
            [not_deg_in_pathway, not_deg_not_in_pathway]
        ])

        # Calculate fold enrichment
        expected = (len(pathway_genes) / total_mito_genes) * total_deg_genes
        fold_enrichment = deg_in_pathway / expected if expected > 0 else 0

        enrichment_results.append({
            'pathway': pathway,
            'total_genes_in_pathway': len(pathway_genes),
            'degs_in_pathway': deg_in_pathway,
            'expected': expected,
            'fold_enrichment': fold_enrichment,
            'odds_ratio': odds_ratio,
            'p_value': p_value
        })

    enrichment_df = pd.DataFrame(enrichment_results)
    enrichment_df = enrichment_df.sort_values('p_value')

    # Apply Bonferroni correction
    enrichment_df['p_adj'] = enrichment_df['p_value'] * len(enrichment_df)
    enrichment_df['p_adj'] = enrichment_df['p_adj'].clip(upper=1.0)

    # Filter significant
    significant = enrichment_df[enrichment_df['p_adj'] < 0.05]

    print(f"Significant enrichments (p_adj < 0.05): {len(significant)}")
    print("\n" + significant.head(15).to_string(index=False))

    return enrichment_df

def run_gprofiler_analysis(mito_degs):
    """
    APPROACH 3: g:Profiler enrichment (fully unbiased)
    Submit mitochondrial DEG symbols to g:Profiler via web API
    """
    print("\n" + "="*70)
    print("APPROACH 3: g:Profiler Enrichment (UNBIASED)")
    print("="*70)

    try:
        import requests

        # Get gene symbols
        gene_list = mito_degs['GeneSymbol'].dropna().unique().tolist()
        print(f"\nSubmitting {len(gene_list)} genes to g:Profiler web API...")

        # Call g:Profiler API directly
        url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

        payload = {
            "organism": "mmusculus",
            "query": gene_list,
            "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"],
            "user_threshold": 0.05,
            "significance_threshold_method": "fdr",
            "ordered": False,
            "exclude_iea": False
        }

        response = requests.post(url, json=payload)

        if response.status_code == 200:
            data = response.json()

            if 'result' in data and data['result']:
                results_list = data['result']

                # Convert to DataFrame
                results = pd.DataFrame(results_list)
                print(f"\nFound {len(results)} enriched terms (FDR < 0.05)")

                # Filter for mitochondrial-related terms
                mito_keywords = ['mitochondri', 'OXPHOS', 'oxidative phosphorylation',
                               'electron transport', 'respiratory chain', 'ATP synth',
                               'TCA cycle', 'citric acid', 'respiratory']

                if 'name' in results.columns:
                    mito_results = results[
                        results['name'].str.contains('|'.join(mito_keywords), case=False, na=False)
                    ]
                elif 'description' in results.columns:
                    mito_results = results[
                        results['description'].str.contains('|'.join(mito_keywords), case=False, na=False)
                    ]
                else:
                    mito_results = pd.DataFrame()

                print(f"\nMitochondrial-related terms: {len(mito_results)}")
                print("\nTop 20 enriched terms:")

                # Display key columns
                display_cols = []
                for col in ['source', 'native', 'name', 'description', 'p_value',
                           'intersection_size', 'term_size', 'precision', 'recall']:
                    if col in results.columns:
                        display_cols.append(col)

                if display_cols:
                    print(results[display_cols].head(20).to_string(index=False, max_colwidth=60))

                print("\n\nMitochondrial-specific terms:")
                if len(mito_results) > 0:
                    print(mito_results[display_cols].to_string(index=False, max_colwidth=60))
                else:
                    print("No terms matched mitochondrial keywords")

                return results
            else:
                print("\nNo significant enrichments found")
                return pd.DataFrame()
        else:
            raise Exception(f"API request failed with status code {response.status_code}")

    except Exception as e:
        print(f"\n⚠️  Error running g:Profiler API: {e}")
        print("\nAlternative: Use web interface at https://biit.cs.ut.ee/gprofiler/")

        # Save gene list for manual submission
        gene_list_file = RESULTS_DIR / "mitochondrial_genes_for_gprofiler.txt"
        gene_list = mito_degs['GeneSymbol'].dropna().unique()
        with open(gene_list_file, 'w') as f:
            f.write('\n'.join(gene_list))
        print(f"\nGene list saved to: {gene_list_file}")
        print(f"Number of genes: {len(gene_list)}")
        print("\nInstructions:")
        print("1. Go to https://biit.cs.ut.ee/gprofiler/gost")
        print("2. Copy the gene list from the file above")
        print("3. Select organism: Mus musculus")
        print("4. Run enrichment analysis")

        return None

def main():
    """Main analysis pipeline"""
    print("\n" + "="*70)
    print("COMPREHENSIVE MITOCHONDRIAL PATHWAY ANALYSIS")
    print("="*70)

    # Load data
    print("\nLoading data...")
    mito_degs = pd.read_csv(RESULTS_DIR / "mitochondrial_degs.csv")
    mitocarta_all = pd.read_csv(DATA_DIR / "mitocarta.csv")

    print(f"Mitochondrial DEGs: {len(mito_degs)}")
    print(f"Total MitoCarta genes: {len(mitocarta_all)}")

    # APPROACH 1: Unbiased MitoCarta hierarchy
    hierarchy_results = analyze_mitocarta_pathways_unbiased(mito_degs, mitocarta_all)

    # Save hierarchy results
    for level_name, df in hierarchy_results.items():
        output_file = RESULTS_DIR / f"mitocarta_pathways_{level_name}.csv"
        df.to_csv(output_file, index=False)
        print(f"\nSaved {level_name} results to: {output_file}")

    # Also emit level 2 in the "detailed" schema that downstream figure scripts
    # (07_generate_figure3_with_statistics.py, mitocarta_pathway_gene_extraction.py)
    # read: columns = level2, down, up, total, pct_down, pct_up.
    level2_df = hierarchy_results['level_2'].copy()
    detailed = pd.DataFrame({
        'level2': level2_df['pathway'],
        'down': level2_df['downregulated'],
        'up': level2_df['upregulated'],
        'total': level2_df['total_degs'],
        'pct_down': level2_df['pct_downregulated'].round(1),
        'pct_up': (100 - level2_df['pct_downregulated']).round(1),
    })
    detailed_file = RESULTS_DIR / "mitocarta_pathways_level_2_detailed.csv"
    detailed.to_csv(detailed_file, index=False)
    print(f"Saved level_2_detailed to: {detailed_file}")

    # APPROACH 2: Statistical enrichment
    enrichment_results = calculate_enrichment_statistics(mito_degs, mitocarta_all)

    # Save enrichment results
    enrichment_file = RESULTS_DIR / "mitocarta_pathway_enrichment.csv"
    enrichment_results.to_csv(enrichment_file, index=False)
    print(f"\nSaved enrichment results to: {enrichment_file}")

    # APPROACH 3: g:Profiler
    gprofiler_results = run_gprofiler_analysis(mito_degs)

    if gprofiler_results is not None and len(gprofiler_results) > 0:
        gprofiler_file = RESULTS_DIR / "gprofiler_results.csv"
        gprofiler_results.to_csv(gprofiler_file, index=False)
        print(f"\nSaved g:Profiler results to: {gprofiler_file}")

    print("\n" + "="*70)
    print("COMPREHENSIVE PATHWAY ANALYSIS COMPLETE")
    print("="*70)
    print("\nResults saved:")
    print("  1. MitoCarta hierarchy: mitocarta_pathways_level_*.csv")
    print("  2. Enrichment statistics: mitocarta_pathway_enrichment.csv")
    print("  3. g:Profiler: gprofiler_results.csv")

if __name__ == "__main__":
    main()
