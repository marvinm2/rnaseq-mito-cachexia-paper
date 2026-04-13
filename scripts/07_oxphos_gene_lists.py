#!/usr/bin/env python3
"""
MitoCarta OXPHOS Hierarchy Gene Annotation
==========================================

This script parses the MitoCarta 3.0 pathway annotations to extract genes
at all three hierarchical levels of OXPHOS (Oxidative Phosphorylation):
- Level 0: OXPHOS parent (all OXPHOS genes)
- Level 1: Direct children (OXPHOS subunits, assembly factors, Complex I/IV/V)
- Level 2: Grandchildren (CI subunits, CI assembly, CIV subunits, etc.)

The MitoCarta3.0_MitoPathways column uses:
- "|" to separate multiple pathways for a single gene
- ">" to separate hierarchical levels within a pathway

Example: "OXPHOS > Complex I > CI subunits | OXPHOS > OXPHOS subunits"

Output:
-------
Gene lists saved as CSV files for each OXPHOS category at all levels.
Files saved to: results/tables/oxphos_genes/

Author: Bioinformatics Team
Date: December 2025
"""

import pandas as pd
from pathlib import Path
from collections import defaultdict


def parse_mitocarta_pathway_hierarchy(pathway_string):
    """
    Parse MitoCarta pathway string into hierarchical levels.

    Parameters
    ----------
    pathway_string : str
        MitoCarta pathway annotation with "|" and ">" separators
        Example: "OXPHOS > Complex I > CI subunits | OXPHOS > OXPHOS subunits"

    Returns
    -------
    dict
        Dictionary with keys 'level_0', 'level_1', 'level_2'
        Each contains a list of category names at that level
    """
    if pd.isna(pathway_string):
        return {'level_0': [], 'level_1': [], 'level_2': []}

    pathways = pathway_string.split('|')
    levels = {'level_0': [], 'level_1': [], 'level_2': []}

    for pathway in pathways:
        parts = [p.strip() for p in pathway.split('>')]

        if len(parts) >= 1:
            levels['level_0'].append(parts[0])
        if len(parts) >= 2:
            levels['level_1'].append(parts[1])
        if len(parts) >= 3:
            levels['level_2'].append(parts[2])

    return levels


def extract_oxphos_genes(mitocarta_df, min_genes_l2=3):
    """
    Extract genes for all OXPHOS hierarchy levels.

    Parameters
    ----------
    mitocarta_df : pd.DataFrame
        MitoCarta 3.0 data with 'Symbol' and 'MitoCarta3.0_MitoPathways' columns
    min_genes_l2 : int
        Minimum number of genes required for level 2 categories (default: 3)

    Returns
    -------
    dict
        Nested dictionary with structure:
        {
            'level_0': {'OXPHOS': [gene_list]},
            'level_1': {
                'OXPHOS subunits': [gene_list],
                'Complex I': [gene_list],
                ...
            },
            'level_2': {
                'CI subunits': [gene_list],
                'CI assembly factors': [gene_list],
                ...
            }
        }
    """
    # Filter to genes with OXPHOS in their pathway annotation
    oxphos_genes = mitocarta_df[
        mitocarta_df['MitoCarta3.0_MitoPathways'].str.contains('OXPHOS', na=False)
    ].copy()

    print(f"Found {len(oxphos_genes)} genes with OXPHOS pathway annotations")

    # Initialize hierarchy storage
    hierarchy = {
        'level_0': defaultdict(set),
        'level_1': defaultdict(set),
        'level_2': defaultdict(set)
    }

    # Parse each gene's pathway annotations
    for idx, row in oxphos_genes.iterrows():
        gene = row['Symbol']
        pathway_str = row['MitoCarta3.0_MitoPathways']

        levels = parse_mitocarta_pathway_hierarchy(pathway_str)

        # Add gene to each level's categories
        for level_key in ['level_0', 'level_1', 'level_2']:
            for category in levels[level_key]:
                # Only include if the category is OXPHOS-related
                if level_key == 'level_0':
                    if category == 'OXPHOS':
                        hierarchy[level_key][category].add(gene)
                else:
                    # For L1 and L2, the gene must have OXPHOS as L0 parent
                    if 'OXPHOS' in levels['level_0']:
                        hierarchy[level_key][category].add(gene)

    # Filter L2 by minimum gene threshold
    print(f"\nLevel 2 categories before filtering (min_genes={min_genes_l2}):")
    for category, genes in sorted(hierarchy['level_2'].items()):
        print(f"  {category}: {len(genes)} genes")

    hierarchy['level_2'] = {
        k: v for k, v in hierarchy['level_2'].items()
        if len(v) >= min_genes_l2
    }

    # Convert sets to sorted lists
    for level in hierarchy.values():
        for category in level:
            level[category] = sorted(list(level[category]))

    return hierarchy


def save_oxphos_gene_lists(hierarchy, output_dir):
    """
    Save gene lists for each OXPHOS category to CSV files.

    Parameters
    ----------
    hierarchy : dict
        Nested dictionary from extract_oxphos_genes()
    output_dir : Path
        Directory to save output CSV files

    Returns
    -------
    dict
        Summary statistics for each level
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = {}

    for level_name, categories in hierarchy.items():
        level_files = []

        for category_name, genes in categories.items():
            # Create safe filename
            safe_name = category_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
            filename = f"oxphos_{level_name}_{safe_name}_genes.csv"
            filepath = output_dir / filename

            # Create DataFrame with gene list
            df = pd.DataFrame({
                'gene_symbol': genes,
                'level': level_name,
                'category': category_name,
                'n_genes': len(genes)
            })

            # Save to CSV
            df.to_csv(filepath, index=False)
            level_files.append(filename)

            print(f"Saved {len(genes)} genes to {filename}")

        summary[level_name] = {
            'n_categories': len(categories),
            'files': level_files
        }

    return summary


def print_hierarchy_summary(hierarchy):
    """
    Print summary statistics for the OXPHOS hierarchy.

    Parameters
    ----------
    hierarchy : dict
        Nested dictionary from extract_oxphos_genes()
    """
    print("\n" + "="*70)
    print("OXPHOS Hierarchy Summary")
    print("="*70)

    for level_name in ['level_0', 'level_1', 'level_2']:
        categories = hierarchy[level_name]
        print(f"\n{level_name.upper().replace('_', ' ')}:")
        print(f"  Number of categories: {len(categories)}")

        for category, genes in sorted(categories.items(), key=lambda x: -len(x[1])):
            print(f"    {category}: {len(genes)} genes")

    # Total unique genes across all levels
    all_genes = set()
    for level in hierarchy.values():
        for genes in level.values():
            all_genes.update(genes)

    print(f"\nTotal unique OXPHOS genes: {len(all_genes)}")
    print("="*70)


def main():
    """Main execution function."""

    # Define paths
    base_dir = Path(__file__).parent.parent
    mitocarta_file = base_dir / "data" / "mitocarta.csv"
    output_dir = base_dir / "results" / "tables" / "oxphos_genes"

    print("="*70)
    print("MitoCarta OXPHOS Hierarchy Gene Annotation")
    print("="*70)
    print(f"\nInput file: {mitocarta_file}")
    print(f"Output directory: {output_dir}")

    # Load MitoCarta data
    print("\nLoading MitoCarta 3.0 data...")
    mitocarta_df = pd.read_csv(mitocarta_file)
    print(f"Loaded {len(mitocarta_df)} total mitochondrial genes")

    # Extract OXPHOS genes at all hierarchy levels
    print("\nExtracting OXPHOS hierarchy...")
    hierarchy = extract_oxphos_genes(mitocarta_df, min_genes_l2=3)

    # Print summary
    print_hierarchy_summary(hierarchy)

    # Save gene lists
    print(f"\nSaving gene lists to {output_dir}...")
    summary = save_oxphos_gene_lists(hierarchy, output_dir)

    print(f"\n{'='*70}")
    print("OXPHOS gene annotation complete!")
    print(f"{'='*70}")
    print(f"\nTotal files created:")
    for level_name, info in summary.items():
        print(f"  {level_name}: {info['n_categories']} files")

    print(f"\nAll gene lists saved to: {output_dir}/")


if __name__ == "__main__":
    main()
