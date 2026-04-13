#!/usr/bin/env python3
"""
Script 16: Fully Automated GO Tree-Based Gene Annotation

Uses GO ontology tree structure to automatically identify proteolysis
and myogenesis gene subcategories based on GO hierarchy.

Approach:
- Loads GO ontology (go-basic.obo) and queries tree for child terms
- Uses QuickGO REST API to get gene-to-GO annotations from EBI database
- Maps genes to GO child terms automatically using official GO annotations
- NO manual gene lists - everything derived from GO tree + QuickGO database

Primary GO Terms (only hardcoded values):
1. GO:0006508 - Proteolysis
2. GO:0042692 - Muscle cell differentiation

Data Sources:
- GO Ontology: Gene Ontology Consortium
- Gene Annotations: QuickGO (EBI) - http://www.ebi.ac.uk/QuickGO

Author: Analysis pipeline
Date: 2025-12-04
"""

import pandas as pd
import numpy as np
import pronto
import requests
import json
import io
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Setup paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
GO_DIR = DATA_DIR / "go_ontology"
RESULTS_DIR = BASE_DIR / "results" / "tables"

# GO Terms (ONLY hardcoded values)
GO_PROTEOLYSIS = 'GO:0006508'
GO_MYOGENESIS = 'GO:0042692'
GO_SKELETAL_MUSCLE = 'GO:0007519'  # Skeletal muscle tissue development

# QuickGO API endpoint
QUICKGO_URL = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch"


def load_go_ontology():
    """Load GO ontology from OBO file"""
    print("\n" + "="*70)
    print("LOADING GO ONTOLOGY")
    print("="*70)

    obo_file = GO_DIR / "go-basic.obo"
    if not obo_file.exists():
        raise FileNotFoundError(f"GO ontology file not found: {obo_file}")

    print(f"\nLoading: {obo_file}")
    ont = pronto.Ontology(str(obo_file))

    print(f"✓ Loaded {len(ont)} GO terms")

    return ont


def query_quickgo_annotations(taxon_id='10090', timeout=180):
    """
    Query QuickGO REST API for mouse GO annotations

    Parameters:
    -----------
    taxon_id : str
        NCBI taxonomy ID (10090 for Mus musculus)
    timeout : int
        Request timeout in seconds

    Returns:
    --------
    pandas.DataFrame with columns: GENE_PRODUCT_ID, SYMBOL, GO_TERM, GO_NAME, GO_ASPECT, EVIDENCE_CODE
    """
    print("\n" + "="*70)
    print("QUERYING QUICKGO FOR MOUSE GO ANNOTATIONS")
    print("="*70)

    print(f"\nTaxon ID: {taxon_id} (Mus musculus)")
    print("Querying QuickGO API...")
    print("This may take 1-2 minutes for large datasets...")

    # Build query parameters - simplified to avoid API issues
    # QuickGO returns default columns including symbol, GO ID, GO name, etc.
    params = {
        'taxonId': taxon_id,
        'aspect': 'biological_process',  # Only Biological Process
        'downloadLimit': 2000000  # Max limit (may return 100k due to API limitation)
    }

    try:
        # Query API
        response = requests.get(
            QUICKGO_URL,
            params=params,
            headers={'Accept': 'text/tsv'},
            timeout=timeout
        )

        response.raise_for_status()  # Raise error for bad status codes

        # Parse TSV response
        print("✓ Received response, parsing...")
        annotations_df = pd.read_csv(io.StringIO(response.text), sep='\t', comment='!')

        # Standardize column names
        column_map = {
            'GENE PRODUCT ID': 'GENE_PRODUCT_ID',
            'SYMBOL': 'SYMBOL',
            'GO TERM': 'GO_TERM',
            'GO NAME': 'GO_NAME',
            'GO ASPECT': 'GO_ASPECT',
            'EVIDENCE CODE': 'EVIDENCE_CODE'
        }
        annotations_df.rename(columns=column_map, inplace=True)

        print(f"✓ Retrieved {len(annotations_df)} annotations")
        print(f"✓ Unique genes: {annotations_df['SYMBOL'].nunique()}")
        print(f"✓ Unique GO terms: {annotations_df['GO_TERM'].nunique()}")

        return annotations_df

    except requests.exceptions.RequestException as e:
        print(f"ERROR: QuickGO API request failed: {e}")
        raise
    except Exception as e:
        print(f"ERROR: Failed to parse QuickGO response: {e}")
        raise


def create_gene_to_go_mapping(annotations_df):
    """
    Create gene-to-GO mapping from QuickGO annotations

    Returns:
    --------
    dict: {gene_symbol: set(GO_IDs)}
    """
    print("\n" + "="*70)
    print("CREATING GENE-TO-GO MAPPING")
    print("="*70)

    gene_to_go = defaultdict(set)

    for _, row in annotations_df.iterrows():
        symbol = str(row['SYMBOL']).upper()
        go_term = row['GO_TERM']

        if pd.notna(symbol) and pd.notna(go_term):
            gene_to_go[symbol].add(go_term)

    print(f"✓ Mapped {len(gene_to_go)} genes to GO terms")

    return dict(gene_to_go)


def get_all_descendants(ont, parent_go_id, include_self=False):
    """
    Get all descendant terms of a GO term

    Returns:
    --------
    set of GO IDs
    """
    if parent_go_id not in ont:
        raise ValueError(f"GO term not found: {parent_go_id}")

    parent_term = ont[parent_go_id]
    descendants = set()

    for term in parent_term.subclasses(with_self=include_self):
        descendants.add(term.id)

    return descendants


def get_direct_children(ont, parent_go_id):
    """
    Get direct children (distance=1) of a GO term

    Returns:
    --------
    list of tuples: [(go_id, go_name, go_def), ...]
    """
    if parent_go_id not in ont:
        raise ValueError(f"GO term not found: {parent_go_id}")

    parent_term = ont[parent_go_id]
    direct_children = []

    for term in parent_term.subclasses(distance=1, with_self=False):
        direct_children.append((
            term.id,
            term.name,
            term.definition if hasattr(term, 'definition') else ''
        ))

    return sorted(direct_children, key=lambda x: x[1])


def get_children_by_relationship(ont, parent_go_id, relationship_type='is_a'):
    """
    Get direct children of a GO term by specific relationship type.

    Parameters:
    -----------
    ont : pronto.Ontology
        Pronto Ontology object
    parent_go_id : str
        Parent GO term ID (e.g., "GO:0042692")
    relationship_type : str
        Relationship type: 'is_a', 'part_of', 'capable_of_part_of', etc.

    Returns:
    --------
    list of tuples: [(child_id, child_name, relationship_type, definition), ...]
    """
    if parent_go_id not in ont:
        raise ValueError(f"GO term not found: {parent_go_id}")

    parent_term = ont[parent_go_id]
    children = []

    if relationship_type == 'is_a':
        # Use existing .subclasses() method for is_a relationships
        for term in parent_term.subclasses(distance=1, with_self=False):
            children.append((
                term.id,
                term.name,
                'is_a',
                term.definition if hasattr(term, 'definition') else ''
            ))
    else:
        # For part_of and other relationships, iterate through all terms
        # Find terms that have the parent in their relationships
        for term in ont.terms():
            if hasattr(term, 'relationships') and term.relationships:
                for rel, targets in term.relationships.items():
                    # Check if this relationship matches and parent is in targets
                    if rel.id == relationship_type and parent_term in targets:
                        children.append((
                            term.id,
                            term.name,
                            relationship_type,
                            term.definition if hasattr(term, 'definition') else ''
                        ))
                        break

    return sorted(children, key=lambda x: x[1])


def map_genes_to_go_children(ont, parent_go_id, gene_to_go, min_genes=3):
    """
    Map genes to GO child terms using GO tree and QuickGO annotations

    Returns:
    --------
    dict: {go_id: {'name': ..., 'genes': [...], 'definition': ...}}
    """
    print(f"\n  Mapping genes to children of {parent_go_id}...")

    # Get direct children
    direct_children = get_direct_children(ont, parent_go_id)
    print(f"    Direct children: {len(direct_children)}")

    # For each direct child, find genes annotated to it or its descendants
    child_mappings = {}

    for child_id, child_name, child_def in direct_children:
        # Get all descendants of this child (including itself)
        child_descendants = get_all_descendants(ont, child_id, include_self=True)

        # Find genes annotated to this child or any of its descendants
        genes_in_child = set()
        for gene, go_ids in gene_to_go.items():
            if child_descendants & go_ids:  # Intersection
                genes_in_child.add(gene)

        if len(genes_in_child) >= min_genes:
            child_mappings[child_id] = {
                'name': child_name,
                'definition': child_def,
                'genes': sorted(list(genes_in_child))
            }

    print(f"    Subcategories with ≥{min_genes} genes: {len(child_mappings)}")

    return child_mappings


def map_genes_to_go_children_multi_relationship(ont, parent_go_id, gene_to_go, relationship_types=['is_a', 'part_of'], min_genes=3):
    """
    Map genes to GO children using multiple relationship types.

    Parameters:
    -----------
    ont : pronto.Ontology
        Pronto Ontology object
    parent_go_id : str
        Parent GO term ID
    gene_to_go : dict
        Mapping of gene symbols to sets of GO IDs
    relationship_types : list
        List of relationship types to consider (e.g., ['is_a', 'part_of'])
    min_genes : int
        Minimum number of genes required to include a subcategory

    Returns:
    --------
    dict: {go_id: {'name': ..., 'relationship': ..., 'definition': ..., 'genes': [...]}}
    """
    print(f"\n  Mapping genes to children of {parent_go_id} (multi-relationship)...")
    print(f"    Relationship types: {relationship_types}")

    child_mappings = {}

    # Get children for each relationship type
    for rel_type in relationship_types:
        children = get_children_by_relationship(ont, parent_go_id, rel_type)
        print(f"    {rel_type} children: {len(children)}")

        for child_id, child_name, relationship, child_def in children:
            # Get all descendants of this child (including itself)
            child_descendants = get_all_descendants(ont, child_id, include_self=True)

            # Find genes annotated to this subtree
            genes_in_child = set()
            for gene, go_ids in gene_to_go.items():
                if child_descendants & go_ids:
                    genes_in_child.add(gene)

            # Filter by minimum gene threshold
            if len(genes_in_child) >= min_genes:
                child_mappings[child_id] = {
                    'name': child_name,
                    'relationship': relationship,
                    'definition': child_def,
                    'genes': sorted(list(genes_in_child))
                }

    print(f"    Total subcategories with ≥{min_genes} genes: {len(child_mappings)}")

    return child_mappings


def create_gene_annotations(ont, gene_to_go):
    """
    Create comprehensive GO-based gene lists using tree structure

    Returns:
    --------
    tuple: (proteolysis_df, myogenesis_df, proteolysis_subcats, myogenesis_subcats)
    """
    print("\n" + "="*70)
    print("CREATING GO TREE-BASED GENE LISTS")
    print("="*70)

    # Map genes to proteolysis children
    print(f"\nParent term: {GO_PROTEOLYSIS}")
    prot_parent = ont[GO_PROTEOLYSIS]
    print(f"  Name: {prot_parent.name}")

    prot_mappings = map_genes_to_go_children(ont, GO_PROTEOLYSIS, gene_to_go, min_genes=3)

    # Map genes to myogenesis children
    print(f"\nParent term: {GO_MYOGENESIS}")
    myo_parent = ont[GO_MYOGENESIS]
    print(f"  Name: {myo_parent.name}")

    myo_mappings = map_genes_to_go_children(ont, GO_MYOGENESIS, gene_to_go, min_genes=3)

    # Create proteolysis annotations
    print("\n" + "="*70)
    print("PROTEOLYSIS SUBCATEGORIES:")
    print("="*70)

    prot_rows = []
    for go_id, mapping in prot_mappings.items():
        print(f"  {go_id} - {mapping['name']}: {len(mapping['genes'])} genes")

        for gene in mapping['genes']:
            prot_rows.append({
                'GeneSymbol': gene,
                'Subcategory': mapping['name'],
                'GO_Terms': go_id,
                'Description': mapping['definition']
            })

    prot_subcat_df = pd.DataFrame(prot_rows)

    # Create complete proteolysis gene list
    all_prot_genes = sorted(list(set([row['GeneSymbol'] for row in prot_rows])))
    prot_df = pd.DataFrame({
        'GeneSymbol': all_prot_genes,
        'GO_Term': GO_PROTEOLYSIS,
        'GO_Name': prot_parent.name
    })

    print(f"\nTotal unique proteolysis genes: {len(all_prot_genes)}")

    # Create myogenesis annotations
    print("\n" + "="*70)
    print("MYOGENESIS SUBCATEGORIES:")
    print("="*70)

    myo_rows = []
    for go_id, mapping in myo_mappings.items():
        print(f"  {go_id} - {mapping['name']}: {len(mapping['genes'])} genes")

        for gene in mapping['genes']:
            myo_rows.append({
                'GeneSymbol': gene,
                'Subcategory': mapping['name'],
                'GO_Terms': go_id,
                'Description': mapping['definition']
            })

    myo_subcat_df = pd.DataFrame(myo_rows)

    # Create complete myogenesis gene list
    all_myo_genes = sorted(list(set([row['GeneSymbol'] for row in myo_rows])))
    myo_df = pd.DataFrame({
        'GeneSymbol': all_myo_genes,
        'GO_Term': GO_MYOGENESIS,
        'GO_Name': myo_parent.name
    })

    print(f"\nTotal unique myogenesis genes: {len(all_myo_genes)}")

    # ==== NEW: Process myogenesis with multi-relationship approach ====
    print("\n" + "="*70)
    print("EXPANDED MYOGENESIS ANALYSIS (Multi-relationship)")
    print("="*70)

    # Process GO:0042692 (muscle cell differentiation) with both is_a and part_of
    print(f"\nParent term 1: {GO_MYOGENESIS}")
    myo_parent_1 = ont[GO_MYOGENESIS]
    print(f"  Name: {myo_parent_1.name}")
    myo_children_1 = map_genes_to_go_children_multi_relationship(
        ont,
        GO_MYOGENESIS,
        gene_to_go,
        relationship_types=['is_a', 'part_of'],
        min_genes=3
    )

    # Process GO:0007519 (skeletal muscle tissue development)
    print(f"\nParent term 2: {GO_SKELETAL_MUSCLE}")
    myo_parent_2 = ont[GO_SKELETAL_MUSCLE]
    print(f"  Name: {myo_parent_2.name}")
    myo_children_2 = map_genes_to_go_children_multi_relationship(
        ont,
        GO_SKELETAL_MUSCLE,
        gene_to_go,
        relationship_types=['is_a', 'part_of'],
        min_genes=3
    )

    # Combine all children
    all_myo_children = {**myo_children_1, **myo_children_2}

    # Save individual GO term gene lists with relationship metadata
    print("\nSaving individual GO term gene lists...")
    for child_id, data in all_myo_children.items():
        output_file = RESULTS_DIR / f'go_{child_id.replace(":", "_")}_genes.csv'
        gene_df = pd.DataFrame({
            'gene_symbol': data['genes'],
            'go_id': child_id,
            'go_name': data['name'],
            'relationship': data['relationship']
        })
        gene_df.to_csv(output_file, index=False)
        print(f"  {child_id}: {len(data['genes'])} genes ({data['relationship']}) -> {output_file.name}")

    print(f"\nTotal GO terms processed: {len(all_myo_children)}")

    return prot_df, myo_df, prot_subcat_df, myo_subcat_df


def save_results(prot_df, myo_df, prot_subcat_df, myo_subcat_df, ont):
    """Save all results to CSV and JSON"""

    print("\n" + "="*70)
    print("SAVING RESULTS")
    print("="*70)

    # Save main gene lists
    prot_file = RESULTS_DIR / "go_proteolysis_genes.csv"
    prot_df.to_csv(prot_file, index=False)
    print(f"\nProteolysis genes: {prot_file}")
    print(f"  Total: {len(prot_df)}")

    myo_file = RESULTS_DIR / "go_myogenesis_genes.csv"
    myo_df.to_csv(myo_file, index=False)
    print(f"\nMyogenesis genes: {myo_file}")
    print(f"  Total: {len(myo_df)}")

    # Save subcategory mappings
    prot_subcat_file = RESULTS_DIR / "go_proteolysis_subcats.csv"
    prot_subcat_df.to_csv(prot_subcat_file, index=False)
    print(f"\nProteolysis subcategories: {prot_subcat_file}")
    print(f"  Subcategories: {prot_subcat_df['Subcategory'].nunique()}")

    myo_subcat_file = RESULTS_DIR / "go_myogenesis_subcats.csv"
    myo_subcat_df.to_csv(myo_subcat_file, index=False)
    print(f"\nMyogenesis subcategories: {myo_subcat_file}")
    print(f"  Subcategories: {myo_subcat_df['Subcategory'].nunique()}")

    # Save metadata
    metadata = {
        'analysis_name': 'Fully Automated GO Tree-Based Gene Annotation',
        'creation_date': datetime.now().isoformat(),
        'annotation_approach': 'GO ontology tree structure + QuickGO REST API annotations',
        'method': 'Fully automated - no manual gene lists',
        'go_ontology': {
            'file': 'go-basic.obo',
            'source': 'http://purl.obolibrary.org/obo/go/go-basic.obo',
            'download_date': datetime.now().isoformat(),
            'format_version': 'OBO 1.2',
            'total_terms': len(ont)
        },
        'annotation_source': {
            'tool': 'QuickGO REST API',
            'organism': 'Mus musculus (taxonId=10090)',
            'database': 'GO:BP (Biological Process)',
            'url': QUICKGO_URL,
            'query_date': datetime.now().isoformat()
        },
        'go_terms': {
            'proteolysis': {
                'go_id': GO_PROTEOLYSIS,
                'go_name': ont[GO_PROTEOLYSIS].name,
                'total_genes': len(prot_df),
                'subcategories': prot_subcat_df['Subcategory'].nunique(),
                'direct_children': len(list(ont[GO_PROTEOLYSIS].subclasses(distance=1, with_self=False)))
            },
            'myogenesis': {
                'go_id': GO_MYOGENESIS,
                'go_name': ont[GO_MYOGENESIS].name,
                'total_genes': len(myo_df),
                'subcategories': myo_subcat_df['Subcategory'].nunique(),
                'direct_children': len(list(ont[GO_MYOGENESIS].subclasses(distance=1, with_self=False)))
            }
        },
        'filtering': {
            'min_genes_per_subcategory': 3,
            'note': 'Only direct children with ≥3 annotated genes included'
        },
        'subcategory_definitions': {
            'proteolysis': {
                row['Subcategory']: {
                    'go_terms': row['GO_Terms'],
                    'description': row['Description'],
                    'gene_count': int((prot_subcat_df['Subcategory'] == row['Subcategory']).sum())
                }
                for _, row in prot_subcat_df[['Subcategory', 'GO_Terms', 'Description']].drop_duplicates().iterrows()
            },
            'myogenesis': {
                row['Subcategory']: {
                    'go_terms': row['GO_Terms'],
                    'description': row['Description'],
                    'gene_count': int((myo_subcat_df['Subcategory'] == row['Subcategory']).sum())
                }
                for _, row in myo_subcat_df[['Subcategory', 'GO_Terms', 'Description']].drop_duplicates().iterrows()
            }
        }
    }

    metadata_file = RESULTS_DIR / "GO_Analysis_Metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"\nMetadata: {metadata_file}")

    return metadata


def main():
    """Main pipeline"""
    print("\n" + "="*70)
    print("FULLY AUTOMATED GO TREE-BASED GENE ANNOTATION")
    print("="*70)
    print(f"\nHardcoded values (only 3):")
    print(f"  1. {GO_PROTEOLYSIS} - Proteolysis")
    print(f"  2. {GO_MYOGENESIS} - Muscle cell differentiation")
    print(f"  3. {GO_SKELETAL_MUSCLE} - Skeletal muscle tissue development")
    print("\nAll genes and subcategories derived from:")
    print("  - GO ontology tree structure (direct children)")
    print("  - QuickGO REST API (EBI) for gene-GO annotations")
    print("  - Multiple relationship types (is_a, part_of)")
    print("  - Minimum 3 genes per subcategory filter")

    # Load GO ontology
    ont = load_go_ontology()

    # Query QuickGO for mouse GO annotations
    annotations_df = query_quickgo_annotations(taxon_id='10090')

    # Create gene-to-GO mapping
    gene_to_go = create_gene_to_go_mapping(annotations_df)

    # Create annotations using GO tree structure
    prot_df, myo_df, prot_subcat_df, myo_subcat_df = create_gene_annotations(ont, gene_to_go)

    # Save results
    metadata = save_results(prot_df, myo_df, prot_subcat_df, myo_subcat_df, ont)

    # Summary
    print("\n" + "="*70)
    print("GO ANNOTATION COMPLETE")
    print("="*70)
    print("\nSummary:")
    print(f"  Proteolysis: {len(prot_df)} genes across {prot_subcat_df['Subcategory'].nunique()} GO subcategories")
    print(f"  Myogenesis: {len(myo_df)} genes across {myo_subcat_df['Subcategory'].nunique()} GO subcategories")
    print("\nSubcategories are GO term names from tree structure")
    print("Genes mapped automatically via QuickGO annotations")
    print("\nOutputs:")
    print(f"  1. {RESULTS_DIR / 'go_proteolysis_genes.csv'}")
    print(f"  2. {RESULTS_DIR / 'go_myogenesis_genes.csv'}")
    print(f"  3. {RESULTS_DIR / 'go_proteolysis_subcats.csv'}")
    print(f"  4. {RESULTS_DIR / 'go_myogenesis_subcats.csv'}")
    print(f"  5. {RESULTS_DIR / 'GO_Analysis_Metadata.json'}")
    print("\nNext steps:")
    print("  Run Script 17 to match DEGs and perform statistical analysis")


if __name__ == "__main__":
    main()
