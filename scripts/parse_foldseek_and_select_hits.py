#!/usr/bin/env python3
"""
Script 3: parse_foldseek_and_select_hits.py

Parse Foldseek output and select candidate Non-CP homolog pairs.

Reads:
  - foldseek_results.tsv: Foldseek search output (query, target, evalue, bitscore)
  - scope40_domains.tsv: Domain metadata (SCCS, species)
  - selected_queries.tsv: Selected query domains

Filters:
  - Remove self-hits (query == target)
  - Keep only same superfamily: sccs(query) == sccs(target)
  - Optionally different species: species_id(query) != species_id(target)
  - Top K hits per query (ranked by bitscore)

Outputs:
  - candidate_noncp_pairs.tsv: Candidate pairs with metadata
"""

import pandas as pd
from pathlib import Path
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(message)s')

# Default parameters
DEFAULT_K_HITS = 2
DEFAULT_REQUIRE_DIFFERENT_SPECIES = True

def main(args):
    logging.info("="*80)
    logging.info("Script 3: Parse Foldseek Results and Select Candidate Pairs")
    logging.info("="*80)
    
    foldseek_file = Path(args.foldseek_results)
    domains_file = Path(args.scope40_domains)
    queries_file = Path(args.selected_queries)
    output_file = Path(args.output)
    k_hits = args.k_hits
    require_diff_species = args.require_different_species
    
    # Load Foldseek results
    logging.info(f"\nLoading Foldseek results from: {foldseek_file}")
    foldseek_df = pd.read_csv(foldseek_file, sep='\t', header=None, 
                               names=['query', 'target', 'qlen', 'tlen', 'evalue', 'bits'])
    logging.info(f"  Loaded {len(foldseek_df)} hits")
    logging.info(f"  Columns: {foldseek_df.columns.tolist()}")
    
    # Load domain metadata
    logging.info(f"\nLoading domain metadata from: {domains_file}")
    domains_df = pd.read_csv(domains_file, sep='\t')
    domain_dict = {}
    for _, row in domains_df.iterrows():
        domain_dict[row['domain_id']] = {
            'sccs': row['sccs'],
            'superfamily': row['superfamily'],
            'species_id': row['species_id'],
            'length': row['length']
        }
    logging.info(f"  Loaded {len(domain_dict)} domains")
    
    # Load selected queries
    logging.info(f"\nLoading selected queries from: {queries_file}")
    queries_df = pd.read_csv(queries_file, sep='\t')
    query_set = set(queries_df['domain_id'].values)
    logging.info(f"  Loaded {len(query_set)} queries")
    
    # Filter and process hits
    logging.info(f"\nProcessing Foldseek hits...")
    
    # Filter: only selected queries
    filtered_df = foldseek_df[foldseek_df['query'].isin(query_set)].copy()
    logging.info(f"  After filtering to selected queries: {len(filtered_df)} hits")
    
    # Remove self-hits
    filtered_df = filtered_df[filtered_df['query'] != filtered_df['target']]
    logging.info(f"  After removing self-hits: {len(filtered_df)} hits")
    
    # Add metadata
    logging.info(f"  Adding metadata...")
    filtered_df['query_sccs'] = filtered_df['query'].map(lambda x: domain_dict.get(x, {}).get('sccs', None))
    filtered_df['target_sccs'] = filtered_df['target'].map(lambda x: domain_dict.get(x, {}).get('sccs', None))
    filtered_df['query_superfamily'] = filtered_df['query'].map(lambda x: domain_dict.get(x, {}).get('superfamily', None))
    filtered_df['target_superfamily'] = filtered_df['target'].map(lambda x: domain_dict.get(x, {}).get('superfamily', None))
    filtered_df['query_species'] = filtered_df['query'].map(lambda x: domain_dict.get(x, {}).get('species_id', None))
    filtered_df['target_species'] = filtered_df['target'].map(lambda x: domain_dict.get(x, {}).get('species_id', None))
    filtered_df['query_length'] = filtered_df['query'].map(lambda x: domain_dict.get(x, {}).get('length', None))
    filtered_df['target_length'] = filtered_df['target'].map(lambda x: domain_dict.get(x, {}).get('length', None))
    
    # Filter: same superfamily
    filtered_df = filtered_df[filtered_df['query_superfamily'] == filtered_df['target_superfamily']]
    logging.info(f"  After same-SF filter: {len(filtered_df)} hits")
    
    # Filter: different species (optional)
    if require_diff_species:
        filtered_df = filtered_df[filtered_df['query_species'] != filtered_df['target_species']]
        logging.info(f"  After different-species filter: {len(filtered_df)} hits")
    
    # Sort by bitscore (descending) and take top K per query
    logging.info(f"  Selecting top {k_hits} hits per query...")
    # Sort by bits (descending, since higher is better for bitscore)
    filtered_df = filtered_df.sort_values('bits', ascending=False)
    
    # Group by query and take top K
    candidate_list = []
    for query_id, group in filtered_df.groupby('query'):
        top_k = group.head(k_hits)
        candidate_list.append(top_k)
    
    if candidate_list:
        candidate_df = pd.concat(candidate_list, ignore_index=True)
    else:
        candidate_df = pd.DataFrame()
    
    logging.info(f"  ✓ Selected {len(candidate_df)} candidate pairs")
    
    # Prepare output columns
    output_df = candidate_df[[
        'query', 'target',
        'query_sccs', 'target_sccs',
        'query_species', 'target_species',
        'query_length', 'target_length',
        'bits', 'evalue'
    ]].copy()
    
    output_df.columns = [
        'query_id', 'target_id',
        'query_sccs', 'target_sccs',
        'query_species', 'target_species',
        'query_length', 'target_length',
        'foldseek_bitscore', 'foldseek_evalue'
    ]
    
    # Save
    logging.info(f"\nSaving to: {output_file}")
    output_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"  ✓ Saved {len(output_df)} pairs")
    
    # Summary
    logging.info("\n" + "-"*80)
    logging.info("Summary:")
    logging.info(f"  Total candidate pairs: {len(output_df)}")
    logging.info(f"  Unique queries: {len(output_df['query_id'].unique())}")
    logging.info(f"  Pairs per query (mean): {len(output_df) / len(output_df['query_id'].unique()):.1f}")
    logging.info(f"  Bitscore range: {output_df['foldseek_bitscore'].min():.1f} - {output_df['foldseek_bitscore'].max():.1f}")
    logging.info("-"*80)
    
    logging.info("\n" + "="*80)
    logging.info("✓ Foldseek parsing completed!")
    logging.info("="*80)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse Foldseek results and select candidate Non-CP pairs"
    )
    parser.add_argument(
        "--foldseek_results",
        type=str,
        default="foldseek_work/foldseek_results.tsv",
        help="Path to Foldseek results TSV"
    )
    parser.add_argument(
        "--scope40_domains",
        type=str,
        default="scope40_domains.tsv",
        help="Path to scope40_domains.tsv"
    )
    parser.add_argument(
        "--selected_queries",
        type=str,
        default="selected_queries.tsv",
        help="Path to selected queries TSV"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="candidate_noncp_pairs.tsv",
        help="Output file for candidate pairs"
    )
    parser.add_argument(
        "--k_hits",
        type=int,
        default=DEFAULT_K_HITS,
        help=f"Maximum hits per query (default: {DEFAULT_K_HITS})"
    )
    parser.add_argument(
        "--require_different_species",
        type=bool,
        default=DEFAULT_REQUIRE_DIFFERENT_SPECIES,
        help=f"Require query and target from different species (default: {DEFAULT_REQUIRE_DIFFERENT_SPECIES})"
    )
    
    args = parser.parse_args()
    main(args)
