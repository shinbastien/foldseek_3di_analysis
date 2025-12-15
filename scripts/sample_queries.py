#!/usr/bin/env python3
"""
Script 1: sample_queries.py

Select query domains for Non-CP homolog discovery using stratified sampling.

Reads:
  - scope40_domains.tsv: All SCOPe40 domains with length, SCCS, species info
    - input_data/datasets/cp_positive_pairs.tsv: CP-positive pairs (to exclude these domains)

Performs:
  - Exclude domains appearing in CP-positive pairs
  - Stratified sampling over length bins to match CP-positive length distribution
  
Outputs:
  - selected_queries.tsv: Selected query domains with metadata
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(message)s')

# Default parameters
DEFAULT_N_QUERY = 2000

# CP-positive length distribution (from analysis)
CP_LENGTH_BINS = {
    '<50': 0.000,
    '50-100': 0.366,
    '100-150': 0.350,
    '150-200': 0.099,
    '200-300': 0.136,
    '300-500': 0.050,
    '>500': 0.000,
}

def define_length_bin(length):
    """Assign length to a bin."""
    if length < 50:
        return '<50'
    elif length < 100:
        return '50-100'
    elif length < 150:
        return '100-150'
    elif length < 200:
        return '150-200'
    elif length < 300:
        return '200-300'
    elif length < 500:
        return '300-500'
    else:
        return '>500'

def main(args):
    logging.info("="*80)
    logging.info("Script 1: Sample Query Domains for Non-CP Set")
    logging.info("="*80)
    
    domains_file = Path(args.scope40_domains)
    cp_pairs_file = Path(args.cp_positive_pairs)
    output_file = Path(args.output)
    n_query = args.n_query
    
    # Load domains
    logging.info(f"\nLoading domains from: {domains_file}")
    domains_df = pd.read_csv(domains_file, sep='\t')
    logging.info(f"  Loaded {len(domains_df)} domains")
    
    # Load CP-positive pairs
    logging.info(f"\nLoading CP-positive pairs from: {cp_pairs_file}")
    cp_pairs_df = pd.read_csv(cp_pairs_file, sep='\t')
    logging.info(f"  Loaded {len(cp_pairs_df)} pairs")
    
    # Build set of CP domains to exclude
    cp_domains = set()
    cp_domains.update(cp_pairs_df['domA'].values)
    cp_domains.update(cp_pairs_df['domB'].values)
    logging.info(f"  Unique domains in CP set: {len(cp_domains)}")
    
    # Filter out CP domains
    logging.info(f"\nFiltering out CP domains...")
    candidate_domains = domains_df[~domains_df['domain_id'].isin(cp_domains)].copy()
    logging.info(f"  Candidate query domains: {len(candidate_domains)}")
    
    # Assign length bins
    candidate_domains['length_bin'] = candidate_domains['length'].apply(define_length_bin)
    
    # Stratified sampling
    logging.info(f"\nPerforming stratified sampling (target N={n_query})...")
    logging.info(f"  Target distribution based on CP-positive set:")
    for bin_name, proportion in CP_LENGTH_BINS.items():
        if proportion > 0:
            count = int(n_query * proportion)
            logging.info(f"    {bin_name}: {proportion*100:.1f}% ({count} domains)")
    
    sampled_list = []
    for bin_name, proportion in CP_LENGTH_BINS.items():
        if proportion == 0:
            continue
        
        target_count = int(n_query * proportion)
        bin_domains = candidate_domains[candidate_domains['length_bin'] == bin_name]
        
        if len(bin_domains) == 0:
            logging.warning(f"  No domains in bin {bin_name}")
            continue
        
        if len(bin_domains) < target_count:
            # Sample all available
            sampled = bin_domains
            logging.warning(f"  {bin_name}: Requested {target_count}, but only {len(bin_domains)} available. Using all.")
        else:
            # Random sampling without replacement
            sampled = bin_domains.sample(n=target_count, random_state=42)
        
        sampled_list.append(sampled)
    
    selected_df = pd.concat(sampled_list, ignore_index=True)
    selected_df = selected_df.drop(columns=['length_bin'])
    
    logging.info(f"\n✓ Selected {len(selected_df)} query domains")
    
    # Check length distribution
    logging.info(f"\nActual selected distribution:")
    selected_df['_bin'] = selected_df['length'].apply(define_length_bin)
    for bin_name in ['<50', '50-100', '100-150', '150-200', '200-300', '300-500', '>500']:
        count = len(selected_df[selected_df['_bin'] == bin_name])
        if count > 0:
            pct = 100 * count / len(selected_df)
            logging.info(f"  {bin_name}: {count} ({pct:.1f}%)")
    selected_df = selected_df.drop(columns=['_bin'])
    
    # Save
    logging.info(f"\nSaving to: {output_file}")
    selected_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"  ✓ Saved {len(selected_df)} domains")
    
    logging.info("\n" + "="*80)
    logging.info("✓ Query sampling completed!")
    logging.info("="*80)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sample query domains for Non-CP homolog discovery"
    )
    parser.add_argument(
        "--scope40_domains",
        type=str,
        default="scope40_domains.tsv",
        help="Path to scope40_domains.tsv"
    )
    parser.add_argument(
        "--cp_positive_pairs",
        type=str,
        default="input_data/datasets/cp_positive_pairs.tsv",
        help="Path to CP-positive pairs TSV (default: input_data/datasets/cp_positive_pairs.tsv)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="selected_queries.tsv",
        help="Output file for selected queries"
    )
    parser.add_argument(
        "--n_query",
        type=int,
        default=DEFAULT_N_QUERY,
        help=f"Number of query domains to select (default: {DEFAULT_N_QUERY})"
    )
    
    args = parser.parse_args()
    main(args)
