#!/usr/bin/env python3
"""
Script 4: run_tmalign_and_filter_noncp.py

Run TM-align on candidate pairs and filter for high-confidence Non-CP homologs.

Reads:
  - candidate_noncp_pairs.tsv: Candidate pairs from Foldseek
  - PDB files: scope_pdb/<domain_id>.pdb

Runs:
  - TM-align A.pdb B.pdb (normal)
  - TM-align A.pdb B.pdb -cp (with CP option)
  - Compute ΔCP = TM_cp - TM_nocp

Filters:
  - TM_nocp >= TM_THRESHOLD (default 0.5)
  - |ΔCP| <= DELTA_CP_THRESHOLD (default 0.02)

Downsamples:
  - To N_FINAL pairs (default 2000)
  - Max MAX_PAIRS_PER_QUERY per query (default 2)

Outputs:
  - noncp_with_tm.tsv: All pairs with TM-align results
  - noncp_homolog_pairs.tsv: Filtered final set
"""

import pandas as pd
import subprocess
from pathlib import Path
import logging
import argparse
import re
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(message)s')

# Default thresholds
DEFAULT_TM_THRESHOLD = 0.5
DEFAULT_DELTA_CP_THRESHOLD = 0.02
DEFAULT_N_FINAL = 2000
DEFAULT_MAX_PAIRS_PER_QUERY = 2

def parse_tm_score(tmalign_output):
    """Parse TM-score from TM-align output."""
    # Look for "TM-score" line
    for line in tmalign_output.split('\n'):
        if 'TM-score' in line and '=' in line:
            # Format: "TM-score= X (if normalized by L_max)"
            match = re.search(r'TM-score=\s*([\d.]+)', line)
            if match:
                return float(match.group(1))
    return None

def run_tmalign(pdb1, pdb2, cp_mode=False):
    """
    Run TM-align and return TM-score.
    
    Args:
        pdb1, pdb2: Path to PDB files
        cp_mode: If True, use -cp option
    
    Returns:
        TM-score (float) or None if failed
    """
    try:
        # Use full path to TMalign
        tmalign_path = "/home/jugipalace/miniconda3/envs/torchclean/bin/TMalign"
        cmd = [tmalign_path, str(pdb1), str(pdb2)]
        if cp_mode:
            cmd.append("-cp")
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode != 0:
            return None
        
        tm_score = parse_tm_score(result.stdout)
        return tm_score
    
    except subprocess.TimeoutExpired:
        logging.warning(f"  TM-align timeout for {pdb1} vs {pdb2}")
        return None
    except Exception as e:
        logging.warning(f"  TM-align error for {pdb1} vs {pdb2}: {e}")
        return None

def main(args):
    logging.info("="*80)
    logging.info("Script 4: Run TM-align and Filter Non-CP Homolog Pairs")
    logging.info("="*80)
    
    candidates_file = Path(args.candidate_pairs)
    pdb_dir = Path(args.pdb_dir)
    output_tm_file = Path(args.output_with_tm)
    output_final_file = Path(args.output_final)
    
    tm_threshold = args.tm_threshold
    delta_cp_threshold = args.delta_cp_threshold
    n_final = args.n_final
    max_pairs_per_query = args.max_pairs_per_query
    
    logging.info(f"\nThresholds:")
    logging.info(f"  TM_nocp >= {tm_threshold}")
    logging.info(f"  |ΔCP| <= {delta_cp_threshold}")
    logging.info(f"  Final target: {n_final} pairs")
    logging.info(f"  Max pairs per query: {max_pairs_per_query}")
    
    # Load candidates
    logging.info(f"\nLoading candidate pairs from: {candidates_file}")
    candidates_df = pd.read_csv(candidates_file, sep='\t')
    logging.info(f"  Loaded {len(candidates_df)} candidates")
    
    # Run TM-align
    logging.info(f"\nRunning TM-align on {len(candidates_df)} pairs...")
    logging.info(f"  (This will take a while...)")
    
    tm_results = []
    failed_count = 0
    
    for idx, row in candidates_df.iterrows():
        query_id = row['query_id']
        target_id = row['target_id']
        
        pdb1 = pdb_dir / f"{query_id}.pdb"
        pdb2 = pdb_dir / f"{target_id}.pdb"
        
        if not pdb1.exists() or not pdb2.exists():
            failed_count += 1
            continue
        
        # Run normal TM-align
        tm_nocp = run_tmalign(pdb1, pdb2, cp_mode=False)
        if tm_nocp is None:
            failed_count += 1
            continue
        
        # Run CP TM-align
        tm_cp = run_tmalign(pdb1, pdb2, cp_mode=True)
        if tm_cp is None:
            failed_count += 1
            continue
        
        # Compute ΔCP
        delta_cp = tm_cp - tm_nocp
        
        tm_results.append({
            'query_id': query_id,
            'target_id': target_id,
            'TM_nocp': tm_nocp,
            'TM_cp': tm_cp,
            'ΔCP': delta_cp,
            'foldseek_bitscore': row.get('foldseek_bitscore', None),
            'foldseek_evalue': row.get('foldseek_evalue', None),
            'query_sccs': row.get('query_sccs', None),
            'target_sccs': row.get('target_sccs', None),
            'query_species': row.get('query_species', None),
            'target_species': row.get('target_species', None),
            'query_length': row.get('query_length', None),
            'target_length': row.get('target_length', None),
        })
        
        if (idx + 1) % 100 == 0:
            logging.info(f"  Processed {idx + 1}/{len(candidates_df)} pairs...")
    
    tm_df = pd.DataFrame(tm_results)
    logging.info(f"\n  ✓ Completed {len(tm_df)} pairs successfully")
    logging.info(f"  Failed: {failed_count} pairs")
    
    # Save all TM results
    logging.info(f"\nSaving all TM results to: {output_tm_file}")
    tm_df.to_csv(output_tm_file, sep='\t', index=False)
    logging.info(f"  ✓ Saved {len(tm_df)} pairs")
    
    # Apply filters
    logging.info(f"\nApplying filters...")
    logging.info(f"  Filter 1: TM_nocp >= {tm_threshold}")
    filtered_df = tm_df[tm_df['TM_nocp'] >= tm_threshold].copy()
    logging.info(f"    Remaining: {len(filtered_df)} pairs")
    
    logging.info(f"  Filter 2: |ΔCP| <= {delta_cp_threshold}")
    filtered_df = filtered_df[filtered_df['ΔCP'].abs() <= delta_cp_threshold].copy()
    logging.info(f"    Remaining: {len(filtered_df)} pairs")
    
    # Check length distribution
    logging.info(f"\n  Length distribution of filtered pairs:")
    for min_len, max_len, label in [
        (50, 100, '50-100'), (100, 150, '100-150'), (150, 200, '150-200'),
        (200, 300, '200-300'), (300, 500, '300-500')
    ]:
        count = len(filtered_df[
            (filtered_df['query_length'] >= min_len) & (filtered_df['query_length'] < max_len)
        ])
        if count > 0:
            pct = 100 * count / len(filtered_df)
            logging.info(f"    {label}: {count} ({pct:.1f}%)")
    
    # Downsample if needed
    if len(filtered_df) > n_final:
        logging.info(f"\nDownsampling to {n_final} pairs...")
        
        # Limit pairs per query
        final_list = []
        query_counts = {}
        
        for _, row in filtered_df.iterrows():
            query_id = row['query_id']
            
            if query_id not in query_counts:
                query_counts[query_id] = 0
            
            if query_counts[query_id] < max_pairs_per_query:
                final_list.append(row)
                query_counts[query_id] += 1
            
            if len(final_list) >= n_final:
                break
        
        final_df = pd.DataFrame(final_list)
        logging.info(f"  ✓ Downsampled to {len(final_df)} pairs")
    else:
        final_df = filtered_df
        logging.info(f"  No downsampling needed ({len(final_df)} < {n_final})")
    
    # Save final set
    logging.info(f"\nSaving final Non-CP homolog set to: {output_final_file}")
    final_df.to_csv(output_final_file, sep='\t', index=False)
    logging.info(f"  ✓ Saved {len(final_df)} pairs")
    
    # Final summary
    logging.info("\n" + "="*80)
    logging.info("Final Non-CP Homolog Set Summary")
    logging.info("="*80)
    logging.info(f"Total pairs in final set: {len(final_df)}")
    logging.info(f"Unique queries: {len(final_df['query_id'].unique())}")
    logging.info(f"Mean TM_nocp: {final_df['TM_nocp'].mean():.3f} (±{final_df['TM_nocp'].std():.3f})")
    logging.info(f"Mean ΔCP: {final_df['ΔCP'].mean():.4f} (±{final_df['ΔCP'].std():.4f})")
    logging.info(f"Mean query length: {final_df['query_length'].mean():.0f}")
    logging.info(f"Mean target length: {final_df['target_length'].mean():.0f}")
    logging.info("="*80)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run TM-align and filter Non-CP homolog pairs"
    )
    parser.add_argument(
        "--candidate_pairs",
        type=str,
        default="candidate_noncp_pairs.tsv",
        help="Path to candidate pairs TSV"
    )
    parser.add_argument(
        "--pdb_dir",
        type=str,
        default="scope_pdb",
        help="Path to PDB directory"
    )
    parser.add_argument(
        "--output_with_tm",
        type=str,
        default="noncp_with_tm.tsv",
        help="Output file with all TM-align results"
    )
    parser.add_argument(
        "--output_final",
        type=str,
        default="noncp_homolog_pairs.tsv",
        help="Output file for final Non-CP homolog set"
    )
    parser.add_argument(
        "--tm_threshold",
        type=float,
        default=DEFAULT_TM_THRESHOLD,
        help=f"Minimum TM_nocp score (default: {DEFAULT_TM_THRESHOLD})"
    )
    parser.add_argument(
        "--delta_cp_threshold",
        type=float,
        default=DEFAULT_DELTA_CP_THRESHOLD,
        help=f"Maximum |ΔCP| (default: {DEFAULT_DELTA_CP_THRESHOLD})"
    )
    parser.add_argument(
        "--n_final",
        type=int,
        default=DEFAULT_N_FINAL,
        help=f"Target number of final pairs (default: {DEFAULT_N_FINAL})"
    )
    parser.add_argument(
        "--max_pairs_per_query",
        type=int,
        default=DEFAULT_MAX_PAIRS_PER_QUERY,
        help=f"Max pairs per query in final set (default: {DEFAULT_MAX_PAIRS_PER_QUERY})"
    )
    
    args = parser.parse_args()
    main(args)
