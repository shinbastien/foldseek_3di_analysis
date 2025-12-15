#!/usr/bin/env python3
"""
Script 2: build_foldseek_dbs.py

Build Foldseek databases and run search.

This script:
1. Converts scope40_domains.tsv to FASTA (queries and targets)
2. Runs foldseek createdb
3. Runs foldseek search

Usage:
  python3 build_foldseek_dbs.py --queries selected_queries.tsv --scope40_domains scope40_domains.tsv
"""

import pandas as pd
from pathlib import Path
import subprocess
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(message)s')

def run_command(cmd, description):
    """Run shell command and log output."""
    logging.info(f"\n{description}...")
    logging.info(f"  Command: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logging.error(f"  ✗ Command failed!")
        logging.error(f"  stderr: {result.stderr}")
        return False
    else:
        logging.info(f"  ✓ Success")
        if result.stdout:
            for line in result.stdout.split('\n')[:5]:  # First 5 lines
                if line.strip():
                    logging.info(f"    {line}")
        return True

def main(args):
    logging.info("="*80)
    logging.info("Script 2: Build Foldseek Databases and Run Search")
    logging.info("="*80)
    
    pdb_dir = Path(args.pdb_dir)
    work_dir = Path(args.work_dir)
    work_dir.mkdir(exist_ok=True)
    
    # Load selected queries
    logging.info(f"\nLoading selected queries from: {args.queries}")
    queries_df = pd.read_csv(args.queries, sep='\t')
    query_ids = queries_df['domain_id'].tolist()
    logging.info(f"  Loaded {len(query_ids)} query domains")
    
    # Load all domains (targets)
    logging.info(f"\nLoading all domains from: {args.scope40_domains}")
    all_domains_df = pd.read_csv(args.scope40_domains, sep='\t')
    target_ids = all_domains_df['domain_id'].tolist()
    logging.info(f"  Loaded {len(target_ids)} target domains")
    
    # Build databases (use PDB directory directly)
    logging.info(f"\nBuilding Foldseek databases...")
    logging.info(f"  Note: Using PDB directory directly for both query and target")
    
    queryDB = work_dir / "queryDB"
    cmd = ["foldseek", "createdb", str(pdb_dir), str(queryDB)]
    if not run_command(cmd, "Creating query DB from PDB directory"):
        logging.error("Failed to create query DB")
        return False
    
    targetDB = work_dir / "targetDB"
    # For target DB, we can reuse queryDB or create a separate one
    # For efficiency, we'll create it the same way
    cmd = ["foldseek", "createdb", str(pdb_dir), str(targetDB)]
    if not run_command(cmd, "Creating target DB from PDB directory"):
        logging.error("Failed to create target DB")
        return False
    
    # Run search
    logging.info(f"\nRunning Foldseek search...")
    out = work_dir / "foldseek_out"
    tmp = work_dir / "tmp"
    
    cmd = ["foldseek", "search", str(queryDB), str(targetDB), str(out), str(tmp)]
    logging.info(f"  Command: {' '.join(cmd)}")
    logging.info(f"  (This may take a while...)")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"  ✗ Foldseek search failed!")
        logging.error(f"  stderr: {result.stderr}")
        return False
    else:
        logging.info(f"  ✓ Search completed")
    
    # Convert output to TSV
    logging.info(f"\nConverting Foldseek output to TSV...")
    out_tsv = work_dir / "foldseek_results.tsv"
    cmd = ["foldseek", "convertalis", str(queryDB), str(targetDB), str(out), str(out_tsv), "--format-output", "query,target,qlen,tlen,evalue,bits"]
    
    if not run_command(cmd, "Converting to TSV"):
        logging.error("Failed to convert output")
        return False
    
    logging.info("\n" + "="*80)
    logging.info("✓ Foldseek databases built and search completed!")
    logging.info(f"  Results: {out_tsv}")
    logging.info("="*80)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build Foldseek databases and run search"
    )
    parser.add_argument(
        "--queries",
        type=str,
        default="selected_queries.tsv",
        help="Path to selected queries TSV"
    )
    parser.add_argument(
        "--scope40_domains",
        type=str,
        default="scope40_domains.tsv",
        help="Path to all scope40 domains TSV"
    )
    parser.add_argument(
        "--pdb_dir",
        type=str,
        default="scope_pdb",
        help="Path to PDB directory"
    )
    parser.add_argument(
        "--work_dir",
        type=str,
        default="foldseek_work",
        help="Working directory for Foldseek outputs"
    )
    
    args = parser.parse_args()
    main(args)
