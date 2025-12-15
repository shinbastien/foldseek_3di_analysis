#!/usr/bin/env python3
"""
Build the CP candidate (positive) set by combining benchmark and candidate pairs.

Steps:
1. Flatten scope40/pdbstyle-2.08 to scope_pdb with .pdb extension
2. Load benchmark and candidate CP pairs
3. Normalize pair order (alphabetically)
4. Combine and deduplicate
5. Create CP positive set with PDB files and metadata
"""

import os
import shutil
import pandas as pd
import logging
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Paths
BASE_DIR = Path("/mnt/scratch/jugipalace/foldseek_new_3di")
PDBSTYLE_DIR = BASE_DIR / "scope40" / "pdbstyle-2.08"
SCOPE_PDB_DIR = BASE_DIR / "scope_pdb"
BENCHMARK_FILE = BASE_DIR / "cirpin" / "scope40" / "benchmark" / "cp_benchmark_db.csv"
CANDIDATE_FILE = BASE_DIR / "cirpin" / "cp_pairs_scope40_dedup_asym.tsv"
CP_POSITIVE_DIR = BASE_DIR / "cp_positive"
CP_POSITIVE_PAIRS_FILE = CP_POSITIVE_DIR / "cp_pairs.tsv"
CP_POSITIVE_PDB_DIR = CP_POSITIVE_DIR / "pdb"

def flatten_pdbstyle_to_pdb():
    """Flatten pdbstyle-2.08 directory and convert .ent to .pdb in scope_pdb."""
    logging.info("=" * 80)
    logging.info("Step 1: Flatten pdbstyle-2.08 to scope_pdb with .pdb extension")
    logging.info("=" * 80)
    
    SCOPE_PDB_DIR.mkdir(parents=True, exist_ok=True)
    
    count = 0
    for subdirs in PDBSTYLE_DIR.rglob("*.ent"):
        try:
            # Extract domain ID (filename without extension)
            domain_id = subdirs.stem  # e.g., "d12asa_" from "d12asa_.ent"
            output_file = SCOPE_PDB_DIR / f"{domain_id}.pdb"
            
            # Skip if already exists
            if output_file.exists():
                continue
            
            # Copy .ent file as .pdb
            shutil.copy2(subdirs, output_file)
            count += 1
            
            if count % 100 == 0:
                logging.info(f"  Copied {count} files...")
        
        except Exception as e:
            logging.error(f"Error processing {subdirs}: {e}")
    
    logging.info(f"✓ Flattened {count} files from pdbstyle-2.08 to scope_pdb")
    logging.info(f"  Total PDB files in scope_pdb: {len(list(SCOPE_PDB_DIR.glob('*.pdb')))}")

def normalize_pair_order(domA, domB):
    """Normalize pair order so (A,B) and (B,A) are treated the same."""
    # Convert to string to handle NaN values
    domA = str(domA) if pd.notna(domA) else ""
    domB = str(domB) if pd.notna(domB) else ""
    return tuple(sorted([domA, domB]))

def load_and_combine_pairs():
    """Load benchmark and candidate files, combine and deduplicate."""
    logging.info("\n" + "=" * 80)
    logging.info("Step 2: Load and combine benchmark and candidate pairs")
    logging.info("=" * 80)
    
    # Load benchmark pairs
    logging.info(f"Loading benchmark from: {BENCHMARK_FILE}")
    benchmark_df = pd.read_csv(BENCHMARK_FILE)
    logging.info(f"  Benchmark pairs: {len(benchmark_df)}")
    logging.info(f"  Columns: {benchmark_df.columns.tolist()}")
    
    # Extract relevant columns (Query and Target)
    # Based on the CSV structure: Query, ..., Target columns
    benchmark_pairs = benchmark_df[['Query', 'Target']].copy()
    benchmark_pairs.columns = ['domA', 'domB']
    benchmark_pairs['source'] = 'benchmark'
    
    logging.info(f"  Extracted pairs: {len(benchmark_pairs)}")
    
    # Load candidate pairs
    logging.info(f"\nLoading candidates from: {CANDIDATE_FILE}")
    candidate_df = pd.read_csv(CANDIDATE_FILE, sep='\t')
    logging.info(f"  Candidate pairs: {len(candidate_df)}")
    logging.info(f"  Columns: {candidate_df.columns.tolist()}")
    
    # Extract relevant columns (query, target)
    candidate_pairs = candidate_df[['query', 'target']].copy()
    candidate_pairs.columns = ['domA', 'domB']
    candidate_pairs['source'] = 'candidate'
    
    logging.info(f"  Extracted pairs: {len(candidate_pairs)}")
    
    # Combine
    combined_df = pd.concat([benchmark_pairs, candidate_pairs], ignore_index=True)
    logging.info(f"\n✓ Combined pairs: {len(combined_df)}")
    
    # Remove rows with NaN values
    initial_len = len(combined_df)
    combined_df = combined_df.dropna(subset=['domA', 'domB'])
    removed_nan = initial_len - len(combined_df)
    if removed_nan > 0:
        logging.info(f"Removed {removed_nan} pairs with missing domain IDs")
    
    # Normalize pair order
    logging.info("Normalizing pair order (alphabetically)...")
    combined_df[['domA', 'domB']] = combined_df.apply(
        lambda row: pd.Series(normalize_pair_order(row['domA'], row['domB'])),
        axis=1
    )
    
    # Drop duplicates (keeping first occurrence)
    logging.info("Dropping duplicates...")
    before_dedup = len(combined_df)
    combined_df = combined_df.drop_duplicates(subset=['domA', 'domB'], keep='first')
    after_dedup = len(combined_df)
    logging.info(f"  Duplicates removed: {before_dedup - after_dedup}")
    logging.info(f"  Final unique pairs: {after_dedup}")
    
    return combined_df

def build_cp_positive_set(pairs_df):
    """Build CP positive set with PDB files and metadata."""
    logging.info("\n" + "=" * 80)
    logging.info("Step 3: Build CP positive set")
    logging.info("=" * 80)
    
    # Create output directories
    CP_POSITIVE_DIR.mkdir(parents=True, exist_ok=True)
    CP_POSITIVE_PDB_DIR.mkdir(parents=True, exist_ok=True)
    
    logging.info(f"Output directory: {CP_POSITIVE_DIR}")
    logging.info(f"PDB directory: {CP_POSITIVE_PDB_DIR}")
    
    # Copy PDB files
    logging.info("\nCopying PDB files...")
    all_domains = set()
    all_domains.update(pairs_df['domA'].values)
    all_domains.update(pairs_df['domB'].values)
    
    missing_pdbs = []
    copied_count = 0
    
    for domain in sorted(all_domains):
        # Try with and without trailing underscore
        domain_id = domain if domain.endswith('_') else domain + '_'
        pdb_file = SCOPE_PDB_DIR / f"{domain_id}.pdb"
        
        if not pdb_file.exists():
            # Try without underscore
            domain_id = domain.rstrip('_')
            pdb_file = SCOPE_PDB_DIR / f"{domain_id}.pdb"
        
        if pdb_file.exists():
            output_pdb = CP_POSITIVE_PDB_DIR / pdb_file.name
            shutil.copy2(pdb_file, output_pdb)
            copied_count += 1
        else:
            missing_pdbs.append(domain)
    
    logging.info(f"✓ Copied {copied_count} PDB files")
    if missing_pdbs:
        logging.warning(f"  Missing PDB files: {len(missing_pdbs)}")
        for domain in missing_pdbs[:10]:  # Show first 10
            logging.warning(f"    - {domain}")
        if len(missing_pdbs) > 10:
            logging.warning(f"    ... and {len(missing_pdbs) - 10} more")
    
    # Save metadata
    logging.info("\nSaving metadata...")
    pairs_df.to_csv(CP_POSITIVE_PAIRS_FILE, sep='\t', index=False)
    logging.info(f"✓ Saved pairs metadata: {CP_POSITIVE_PAIRS_FILE}")
    logging.info(f"  Total pairs: {len(pairs_df)}")
    
    # Print summary
    logging.info("\n" + "=" * 80)
    logging.info("CP Positive Set Summary")
    logging.info("=" * 80)
    logging.info(f"Total pairs: {len(pairs_df)}")
    logging.info(f"Total unique domains: {len(all_domains)}")
    logging.info(f"PDB files copied: {copied_count}")
    logging.info(f"PDB files missing: {len(missing_pdbs)}")
    logging.info(f"Coverage: {100 * copied_count / (copied_count + len(missing_pdbs)):.1f}%")
    
    # Print source distribution
    logging.info("\nSource distribution:")
    for source, count in pairs_df['source'].value_counts().items():
        logging.info(f"  {source}: {count}")

def main():
    logging.info("Starting CP Positive Set Builder")
    logging.info("=" * 80)
    
    # Step 1: Flatten pdbstyle-2.08
    flatten_pdbstyle_to_pdb()
    
    # Step 2: Load and combine pairs
    pairs_df = load_and_combine_pairs()
    
    # Step 3: Build CP positive set
    build_cp_positive_set(pairs_df)
    
    logging.info("\n" + "=" * 80)
    logging.info("✓ CP Positive Set Builder completed successfully!")
    logging.info("=" * 80)

if __name__ == "__main__":
    main()
