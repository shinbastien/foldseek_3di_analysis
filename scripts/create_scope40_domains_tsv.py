#!/usr/bin/env python3
"""
Create scope40_domains.tsv from existing sources.
Combines SCCS dictionary with PDB length information.
"""

import pickle
import pandas as pd
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')

# Paths
BASE_DIR = Path("/mnt/scratch/jugipalace/foldseek_new_3di")
SCCS_PICKLE = BASE_DIR / "cirpin" / "scope40_sid_sccs_dictionary" / "dir.des.scope.2.08-2023-01-06_sid_sccs_dictionary.pkl"
PDB_DIR = BASE_DIR / "scope_pdb"
OUTPUT_FILE = BASE_DIR / "scope40_domains.tsv"

def get_pdb_length(pdb_file):
    """Get sequence length from PDB file by counting unique residues."""
    try:
        with open(pdb_file, 'r') as f:
            atoms = 0
            prev_resnum = None
            for line in f:
                if line.startswith('ATOM'):
                    resnum = int(line[22:26].strip())
                    if resnum != prev_resnum:
                        atoms += 1
                        prev_resnum = resnum
            return atoms if atoms > 0 else None
    except:
        return None

def extract_superfamily(sccs):
    """Extract superfamily (a.b.c) from SCCS (a.b.c.d)."""
    if pd.isna(sccs) or sccs == '-':
        return None
    parts = str(sccs).split('.')
    if len(parts) >= 3:
        return f"{parts[0]}.{parts[1]}.{parts[2]}"
    return None

def extract_species_from_domain(domain_id):
    """
    Extract species ID from domain ID.
    SCOPe domain format: d<PDB><chain><species>
    Example: d1jx6a_ means PDB=1jx6, chain=a, species_id from first digit
    For simplicity, use PDB code as proxy for species grouping.
    """
    # Extract PDB code (first 4-5 characters after 'd')
    if domain_id.startswith('d') and len(domain_id) >= 5:
        # Try to extract numeric/alpha part
        pdb_code = domain_id[1:5]  # 4-char PDB code
        return pdb_code.lower()
    return domain_id

def main():
    logging.info("="*80)
    logging.info("Creating scope40_domains.tsv")
    logging.info("="*80)
    
    # Load SCCS dictionary
    logging.info(f"\nLoading SCCS from: {SCCS_PICKLE}")
    with open(SCCS_PICKLE, 'rb') as f:
        sccs_dict = pickle.load(f)
    logging.info(f"  Loaded {len(sccs_dict)} domains")
    
    # Build domain info
    logging.info("\nBuilding domain information...")
    domains_list = []
    missing_count = 0
    
    for domain_id, sccs in sccs_dict.items():
        pdb_file = PDB_DIR / f"{domain_id}.pdb"
        
        # Only include if PDB exists
        if not pdb_file.exists():
            missing_count += 1
            continue
        
        # Get length
        length = get_pdb_length(pdb_file)
        if length is None:
            continue
        
        # Get SF
        sf = extract_superfamily(sccs)
        if sf is None:
            continue
        
        # Get species (PDB-based proxy)
        species = extract_species_from_domain(domain_id)
        
        domains_list.append({
            'domain_id': domain_id,
            'length': length,
            'sccs': sccs,
            'superfamily': sf,
            'species_id': species
        })
    
    domains_df = pd.DataFrame(domains_list)
    logging.info(f"  Found {len(domains_df)} domains with PDB and SCCS")
    logging.info(f"  Missing PDB files: {missing_count}")
    
    # Save TSV
    logging.info(f"\nSaving to: {OUTPUT_FILE}")
    domains_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
    logging.info(f"  ✓ Saved {len(domains_df)} domains")
    
    # Summary
    logging.info("\nSummary:")
    logging.info(f"  Total domains: {len(domains_df)}")
    logging.info(f"  Unique SFs: {domains_df['superfamily'].nunique()}")
    logging.info(f"  Length range: {domains_df['length'].min()}-{domains_df['length'].max()}")
    logging.info(f"  Mean length: {domains_df['length'].mean():.1f}")

if __name__ == "__main__":
    main()
