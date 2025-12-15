#!/usr/bin/env python3
"""
Analyze CP positive set for SCCS and length distribution.
This will guide stratified sampling for non-CP query set.
"""

import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

# Setup paths
BASE_DIR = Path("/mnt/scratch/jugipalace/foldseek_new_3di")
CP_PAIRS_FILE = BASE_DIR / "cp_positive" / "cp_pairs.tsv"
SCCS_PICKLE = BASE_DIR / "cirpin" / "scope40_sid_sccs_dictionary" / "dir.des.scope.2.08-2023-01-06_sid_sccs_dictionary.pkl"
PDB_DIR = BASE_DIR / "scope_pdb"
OUTPUT_DIR = BASE_DIR / "analysis"

OUTPUT_DIR.mkdir(exist_ok=True)

def load_sccs_dictionary():
    """Load SCCS dictionary from pickle."""
    print("Loading SCCS dictionary...")
    with open(SCCS_PICKLE, 'rb') as f:
        sccs_dict = pickle.load(f)
    print(f"  Loaded {len(sccs_dict)} domains")
    return sccs_dict

def extract_superfamily(sccs):
    """Extract superfamily (a.b.c) from SCCS (a.b.c.d)."""
    if pd.isna(sccs) or sccs == '-':
        return None
    parts = sccs.split('.')
    if len(parts) >= 3:
        return f"{parts[0]}.{parts[1]}.{parts[2]}"
    return None

def extract_fold(sccs):
    """Extract fold (a.b) from SCCS (a.b.c.d)."""
    if pd.isna(sccs) or sccs == '-':
        return None
    parts = sccs.split('.')
    if len(parts) >= 2:
        return f"{parts[0]}.{parts[1]}"
    return None

def get_pdb_length(pdb_file):
    """Get sequence length from PDB file by counting ATOM records."""
    try:
        with open(pdb_file, 'r') as f:
            atoms = 0
            prev_resnum = None
            for line in f:
                if line.startswith('ATOM'):
                    # Extract residue number (columns 22-26)
                    resnum = int(line[22:26].strip())
                    # Count unique residues (CA atoms would be better but count all ATOM)
                    if resnum != prev_resnum:
                        atoms += 1
                        prev_resnum = resnum
            return atoms if atoms > 0 else None
    except:
        return None

def analyze_cp_set(sccs_dict):
    """Analyze CP positive set distribution."""
    print("\n" + "="*80)
    print("Analyzing CP Positive Set")
    print("="*80)
    
    # Load CP pairs
    print(f"\nLoading CP pairs from: {CP_PAIRS_FILE}")
    cp_df = pd.read_csv(CP_PAIRS_FILE, sep='\t')
    print(f"  Loaded {len(cp_df)} pairs")
    
    # Collect all unique domains
    domains = set()
    domains.update(cp_df['domA'].values)
    domains.update(cp_df['domB'].values)
    domains = sorted(list(domains))
    print(f"  Total unique domains: {len(domains)}")
    
    # Extract SCCS and lengths
    print("\nExtracting SCCS and lengths...")
    domain_info = []
    missing_sccs = 0
    missing_pdb = 0
    
    for domain in domains:
        # Get SCCS
        sccs = sccs_dict.get(domain, None)
        if sccs is None:
            missing_sccs += 1
            continue
        
        # Get length
        pdb_file = PDB_DIR / f"{domain}.pdb"
        if not pdb_file.exists():
            missing_pdb += 1
            continue
        
        length = get_pdb_length(pdb_file)
        if length is None:
            continue
        
        sf = extract_superfamily(sccs)
        fold = extract_fold(sccs)
        
        domain_info.append({
            'domain': domain,
            'sccs': sccs,
            'fold': fold,
            'superfamily': sf,
            'length': length
        })
    
    domain_df = pd.DataFrame(domain_info)
    print(f"  Found SCCS/length for {len(domain_df)} domains")
    print(f"  Missing SCCS: {missing_sccs}")
    print(f"  Missing PDB: {missing_pdb}")
    
    # Statistics
    print("\n" + "-"*80)
    print("SCCS Distribution (Superfamily)")
    print("-"*80)
    sf_counts = domain_df['superfamily'].value_counts()
    print(f"Total unique SFs: {len(sf_counts)}")
    print("\nTop 20 SFs by count:")
    print(sf_counts.head(20))
    
    print("\n" + "-"*80)
    print("Fold Distribution")
    print("-"*80)
    fold_counts = domain_df['fold'].value_counts()
    print(f"Total unique Folds: {len(fold_counts)}")
    print("\nTop 20 Folds by count:")
    print(fold_counts.head(20))
    
    print("\n" + "-"*80)
    print("Length Distribution")
    print("-"*80)
    print(f"Mean length: {domain_df['length'].mean():.1f}")
    print(f"Median length: {domain_df['length'].median():.1f}")
    print(f"Std dev: {domain_df['length'].std():.1f}")
    print(f"Min length: {domain_df['length'].min()}")
    print(f"Max length: {domain_df['length'].max()}")
    print(f"\nPercentiles:")
    for p in [10, 25, 50, 75, 90]:
        val = domain_df['length'].quantile(p/100)
        print(f"  {p}%: {val:.0f}")
    
    # Length bins for stratification
    print("\n" + "-"*80)
    print("Length Stratification (for sampling)")
    print("-"*80)
    bins = [0, 50, 100, 150, 200, 300, 500, 10000]
    bin_labels = ['<50', '50-100', '100-150', '150-200', '200-300', '300-500', '>500']
    domain_df['length_bin'] = pd.cut(domain_df['length'], bins=bins, labels=bin_labels, right=False)
    
    length_bin_counts = domain_df['length_bin'].value_counts().sort_index()
    print("\nLength bin distribution:")
    for bin_label, count in length_bin_counts.items():
        pct = 100 * count / len(domain_df)
        print(f"  {bin_label}: {count} ({pct:.1f}%)")
    
    # SF + Length bins for stratification
    print("\n" + "-"*80)
    print("SF + Length Stratification (for stratified sampling)")
    print("-"*80)
    
    domain_df['strata'] = domain_df['superfamily'] + "_" + domain_df['length_bin'].astype(str)
    strata_counts = domain_df['strata'].value_counts()
    print(f"Total strata combinations: {len(strata_counts)}")
    print(f"Strata sizes (min={strata_counts.min()}, max={strata_counts.max()}, mean={strata_counts.mean():.1f})")
    
    # Show top 10 strata
    print("\nTop 10 strata by size:")
    for strata, count in strata_counts.head(10).items():
        pct = 100 * count / len(domain_df)
        print(f"  {strata}: {count} ({pct:.1f}%)")
    
    # Visualization
    print("\n" + "-"*80)
    print("Summary statistics saved")
    print("-"*80)
    
    # Save detailed report
    report_file = OUTPUT_DIR / "cp_set_analysis.tsv"
    domain_df.to_csv(report_file, sep='\t', index=False)
    print(f"  Saved: {report_file}")
    
    return domain_df

def main():
    print("="*80)
    print("CP Positive Set Analysis")
    print("="*80)
    
    # Load SCCS dictionary
    sccs_dict = load_sccs_dictionary()
    
    # Analyze CP set
    domain_df = analyze_cp_set(sccs_dict)
    
    print("\n" + "="*80)
    print("✓ Analysis completed!")
    print("="*80)
    print(f"\nFor stratified sampling:")
    print(f"  - Total domains in CP set: {len(domain_df)}")
    print(f"  - Unique superfamilies: {domain_df['superfamily'].nunique()}")
    print(f"  - Unique folds: {domain_df['fold'].nunique()}")
    print(f"  - Length range: {domain_df['length'].min()}-{domain_df['length'].max()}")
    
if __name__ == "__main__":
    main()
