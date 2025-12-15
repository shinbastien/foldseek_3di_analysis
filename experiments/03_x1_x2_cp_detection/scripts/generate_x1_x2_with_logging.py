#!/usr/bin/env python3
"""
Generate X1 and X2 3Di sequences for each alphabet (8f, 10f) with detailed logging.

Outputs:
- x1_x2_compare_noncp_random/8f/3di_sequences/X1_sequences.fasta
- x1_x2_compare_noncp_random/8f/3di_sequences/X2_sequences.fasta
- x1_x2_compare_noncp_random/10f/3di_sequences/X1_sequences.fasta
- x1_x2_compare_noncp_random/10f/3di_sequences/X2_sequences.fasta
- x1_x2_compare_noncp_random/processing.log
"""

import pandas as pd
from pathlib import Path
import logging
import sys

# Setup logging
log_file = Path('x1_x2_compare_noncp_random/processing.log')
log_file.parent.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def load_fasta_dict(fasta_path):
    """Load FASTA file into dictionary."""
    seqs = {}
    with open(fasta_path) as f:
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    seqs[current_header] = "".join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            seqs[current_header] = "".join(current_seq)
    return seqs

def generate_x1_x2(alphabet, fasta_path, pairs_df, out_dir):
    """Generate X1 and X2 FASTA files for given alphabet."""
    logger.info(f"\n{'='*80}")
    logger.info(f"Processing {alphabet.upper()}")
    logger.info(f"{'='*80}")
    
    # Load sequences
    logger.info(f"Loading {alphabet.upper()} sequences from {fasta_path}")
    seqs = load_fasta_dict(fasta_path)
    logger.info(f"  ✓ Loaded {len(seqs)} sequences")
    
    # Prepare directories
    out_dir_path = Path(out_dir)
    out_dir_path.mkdir(parents=True, exist_ok=True)
    
    # Generate X1: each domain as-is
    logger.info(f"\nGenerating X1 (single) sequences...")
    x1_file = out_dir_path / "X1_sequences.fasta"
    all_domains = sorted(pd.concat([pairs_df['query_id'], pairs_df['target_id']]).unique())
    
    x1_count = 0
    with open(x1_file, 'w') as f:
        for domain_id in all_domains:
            seq = seqs.get(domain_id)
            if seq:
                f.write(f">{domain_id}\n{seq}\n")
                x1_count += 1
            else:
                logger.warning(f"  ⚠ Missing sequence for {domain_id}")
    
    logger.info(f"  ✓ Wrote {x1_count} X1 sequences to {x1_file}")
    
    # Generate X2: query domains doubled
    logger.info(f"\nGenerating X2 (query doubled) sequences...")
    x2_file = out_dir_path / "X2_sequences.fasta"
    
    x2_count = 0
    with open(x2_file, 'w') as f:
        for query_id in pairs_df['query_id'].unique():
            seq = seqs.get(query_id)
            if seq:
                doubled_seq = seq + seq
                f.write(f">{query_id}_x2\n{doubled_seq}\n")
                x2_count += 1
            else:
                logger.warning(f"  ⚠ Missing sequence for {query_id}")
    
    logger.info(f"  ✓ Wrote {x2_count} X2 sequences to {x2_file}")
    
    return x1_file, x2_file

def main():
    logger.info("="*80)
    logger.info("X1/X2 3Di Sequence Generation with Logging")
    logger.info("="*80)
    
    # Load pairs
    pairs_file = 'noncp_work/noncp_homolog_pairs.tsv'
    logger.info(f"\nLoading pairs from {pairs_file}")
    pairs_df = pd.read_csv(pairs_file, sep='\t')
    logger.info(f"  ✓ Loaded {len(pairs_df)} pairs")
    logger.info(f"  Unique queries: {len(pairs_df['query_id'].unique())}")
    logger.info(f"  Unique targets: {len(pairs_df['target_id'].unique())}")
    logger.info(f"  Unique domains: {len(pd.concat([pairs_df['query_id'], pairs_df['target_id']]).unique())}")
    
    # Process 8f
    logger.info(f"\n{'-'*80}")
    logger.info("8F PROCESSING")
    logger.info(f"{'-'*80}")
    
    try:
        x1_8f, x2_8f = generate_x1_x2(
            '8f',
            'x1_x2_compare_noncp_random/8f/tmp/sequences_8f.fasta',
            pairs_df,
            'x1_x2_compare_noncp_random/8f/3di_sequences'
        )
        logger.info(f"  X1 file: {x1_8f}")
        logger.info(f"  X2 file: {x2_8f}")
    except Exception as e:
        logger.error(f"ERROR in 8f processing: {e}", exc_info=True)
        return False
    
    # Process 10f
    logger.info(f"\n{'-'*80}")
    logger.info("10F PROCESSING")
    logger.info(f"{'-'*80}")
    
    try:
        x1_10f, x2_10f = generate_x1_x2(
            '10f',
            'x1_x2_compare_noncp_random/10f/tmp/sequences_10f.fasta',
            pairs_df,
            'x1_x2_compare_noncp_random/10f/3di_sequences'
        )
        logger.info(f"  X1 file: {x1_10f}")
        logger.info(f"  X2 file: {x2_10f}")
    except Exception as e:
        logger.error(f"ERROR in 10f processing: {e}", exc_info=True)
        return False
    
    # Summary
    logger.info(f"\n{'='*80}")
    logger.info("SUMMARY")
    logger.info(f"{'='*80}")
    logger.info(f"\n8F Outputs:")
    logger.info(f"  X1: x1_x2_compare_noncp_random/8f/3di_sequences/X1_sequences.fasta")
    logger.info(f"  X2: x1_x2_compare_noncp_random/8f/3di_sequences/X2_sequences.fasta")
    logger.info(f"\n10F Outputs:")
    logger.info(f"  X1: x1_x2_compare_noncp_random/10f/3di_sequences/X1_sequences.fasta")
    logger.info(f"  X2: x1_x2_compare_noncp_random/10f/3di_sequences/X2_sequences.fasta")
    logger.info(f"\nLog file: {log_file}")
    logger.info(f"{'='*80}\n")
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
