#!/usr/bin/env python3
"""
compare_8f_10f_ssw_v2.py

Enhanced version that:
- Uses foldseek createdb for 10f (official standard)
- Uses pdb_to_3di.py for 8f (custom model)
- Saves sequences in organized directory structure
- Saves all SSW results to files
- Logs complete SSW outputs
"""

import os
import pandas as pd
import logging
import sys
import re
import datetime
import argparse
import subprocess
import shutil
from pathlib import Path

# Ensure the directory containing pdb_to_3di.py is in the Python path
sys.path.insert(0, '/mnt/scratch/jugipalace/foldseek_new_3di')

from pdb_to_3di import pdb_to_3di

# Paths and constants
INPUT_FILE = "/mnt/scratch/jugipalace/foldseek_new_3di/cirpin/scopealigngood_not_cp_pairs_scope40.txt"
PDB_DIR = "/mnt/scratch/jugipalace/foldseek_new_3di/scope_pdb/"
OUTPUT_DIR = "/mnt/scratch/jugipalace/foldseek_new_3di/x1_x2_compare"
MODEL_DIR = "/mnt/scratch/jugipalace/foldseek_new_3di/encoders_and_tools/training_3di_gpu_"  # Base directory for models
SSW_BINARY = "/mnt/scratch/jugipalace/foldseek_new_3di/ssw/tmp/ssw/src/ssw_test"
FOLDSEEK_BIN = "foldseek"

# Argument parsing
parser = argparse.ArgumentParser(description="Compare 8f and 10f encodings using SSW scores.")
parser.add_argument("--output_dir", type=str, default=OUTPUT_DIR, help="Custom output directory for results.")
parser.add_argument("--pair_file", type=str, required=True, help="Path to the file containing query-target pairs.")
args = parser.parse_args()

# Use the provided output directory
OUTPUT_DIR = args.output_dir

def setup_output_directories(output_dir):
    """Create organized output directory structure."""
    dirs = [
        os.path.join(output_dir, "logs"),
        os.path.join(output_dir, "ssw_results"),
        os.path.join(output_dir, "3di_sequences", "8f"),
        os.path.join(output_dir, "3di_sequences", "10f")
    ]
    for d in dirs:
        os.makedirs(d, exist_ok=True)
    return dirs

# Setup directories
setup_output_directories(OUTPUT_DIR)

# Configure logging
log_dir = os.path.join(OUTPUT_DIR, "logs")
log_file = os.path.join(log_dir, f"batch_run_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log")

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s', 
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file)
    ]
)

# Utility functions from tmalign_3di_match_pipeline.py
def cleanse_createdb_file(path: Path) -> str:
    """Clean createdb output file."""
    raw = path.read_bytes()
    try:
        text = raw.decode('utf-8')
    except UnicodeDecodeError:
        text = raw.decode('latin-1')
    lines = []
    for ln in text.splitlines():
        s = ln.rstrip('\x00').rstrip()
        s = re.sub(r"(\t|\s)+(0|NULL)$", '', s)
        lines.append(s)
    return '\n'.join(lines) + '\n'

def extract_sequence_from_file_content(content: str) -> str:
    """Extract sequence from file content."""
    lines = [l.strip() for l in content.splitlines() if l.strip()]
    if not lines:
        return ''
    if lines[0].startswith('>'):
        return ''.join(lines[1:])
    return ''.join(lines)

def write_3di_fasta(out_path: Path, header: str, token_seq: str):
    """Write 3Di sequence in FASTA format."""
    out_path.write_text(f">{header}\n{token_seq}\n")

def run_foldseek_createdb(pdb_file: str, output_dir: str) -> tuple:
    """
    Run foldseek createdb and extract 10f 3Di sequence.
    Returns: (token_seq, token_seq_x2, fasta_path, fasta_path_x2)
    """
    pdb_path = Path(pdb_file).resolve()  # Use absolute path
    base_name = pdb_path.stem
    
    # Create temporary directory for createdb (use absolute path)
    temp_dir = Path(output_dir).resolve() / "temp" / base_name
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy PDB to temp directory
    pdb_tmp = temp_dir / pdb_path.name
    shutil.copy2(str(pdb_path), str(pdb_tmp))
    
    # Run foldseek createdb
    prefix = f"{base_name}_pdzdb"
    cmd = [FOLDSEEK_BIN, 'createdb', str(pdb_tmp), prefix]
    
    # Log the command
    logging.info(f"Running Foldseek createdb for {base_name}:")
    logging.info(f"  Command: {' '.join(cmd)}")
    logging.info(f"  Working directory: {temp_dir}")
    
    try:
        result = subprocess.run(cmd, cwd=str(temp_dir), capture_output=True, text=True, check=False)
        
        # Log stdout/stderr
        if result.stdout:
            logging.debug(f"  Foldseek stdout: {result.stdout}")
        if result.stderr:
            logging.debug(f"  Foldseek stderr: {result.stderr}")
        
        if result.returncode != 0:
            logging.warning(f"Foldseek createdb failed for {pdb_file} (exit code {result.returncode})")
            logging.warning(f"  Error: {result.stderr}")
            return None, None, None, None
        
        logging.info(f"  Foldseek createdb completed successfully")
    except Exception as e:
        logging.error(f"Failed to run foldseek createdb: {e}")
        return None, None, None, None
    
    # Find and parse _ss file
    ss_candidates = [
        temp_dir / f"{prefix}_ss",
        temp_dir / f"{prefix}_pdzdb_ss"
    ]
    
    pdz_ss = None
    for candidate in ss_candidates:
        if candidate.exists():
            pdz_ss = candidate
            break
    
    if pdz_ss is None:
        # Search in directory
        for f in temp_dir.iterdir():
            if f.name in (f"{prefix}_ss", f"{prefix}_pdzdb_ss"):
                pdz_ss = f
                break
    
    if pdz_ss is None or not pdz_ss.exists():
        logging.warning(f"Could not find _ss file for {pdb_file}")
        return None, None, None, None
    
    # Extract token sequence
    text = cleanse_createdb_file(pdz_ss)
    token_seq = extract_sequence_from_file_content(text)
    
    if not token_seq:
        logging.warning(f"Empty token sequence for {pdb_file}")
        return None, None, None, None
    
    # Save FASTA files
    encoding_dir = Path(output_dir) / "3di_sequences" / "10f"
    encoding_dir.mkdir(parents=True, exist_ok=True)
    
    fasta_path = encoding_dir / f"{base_name}_10f_3di.fasta"
    fasta_path_x2 = encoding_dir / f"{base_name}_10f_3di_x2.fasta"
    
    write_3di_fasta(fasta_path, f"{base_name}_10f_3di", token_seq)
    write_3di_fasta(fasta_path_x2, f"{base_name}_10f_3di_x2", token_seq + token_seq)
    
    logging.info(f"Created 10f 3Di: {fasta_path} and {fasta_path_x2}")
    
    # Cleanup temp directory
    try:
        shutil.rmtree(temp_dir)
    except Exception:
        pass
    
    return token_seq, token_seq + token_seq, str(fasta_path), str(fasta_path_x2)

def pdb_to_3di_with_output(pdb_file: str, model_dir: str, output_dir: str, encoding: str) -> tuple:
    """
    Generate 3Di sequences for 8f or 10f.
    Returns: (seq, seq_x2, output_path, output_path_x2)
    """
    base_name = os.path.basename(pdb_file).replace(".pdb", "")
    
    if encoding == "10f":
        # Use foldseek createdb
        return run_foldseek_createdb(pdb_file, output_dir)
    else:
        # Use pdb_to_3di.py for 8f
        encoding_dir = os.path.join(output_dir, "3di_sequences", encoding)
        os.makedirs(encoding_dir, exist_ok=True)
        
        output_path = os.path.join(encoding_dir, f"{base_name}_{encoding}_3di.fasta")
        output_path_x2 = os.path.join(encoding_dir, f"{base_name}_{encoding}_3di_x2.fasta")
        
        # Generate 3Di sequence
        seq, _ = pdb_to_3di(pdb_file, model_dir, output_path=output_path)
        
        if not seq:
            return None, None, None, None
        
        # Create x2 version
        seq_x2 = seq + seq
        write_3di_fasta(Path(output_path_x2), f"{base_name}_{encoding}_3di_x2", seq_x2)
        
        logging.info(f"Created {encoding} 3Di: {output_path} and {output_path_x2}")
        
        return seq, seq_x2, output_path, output_path_x2

def compute_ssw_score_with_binary(ssw_binary: str, query_seq: str, target_seq: str,
                                   query_name: str, target_name: str, encoding: str,
                                   is_x2: bool, output_dir: str, gap_open: int, gap_extend: int) -> int:
    """
    Compute SSW alignment score and save detailed results.
    """
    # Save sequences to temporary files
    query_file = os.path.join(output_dir, "query_tmp.fasta")
    target_file = os.path.join(output_dir, "target_tmp.fasta")
    
    with open(query_file, "w") as qf:
        qf.write(f">{query_name}_{encoding}_3di\n{query_seq}")
    
    with open(target_file, "w") as tf:
        tf.write(f">{target_name}_{encoding}_3di\n{target_seq}")
    
    # Run SSW (with -c option to get alignment path)
    cmd = [ssw_binary, '-o', str(gap_open), '-e', str(gap_extend), '-p', '-c', query_file, target_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"SSW binary failed with error: {result.stderr}")
    
    # Parse and display alignment
    x_suffix = "X2" if is_x2 else "X1"
    logging.info(f"\n{'='*80}")
    logging.info(f"SSW {x_suffix} Result: {query_name} vs {target_name} ({encoding})")
    logging.info(f"{'='*80}")
    
    # Parse the SSW output to extract and highlight the alignment
    lines = result.stdout.split('\n')
    in_alignment = False
    alignment_lines = []
    
    for line in lines:
        # Log all output
        logging.info(line)
        
        # Collect alignment section
        if 'target_begin' in line or in_alignment:
            in_alignment = True
            alignment_lines.append(line)
            if line.strip() == '' and len(alignment_lines) > 5:
                in_alignment = False
    
    # Highlight alignment section if found
    if alignment_lines:
        logging.info(f"\n{'-'*80}")
        logging.info("ALIGNMENT VISUALIZATION:")
        logging.info(f"{'-'*80}")
        for line in alignment_lines:
            logging.info(line)
        logging.info(f"{'-'*80}")
    
    logging.info(f"{'='*80}\n")
    
    # Save result to file
    ssw_results_dir = os.path.join(output_dir, "ssw_results")
    x_label = "x2" if is_x2 else "x1"
    output_filename = f"ssw_{query_name}_{x_label}_vs_{target_name}_{encoding}.txt"
    output_path = os.path.join(ssw_results_dir, output_filename)
    
    with open(output_path, "w") as f:
        f.write(result.stdout)
    
    logging.info(f"SSW result saved to: {output_path}")
    
    # Clean up temporary files
    try:
        if os.path.exists(query_file):
            os.remove(query_file)
        if os.path.exists(target_file):
            os.remove(target_file)
    except Exception as e:
        logging.debug(f"Could not remove temporary files: {e}")
    
    # Extract score
    match = re.search(r"optimal_alignment_score:\s*(\d+)", result.stdout)
    if match:
        score = int(match.group(1))
        logging.info(f"Optimal alignment score: {score}")
        return score
    
    raise ValueError("Could not find alignment score in SSW output")

def process_pair(query: str, target: str, pdb_dir: str, model_dir: str, ssw_binary: str) -> list:
    """Process a single pair to compute SSW scores for 8f and 10f encodings."""
    results = []
    
    # Check if both PDB files exist
    query_pdb = os.path.join(pdb_dir, f"{query.lower()}.pdb")
    target_pdb = os.path.join(pdb_dir, f"{target.lower()}.pdb")
    
    if not os.path.exists(query_pdb):
        logging.warning(f"SKIPPED: Query PDB file not found: {query_pdb}")
        return results
    
    if not os.path.exists(target_pdb):
        logging.warning(f"SKIPPED: Target PDB file not found: {target_pdb}")
        return results
    
    for encoding in ["8f", "10f"]:
        try:
            logging.info(f"\n{'#'*60}")
            logging.info(f"Processing pair: {query} vs {target}, encoding: {encoding}")
            logging.info(f"{'#'*60}")
            
            # Construct model directory for the specific encoding
            if encoding == "10f":
                encoding_model_dir = None  # Will use foldseek
            else:
                encoding_model_dir = os.path.join(f"{model_dir}{encoding}", "tmp")
            
            # Generate 3Di sequences
            query_seq, query_seq_x2, query_output, query_output_x2 = pdb_to_3di_with_output(
                query_pdb, encoding_model_dir, OUTPUT_DIR, encoding
            )
            target_seq, target_seq_x2, target_output, target_output_x2 = pdb_to_3di_with_output(
                target_pdb, encoding_model_dir, OUTPUT_DIR, encoding
            )
            
            if not query_seq or not target_seq:
                logging.error(f"Failed to generate 3Di sequences for {query} vs {target} ({encoding})")
                continue
            
            # Log 3Di sequence information
            logging.info(f"Query: {query}, PDB: {query_pdb}")
            logging.info(f"  Length: {len(query_seq)}")
            logging.info(f"  3Di Sequence: {query_seq}")
            logging.info(f"  Output: {query_output}")
            logging.info(f"  Output X2: {query_output_x2}")
            
            logging.info(f"Target: {target}, PDB: {target_pdb}")
            logging.info(f"  Length: {len(target_seq)}")
            logging.info(f"  3Di Sequence: {target_seq}")
            logging.info(f"  Output: {target_output}")
            logging.info(f"  Output X2: {target_output_x2}")
            
            # Compute SSW X1 score (single sequence)
            logging.info(f"\nComputing SSW X1 score...")
            score_x1 = compute_ssw_score_with_binary(
                ssw_binary, query_seq, target_seq,
                query, target, encoding, False, OUTPUT_DIR, 8, 2
            )
            
            # Compute SSW X2 score (doubled query sequence)
            logging.info(f"\nComputing SSW X2 score...")
            score_x2 = compute_ssw_score_with_binary(
                ssw_binary, query_seq_x2, target_seq,
                query, target, encoding, True, OUTPUT_DIR, 8, 2
            )
            
            # Calculate difference
            score_diff = score_x2 - score_x1
            
            logging.info(f"\nSummary for {query} vs {target} ({encoding}):")
            logging.info(f"  Score X1: {score_x1}")
            logging.info(f"  Score X2: {score_x2}")
            logging.info(f"  Difference (X2 - X1): {score_diff}")
            
            results.append({
                "query": query,
                "target": target,
                "encoding": encoding,
                "query_output": query_output,
                "target_output": target_output,
                "query_output_x2": query_output_x2,
                "target_output_x2": target_output_x2,
                "score_x1": score_x1,
                "score_x2": score_x2,
                "score_diff": score_diff
            })
            
        except Exception as e:
            logging.error(f"Error processing pair {query} vs {target} with encoding {encoding}: {str(e)}", exc_info=True)
    
    return results

def read_pairs_from_file(pair_file: str) -> pd.DataFrame:
    """Read query-target pairs from a file, skipping the first row."""
    df = pd.read_csv(pair_file, sep="\t", header=None, usecols=[0, 1], 
                     names=["query", "target"], dtype=str, skiprows=1)
    df = df.drop_duplicates(subset=["query", "target"])
    return df

def main():
    logging.info(f"Starting comparison pipeline")
    logging.info(f"Output directory: {OUTPUT_DIR}")
    logging.info(f"Pair file: {args.pair_file}")
    logging.info(f"PDB directory: {PDB_DIR}")
    logging.info(f"SSW binary: {SSW_BINARY}")
    
    # Read and deduplicate input pairs
    pairs = read_pairs_from_file(args.pair_file)
    logging.info(f"Loaded {len(pairs)} unique pairs")
    
    all_results = []
    for idx, row in pairs.iterrows():
        query, target = row["query"], row["target"]
        logging.info(f"\n{'='*80}")
        logging.info(f"Processing pair {idx+1}/{len(pairs)}: {query} vs {target}")
        logging.info(f"{'='*80}")
        
        results = process_pair(query, target, PDB_DIR, MODEL_DIR, SSW_BINARY)
        all_results.extend(results)
    
    # Save results to CSV
    if all_results:
        results_df = pd.DataFrame(all_results)
        output_csv = os.path.join(OUTPUT_DIR, "comparison_results.csv")
        results_df.to_csv(output_csv, index=False)
        logging.info(f"\nSaved results to: {output_csv}")
        logging.info(f"Total results: {len(results_df)}")
    else:
        logging.warning("No results to save")
    
    logging.info(f"\nPipeline completed successfully!")

if __name__ == "__main__":
    main()
