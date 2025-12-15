import os
import pandas as pd
import logging
import sys
import re
import datetime
import argparse

# Ensure the directory containing pdb_to_3di.py is in the Python path
sys.path.insert(0, '/mnt/scratch/jugipalace/foldseek_new_3di')

from pdb_to_3di import pdb_to_3di  # Import the actual pdb_to_3di function

# Correct the import path for pdb_to_3di
sys.path.insert(0, '/mnt/scratch/jugipalace/foldseek_new_3di')

# Add custom SSW library path
sys.path.insert(0, '/mnt/scratch/jugipalace/foldseek_new_3di/ssw/tmp/src')

# Paths and constants
INPUT_FILE = "/mnt/scratch/jugipalace/foldseek_new_3di/cirpin/scopealigngood_not_cp_pairs_scope40.txt"
PDB_DIR = "/mnt/scratch/jugipalace/foldseek_new_3di/scope_pdb/"
OUTPUT_DIR = "/mnt/scratch/jugipalace/foldseek_new_3di/x1_x2_compare"
MODEL_DIR = "/mnt/scratch/jugipalace/foldseek_new_3di/encoders_and_tools/training_3di_gpu_"  # Base directory for models
SSW_BINARY = "/mnt/scratch/jugipalace/foldseek_new_3di/ssw/tmp/ssw/src/ssw_test"  # Path to the SSW binary

# Add argument parsing for custom output directory
parser = argparse.ArgumentParser(description="Compare 8f and 10f encodings using SSW scores.")
parser.add_argument("--output_dir", type=str, default=OUTPUT_DIR, help="Custom output directory for results.")
# Add argument parsing for input pair file
parser.add_argument("--pair_file", type=str, required=True, help="Path to the file containing query-target pairs.")
args = parser.parse_args()

# Use the provided output directory
OUTPUT_DIR = args.output_dir
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Create a log file in the logs subdirectory of OUTPUT_DIR
log_dir = os.path.join(OUTPUT_DIR, "logs")
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f"batch_run_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log")

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[
    logging.StreamHandler(sys.stdout),
    logging.FileHandler(log_file)
])

# Helper functions
def read_scopealigngood(file_path):
    """Read the scopealigngood file and deduplicate pairs."""
    df = pd.read_csv(file_path, sep="\t", header=None, names=["query", "target"])
    df = df.drop_duplicates(subset=["query", "target"])
    return df

def compute_ssw_score_with_binary(ssw_binary, query_seq, target_seq, gap_open, gap_extend):
    """Compute the SSW alignment score using the ssw_test binary."""
    import subprocess
    # Save sequences to temporary files
    query_file = os.path.join(OUTPUT_DIR, "query_tmp.fasta")
    target_file = os.path.join(OUTPUT_DIR, "target_tmp.fasta")

    with open(query_file, "w") as qf:
        qf.write(f">query\n{query_seq}")

    with open(target_file, "w") as tf:
        tf.write(f">target\n{target_seq}")

    cmd = [ssw_binary, '-o', str(gap_open), '-e', str(gap_extend), '-p', query_file, target_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"SSW binary failed with error: {result.stderr}")

    # Use regex to extract the optimal alignment score
    match = re.search(r"optimal_alignment_score:\s*(\d+)", result.stdout)
    if match:
        return int(match.group(1))

    raise ValueError("Could not find alignment score in SSW output")

def pdb_to_3di_with_output(pdb_file, model_dir, output_dir):
    """Wrapper for pdb_to_3di to save outputs in the specified directory."""
    base_name = os.path.basename(pdb_file).replace(".pdb", "")
    output_path = os.path.join(output_dir, f"{base_name}_3di.txt")
    seq, _ = pdb_to_3di(pdb_file, model_dir, output_path=output_path)
    return seq, None, output_path

def process_pair(query, target, pdb_dir, model_dir, ssw_binary):
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
            # Correctly construct the model directory for the specific encoding
            encoding_model_dir = os.path.join(f"{model_dir}{encoding}", "tmp")

            query_seq, _, query_output = pdb_to_3di_with_output(query_pdb, encoding_model_dir, OUTPUT_DIR)
            target_seq, _, target_output = pdb_to_3di_with_output(target_pdb, encoding_model_dir, OUTPUT_DIR)

            # Log PDB fasta, length, and 3Di sequence
            logging.info(f"Processing pair: {query} vs {target}, encoding: {encoding}")
            logging.info(f"Query PDB: {query_pdb}, Length: {len(query_seq)}, 3Di Sequence: {query_seq}, Output: {query_output}")
            logging.info(f"Target PDB: {target_pdb}, Length: {len(target_seq)}, 3Di Sequence: {target_seq}, Output: {target_output}")

            # Compute SSW scores using the binary
            score_x1 = compute_ssw_score_with_binary(ssw_binary, query_seq, target_seq, gap_open=8, gap_extend=2)
            doubled_query_seq = query_seq * 2  # Double the query sequence for score_x2
            score_x2 = compute_ssw_score_with_binary(ssw_binary, doubled_query_seq, target_seq, gap_open=8, gap_extend=2)

            # Log SSW results with sequence details
            logging.info(f"SSW Results - Score X1: {score_x1}, Query Sequence: {query_seq}, Target Sequence: {target_seq}")
            logging.info(f"SSW Results - Score X2: {score_x2}, Doubled Query Sequence: {doubled_query_seq}, Target Sequence: {target_seq}")

            results.append({
                "query": query,
                "target": target,
                "encoding": encoding,
                "query_output": query_output,
                "target_output": target_output,
                "score_x1": score_x1,
                "score_x2": score_x2,
                "score_diff": score_x2 - score_x1
            })
        except Exception as e:
            logging.error(f"Error processing pair {query} vs {target} with encoding {encoding}: {str(e)}")
    
    return results

# Read pairs from the provided file
def read_pairs_from_file(pair_file):
    """Read query-target pairs from a file, skipping the first row."""
    df = pd.read_csv(pair_file, sep="\t", header=None, usecols=[0, 1], names=["query", "target"], dtype=str, skiprows=1)  # Skip the first row
    df = df.drop_duplicates(subset=["query", "target"])
    return df

def main():
    # Read and deduplicate input pairs
    pairs = read_scopealigngood(INPUT_FILE)
    # Update main function to use the pair file
    pairs = read_pairs_from_file(args.pair_file)

    all_results = []
    for _, row in pairs.iterrows():
        query, target = row["query"], row["target"]
        results = process_pair(query, target, PDB_DIR, MODEL_DIR, SSW_BINARY)
        all_results.extend(results)

    # Save results to CSV
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(os.path.join(OUTPUT_DIR, "comparison_results.csv"), index=False)

if __name__ == "__main__":
    main()