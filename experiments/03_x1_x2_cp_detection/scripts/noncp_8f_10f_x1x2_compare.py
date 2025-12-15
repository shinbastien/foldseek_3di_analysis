#!/usr/bin/env python3
"""
Compute X1/X2 SSW scores for non-CP pairs using pre-generated 3Di sequences (8f and 10f).

Reads:
  - noncp_pairs.tsv: pairs with query_id, target_id
  - Pre-generated 3Di sequences in fasta format (8f and 10f)

Outputs:
  - CSV with query, target, score_x1_8f, score_x2_8f, score_diff_8f,
                              score_x1_10f, score_x2_10f, score_diff_10f
"""

import argparse
import re
import subprocess
import logging
import sys
from pathlib import Path
import pandas as pd


def setup_logger(debug=False, log_file=None):
    """Setup logging with optional file output."""
    logger = logging.getLogger(__name__)
    level = logging.DEBUG if debug else logging.INFO
    logger.setLevel(level)
    
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


def run_cmd(cmd, cwd=None):
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(cmd)}\nSTDERR:{result.stderr}")
    return result.stdout


def write_fasta(path: Path, header: str, seq: str):
    path.write_text(f">{header}\n{seq}\n")


def run_ssw(ssw_bin: str, matrix: Path, query_fa: Path, target_fa: Path, gap_open=8, gap_ext=2) -> int:
    cmd = [ssw_bin, "-o", str(gap_open), "-e", str(gap_ext), "-a", str(matrix), "-p", "-c", str(query_fa), str(target_fa)]
    out = run_cmd(cmd)
    # Parse SAM format output: alignment score in AS:i tag
    m = re.search(r"AS:i:([-]?\d+)", out)
    if not m:
        raise ValueError(f"Could not parse SSW score from output")
    return int(m.group(1))


def load_3di_dict(fasta_path: Path):
    """Load 3Di sequences from fasta file. Remove NULL characters."""
    seqs = {}
    with open(fasta_path, 'rb') as f:
        current_header = None
        current_seq = []
        for line in f:
            # Decode and remove NULL characters
            line = line.decode('utf-8', errors='ignore').strip()
            line = line.replace('\x00', '')  # Remove NULL character
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


def process_pair(qid: str, tid: str, seqs_8f: dict, seqs_10f: dict,
                 tmp_dir: Path, ssw_bin: str, matrix_8f: Path, matrix_10f: Path, 
                 pair_num: int = 0, logger=None, verbose=False):
    """Process one pair with both 8f and 10f. Optional verbose logging."""
    q_seq_8f = seqs_8f.get(qid)
    t_seq_8f = seqs_8f.get(tid)
    q_seq_10f = seqs_10f.get(qid)
    t_seq_10f = seqs_10f.get(tid)

    if not q_seq_8f or not t_seq_8f or not q_seq_10f or not t_seq_10f:
        raise ValueError(f"Missing sequence for {qid} or {tid}")

    results = {}

    if verbose and logger:
        logger.info(f"\n{'='*80}")
        logger.info(f"Pair #{pair_num}: {qid} vs {tid}")
        logger.info(f"{'='*80}")
        logger.info(f"[8f Alphabet]")
        logger.info(f"  X1 Query ({qid}): {q_seq_8f[:60]}{'...' if len(q_seq_8f) > 60 else ''}")
        logger.info(f"  X1 Target ({tid}): {t_seq_8f[:60]}{'...' if len(t_seq_8f) > 60 else ''}")

    # 8f X1
    q_fa = tmp_dir / "query.fasta"
    t_fa = tmp_dir / "target.fasta"
    write_fasta(q_fa, f"{qid}_x1", q_seq_8f)
    write_fasta(t_fa, f"{tid}", t_seq_8f)
    results['score_x1_8f'] = run_ssw(ssw_bin, matrix_8f, q_fa, t_fa)
    
    if verbose and logger:
        logger.info(f"  X1 SSW Score: {results['score_x1_8f']}")

    # 8f X2
    q2_fa = tmp_dir / "query_x2.fasta"
    write_fasta(q2_fa, f"{qid}_x2", q_seq_8f + q_seq_8f)
    results['score_x2_8f'] = run_ssw(ssw_bin, matrix_8f, q2_fa, t_fa)
    results['score_diff_8f'] = results['score_x2_8f'] - results['score_x1_8f']
    
    if verbose and logger:
        logger.info(f"  X2 Query ({qid}_x2): {(q_seq_8f + q_seq_8f)[:60]}{'...' if len(q_seq_8f + q_seq_8f) > 60 else ''}")
        logger.info(f"  X2 SSW Score: {results['score_x2_8f']}")
        logger.info(f"  Diff (X2-X1): {results['score_diff_8f']:+d}")

    # 10f X1
    if verbose and logger:
        logger.info(f"\n  [10f Alphabet]")
        logger.info(f"  X1 Query ({qid}): {q_seq_10f[:60]}{'...' if len(q_seq_10f) > 60 else ''}")
        logger.info(f"  X1 Target ({tid}): {t_seq_10f[:60]}{'...' if len(t_seq_10f) > 60 else ''}")
    
    write_fasta(q_fa, f"{qid}_x1", q_seq_10f)
    write_fasta(t_fa, f"{tid}", t_seq_10f)
    results['score_x1_10f'] = run_ssw(ssw_bin, matrix_10f, q_fa, t_fa)
    
    if verbose and logger:
        logger.info(f"  X1 SSW Score: {results['score_x1_10f']}")

    # 10f X2
    write_fasta(q2_fa, f"{qid}_x2", q_seq_10f + q_seq_10f)
    results['score_x2_10f'] = run_ssw(ssw_bin, matrix_10f, q2_fa, t_fa)
    results['score_diff_10f'] = results['score_x2_10f'] - results['score_x1_10f']
    
    if verbose and logger:
        logger.info(f"  X2 Query ({qid}_x2): {(q_seq_10f + q_seq_10f)[:60]}{'...' if len(q_seq_10f + q_seq_10f) > 60 else ''}")
        logger.info(f"  X2 SSW Score: {results['score_x2_10f']}")
        logger.info(f"  Diff (X2-X1): {results['score_diff_10f']:+d}")

    # Clean
    for p in [q_fa, q2_fa, t_fa]:
        try:
            p.unlink()
        except:
            pass

    return results


def generate_x1_x2_fasta(pairs_df, seqs_8f, seqs_10f, out_base_8f: Path, out_base_10f: Path, logger=None):
    """Generate separate X1 and X2 FASTA files for structured comparison."""
    if logger:
        logger.info(f"Generating X1/X2 FASTA files for 8f alphabet...")
    
    # Collect all unique sequences
    x1_seqs_8f = {}
    x2_seqs_8f = {}
    x1_seqs_10f = {}
    x2_seqs_10f = {}
    
    for _, row in pairs_df.iterrows():
        qid = row["query_id"]
        tid = row["target_id"]
        
        # 8f
        if qid in seqs_8f:
            x1_seqs_8f[qid] = seqs_8f[qid]
            x2_seqs_8f[f"{qid}_x2"] = seqs_8f[qid] + seqs_8f[qid]
        if tid in seqs_8f:
            x1_seqs_8f[tid] = seqs_8f[tid]
        
        # 10f
        if qid in seqs_10f:
            x1_seqs_10f[qid] = seqs_10f[qid]
            x2_seqs_10f[f"{qid}_x2"] = seqs_10f[qid] + seqs_10f[qid]
        if tid in seqs_10f:
            x1_seqs_10f[tid] = seqs_10f[tid]
    
    # Write X1 files
    x1_out_8f = out_base_8f / "X1_sequences.fasta"
    with open(x1_out_8f, 'w') as f:
        for header, seq in sorted(x1_seqs_8f.items()):
            f.write(f">{header}\n{seq}\n")
    if logger:
        logger.info(f"  Wrote {len(x1_seqs_8f)} X1 sequences (8f) to {x1_out_8f}")
    
    x1_out_10f = out_base_10f / "X1_sequences.fasta"
    with open(x1_out_10f, 'w') as f:
        for header, seq in sorted(x1_seqs_10f.items()):
            f.write(f">{header}\n{seq}\n")
    if logger:
        logger.info(f"  Wrote {len(x1_seqs_10f)} X1 sequences (10f) to {x1_out_10f}")
    
    # Write X2 files
    x2_out_8f = out_base_8f / "X2_sequences.fasta"
    with open(x2_out_8f, 'w') as f:
        for header, seq in sorted(x2_seqs_8f.items()):
            f.write(f">{header}\n{seq}\n")
    if logger:
        logger.info(f"  Wrote {len(x2_seqs_8f)} X2 sequences (8f) to {x2_out_8f}")
    
    x2_out_10f = out_base_10f / "X2_sequences.fasta"
    with open(x2_out_10f, 'w') as f:
        for header, seq in sorted(x2_seqs_10f.items()):
            f.write(f">{header}\n{seq}\n")
    if logger:
        logger.info(f"  Wrote {len(x2_seqs_10f)} X2 sequences (10f) to {x2_out_10f}")
    
    if logger:
        logger.info("X1/X2 FASTA generation complete")


def main():
    ap = argparse.ArgumentParser(description="X1/X2 SSW for noncp pairs (8f + 10f)")
    ap.add_argument("--pairs", required=True, help="TSV with query_id, target_id")
    ap.add_argument("--fasta-8f", required=True, help="Fasta file with 8f 3Di sequences")
    ap.add_argument("--fasta-10f", required=True, help="Fasta file with 10f 3Di sequences")
    ap.add_argument("--out", default="noncp_8f_10f_comparison.csv", help="Output CSV")
    ap.add_argument("--matrix-8f", default="encoders_and_tools/training_3di_gpu_8f/s_8f.mat", help="8f substitution matrix")
    ap.add_argument("--matrix-10f", default="encoders_and_tools/training_3di_gpu_10f/s_10f.mat", help="10f substitution matrix")
    ap.add_argument("--ssw-binary", default="ssw/tmp/ssw/src/ssw_test", help="SSW binary")
    ap.add_argument("--tmp-dir", default="tmp_x1x2", help="Temp directory")
    ap.add_argument("--debug", action="store_true", help="Enable debug mode with X1/X2 FASTA generation and logging")
    ap.add_argument("--verbose", action="store_true", help="Enable verbose logging of each pair's sequences and scores")
    args = ap.parse_args()

    # Setup logging if debug mode
    logger = None
    if args.debug or args.verbose:
        log_file = Path(args.out).parent / "processing.log"
        logger = setup_logger(debug=True, log_file=str(log_file))
        logger.info("=== X1/X2 Comparison with Debug Logging ===")
        logger.info(f"Pairs TSV: {args.pairs}")
        logger.info(f"8f FASTA: {args.fasta_8f}")
        logger.info(f"10f FASTA: {args.fasta_10f}")
        logger.info(f"Output CSV: {args.out}")
        logger.info(f"SSW binary: {args.ssw_binary}")
        if args.verbose:
            logger.info(f"Verbose mode: ON (detailed sequence logging)")

    pairs_df = pd.read_csv(args.pairs, sep="\t")
    matrix_8f = Path(args.matrix_8f).resolve()
    matrix_10f = Path(args.matrix_10f).resolve()
    ssw_bin = str(Path(args.ssw_binary).resolve())
    tmp_dir = Path(args.tmp_dir).resolve()
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Load sequences
    msg = f"Loading 8f sequences from {args.fasta_8f}..."
    if logger:
        logger.info(msg)
    print(msg, flush=True)
    seqs_8f = load_3di_dict(Path(args.fasta_8f))
    msg = f"  Loaded {len(seqs_8f)} sequences"
    if logger:
        logger.info(msg)
    print(msg, flush=True)

    msg = f"Loading 10f sequences from {args.fasta_10f}..."
    if logger:
        logger.info(msg)
    print(msg, flush=True)
    seqs_10f = load_3di_dict(Path(args.fasta_10f))
    msg = f"  Loaded {len(seqs_10f)} sequences"
    if logger:
        logger.info(msg)
    print(msg, flush=True)

    # Generate X1/X2 FASTA if debug mode
    if args.debug:
        out_base = Path(args.out).parent
        out_base_8f = out_base / "8f" / "3di_sequences"
        out_base_10f = out_base / "10f" / "3di_sequences"
        out_base_8f.mkdir(parents=True, exist_ok=True)
        out_base_10f.mkdir(parents=True, exist_ok=True)
        generate_x1_x2_fasta(pairs_df, seqs_8f, seqs_10f, out_base_8f, out_base_10f, logger)

    results = []
    for idx, row in pairs_df.iterrows():
        qid = row["query_id"]
        tid = row["target_id"]
        try:
            res = process_pair(qid, tid, seqs_8f, seqs_10f, tmp_dir, ssw_bin, matrix_8f, matrix_10f,
                             pair_num=idx+1, logger=logger, verbose=args.verbose)
            res['query'] = qid
            res['target'] = tid
            results.append(res)
            if logger and (idx + 1) % 100 == 0:
                logger.debug(f"Processed {idx + 1}/{len(pairs_df)}")
        except Exception as e:
            results.append({
                'query': qid,
                'target': tid,
                'error': str(e)
            })
            if logger:
                logger.warning(f"Error processing {qid}-{tid}: {str(e)}")

        if (idx + 1) % 100 == 0:
            print(f"Processed {idx + 1}/{len(pairs_df)}...", flush=True)

    out_df = pd.DataFrame(results)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    msg = f"Wrote {len(out_df)} rows to {out_path}"
    if logger:
        logger.info(msg)
    print(msg, flush=True)


if __name__ == "__main__":
    main()
