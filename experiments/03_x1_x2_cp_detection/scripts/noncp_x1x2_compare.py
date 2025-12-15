#!/usr/bin/env python3
"""
Compute SSW scores for non-CP candidate pairs with X1 (single) vs X2 (query doubled).

Outputs a CSV with columns:
  query,target,score_x1,score_x2,score_diff

Usage:
  python scripts/noncp_x1x2_compare.py \
    --pairs noncp_work/candidate_noncp_pairs.tsv \
    --pdb-dir scope_pdb \
    --out noncp_work/noncp_x1x2_scores.csv \
    --matrix foldseek_10f/s_10f.mat \
    --ssw-binary ssw/tmp/ssw/src/ssw_test \
    --foldseek-bin foldseek

Notes:
- Only query is doubled for X2; target remains single.
- Uses foldseek createdb to extract official 10f 3Di sequences from PDB.
- Gap penalties: open=8, extend=2 (same as prior runs).
"""

import argparse
import os
import re
import subprocess
from pathlib import Path
import pandas as pd


def run_cmd(cmd, cwd=None):
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(cmd)}\nSTDERR:{result.stderr}")
    return result.stdout


def extract_3di_with_foldseek(pdb_path: Path, foldseek_bin: str, tmp_root: Path) -> str:
    """Run foldseek createdb on a single PDB and return the 3Di sequence."""
    tmp_dir = tmp_root / pdb_path.stem
    tmp_dir.mkdir(parents=True, exist_ok=True)

    pdb_tmp = tmp_dir / pdb_path.name
    if pdb_tmp.resolve() != pdb_path.resolve():
        pdb_tmp.write_bytes(pdb_path.read_bytes())

    prefix = f"{pdb_path.stem}_pdzdb"
    cmd = [foldseek_bin, "createdb", str(pdb_tmp), prefix]
    run_cmd(cmd, cwd=tmp_dir)

    ss_candidates = [tmp_dir / f"{prefix}_ss", tmp_dir / f"{prefix}_pdzdb_ss"]
    ss_file = next((p for p in ss_candidates if p.exists()), None)
    if not ss_file:
        raise FileNotFoundError(f"No _ss file found after createdb for {pdb_path}")

    txt = ss_file.read_text().splitlines()
    seq_lines = [ln.strip() for ln in txt if ln and not ln.startswith(">")]
    if not seq_lines:
        raise ValueError(f"Empty 3Di sequence for {pdb_path}")
    return "".join(seq_lines)


def write_fasta(path: Path, header: str, seq: str):
    path.write_text(f">{header}\n{seq}\n")


def run_ssw(ssw_bin: str, matrix: Path, query_fa: Path, target_fa: Path, gap_open=8, gap_ext=2) -> int:
    cmd = [ssw_bin, "-o", str(gap_open), "-e", str(gap_ext), "-a", str(matrix), "-p", "-c", str(query_fa), str(target_fa)]
    out = run_cmd(cmd)
    # ssw_test emits SAM when -p -c is used; alignment score lives in the AS:i tag.
    m = re.search(r"AS:i:([-]?\d+)", out)
    if not m:
        # Fallback to legacy output format if ever used.
        m = re.search(r"optimal_alignment_score:\s*(\d+)", out)
    if not m:
        raise ValueError("Could not parse SSW score")
    return int(m.group(1))


def process_pair(row, pdb_dir: Path, tmp_dir: Path, foldseek_bin: str, ssw_bin: str, matrix: Path):
    qid = row["query_id"]
    tid = row["target_id"]

    q_pdb = pdb_dir / f"{qid.lower()}.pdb"
    t_pdb = pdb_dir / f"{tid.lower()}.pdb"
    if not q_pdb.exists():
        raise FileNotFoundError(f"Missing PDB: {q_pdb}")
    if not t_pdb.exists():
        raise FileNotFoundError(f"Missing PDB: {t_pdb}")

    # Extract 3Di sequences (10f official)
    q_seq = extract_3di_with_foldseek(q_pdb, foldseek_bin, tmp_dir)
    t_seq = extract_3di_with_foldseek(t_pdb, foldseek_bin, tmp_dir)

    # X1
    q_fa = tmp_dir / "query.fasta"
    t_fa = tmp_dir / "target.fasta"
    write_fasta(q_fa, f"{qid}_x1", q_seq)
    write_fasta(t_fa, f"{tid}", t_seq)
    score_x1 = run_ssw(ssw_bin, matrix, q_fa, t_fa)

    # X2 (query doubled)
    q2_fa = tmp_dir / "query_x2.fasta"
    write_fasta(q2_fa, f"{qid}_x2", q_seq + q_seq)
    score_x2 = run_ssw(ssw_bin, matrix, q2_fa, t_fa)

    # Clean small temps
    for p in [q_fa, q2_fa, t_fa]:
        try:
            p.unlink()
        except Exception:
            pass

    return score_x1, score_x2


def main():
    ap = argparse.ArgumentParser(description="SSW X1 vs X2 for non-CP pairs (10f)")
    ap.add_argument("--pairs", required=True, help="TSV with columns query_id, target_id (header required)")
    ap.add_argument("--pdb-dir", default="scope_pdb", help="Directory containing PDB files")
    ap.add_argument("--out", default="noncp_work/noncp_x1x2_scores.csv", help="Output CSV path")
    ap.add_argument("--matrix", default="encoders_and_tools/training_3di_gpu_10f/s_10f.mat", help="SSW substitution matrix")
    ap.add_argument("--ssw-binary", default="ssw/tmp/ssw/src/ssw_test", help="Path to ssw_test binary")
    ap.add_argument("--foldseek-bin", default="foldseek", help="Path to foldseek binary")
    ap.add_argument("--tmp-dir", default="noncp_work/tmp", help="Temp workspace for createdb/fastas")
    args = ap.parse_args()

    pairs_df = pd.read_csv(args.pairs, sep="\t")
    required_cols = {"query_id", "target_id"}
    if not required_cols.issubset(pairs_df.columns):
        raise ValueError(f"Input TSV must contain columns: {required_cols}")

    pdb_dir = Path(args.pdb_dir).resolve()
    tmp_dir = Path(args.tmp_dir).resolve()
    tmp_dir.mkdir(parents=True, exist_ok=True)
    matrix = Path(args.matrix).resolve()
    ssw_bin = str(Path(args.ssw_binary).resolve())
    foldseek_bin = args.foldseek_bin

    results = []
    for _, row in pairs_df.iterrows():
        try:
            score_x1, score_x2 = process_pair(row, pdb_dir, tmp_dir, foldseek_bin, ssw_bin, matrix)
            results.append({
                "query": row["query_id"],
                "target": row["target_id"],
                "score_x1": score_x1,
                "score_x2": score_x2,
                "score_diff": score_x2 - score_x1
            })
        except Exception as e:
            results.append({
                "query": row["query_id"],
                "target": row["target_id"],
                "score_x1": None,
                "score_x2": None,
                "score_diff": None,
                "error": str(e)
            })

    out_df = pd.DataFrame(results)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, index=False)
    print(f"Wrote {len(out_df)} rows to {args.out}")


if __name__ == "__main__":
    main()
