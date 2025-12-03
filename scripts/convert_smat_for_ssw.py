#!/usr/bin/env python3
"""Convert an existing s.mat-like file into a ssw_test-compatible matrix.

Behavior:
 - Reads input matrix that may have a header line and rows like: "A 4 -1 0 ..."
 - Produces a normalized matrix where first line is the header letters separated by single spaces
   (no leading spaces), and subsequent lines are: LETTER <scores...>
 - If --letters is provided, the output header/order will follow that list (only those letters).
 - If input lacks scores for some letter pairs, fallback to simple match/mismatch defaults.

Usage:
  scripts/convert_smat_for_ssw.py -i training_3di_gpu_new/tmp/s.mat -o /tmp/converted_s.mat \
      --letters ABCFGHIJKLMNPQRTX

This helps avoid "Problem of reading the weight matrix file." errors from ssw_test by
ensuring a strict whitespace-normalized format and a header matching the sequences.
"""
import argparse
from pathlib import Path
import sys


def read_matrix(path: Path):
    txt = path.read_text().splitlines()
    lines = [l for l in (t.strip() for t in txt) if l]
    if not lines:
        raise RuntimeError('Empty matrix file')
    first = lines[0].split()
    header = None
    raw_rows = {}
    # Detect header when first tokens are single-letter tokens
    if all(tok.isalpha() and len(tok) == 1 for tok in first):
        header = first
        for ln in lines[1:]:
            parts = ln.split()
            if not parts:
                continue
            if parts[0].isalpha() and len(parts[0]) == 1:
                key = parts[0]
                vals = parts[1:]
                try:
                    raw_rows[key] = [int(x) for x in vals]
                except ValueError:
                    # ignore malformed numeric tokens
                    pass
    else:
        # fallback: try to detect any alphabetic tokens in first line as header
        if any(tok.isalpha() for tok in first):
            header = [tok for tok in first if tok.isalpha() and len(tok) == 1]
        # parse all lines that look like rows
        for ln in lines:
            parts = ln.split()
            if not parts:
                continue
            if parts[0].isalpha() and len(parts[0]) == 1:
                key = parts[0]
                vals = parts[1:]
                try:
                    raw_rows[key] = [int(x) for x in vals]
                except ValueError:
                    pass

    return header or [], raw_rows


def write_matrix(out_path: Path, letters, raw_rows, default_match=4, default_mismatch=-1):
    # Build a dict-of-dicts from raw_rows (which maps row_letter -> list of values aligned with input header)
    # raw_rows may have been parsed with an input header; we need to align columns to letter positions.
    # First try to infer input header length from any row lengths
    inferred_header = None
    for vals in raw_rows.values():
        if inferred_header is None or len(vals) > len(inferred_header or []):
            inferred_header = ['?'] * len(vals)

    # Create matrix mapping a->b->score
    mat = {a: {} for a in letters}

    # If raw_rows contains lists, attempt to map them to letters using best-effort ordering:
    # If raw_rows length equals len(letters), assume columns match 'letters' order; otherwise align by index if possible.
    for r, vals in raw_rows.items():
        # if vals length equals letters length, map directly
        if len(vals) >= len(letters):
            for i, b in enumerate(letters):
                try:
                    mat.setdefault(r, {})[b] = int(vals[i])
                except Exception:
                    pass
        else:
            # shorter rows: map as far as possible
            for i, v in enumerate(vals):
                if i < len(letters):
                    mat.setdefault(r, {})[letters[i]] = int(v)

    # Symmetrize: for each pair (i,j), if one side missing copy the other; if both present but differ, average
    for i in letters:
        for j in letters:
            a_exists = j in mat.get(i, {})
            b_exists = i in mat.get(j, {})
            if a_exists and b_exists:
                ai = mat[i][j]
                bi = mat[j][i]
                if ai != bi:
                    avg = int(round((ai + bi) / 2.0))
                    mat[i][j] = mat[j][i] = avg
            elif a_exists and not b_exists:
                mat.setdefault(j, {})[i] = mat[i][j]
            elif (not a_exists) and b_exists:
                mat.setdefault(i, {})[j] = mat[j][i]
            else:
                # neither present: fill diagonal with match and others with mismatch
                if i == j:
                    mat[i][j] = default_match
                else:
                    mat[i][j] = default_mismatch

    # Write normalized matrix
    with out_path.open('w') as f:
        f.write(' '.join(letters) + '\n')
        for a in letters:
            vals = [str(int(mat[a][b])) for b in letters]
            f.write(a + ' ' + ' '.join(vals) + '\n')


def main():
    p = argparse.ArgumentParser(description='Convert s.mat to ssw_test-friendly matrix')
    p.add_argument('-i', '--input', required=True)
    p.add_argument('-o', '--output', required=False,
                   help='Output path. If omitted or outside the repo, will be written to repo root as converted_s.mat')
    p.add_argument('--pairwise-dir', required=False,
                   help='Optional: path to pairwise run dir where converted matrix should be written. If provided, output will be written inside this dir with an auto-chosen name (converted_8f_s.mat or converted_10f_s.mat).')
    p.add_argument('--letters', type=str, default=None,
                   help='Optional string of letters to include (e.g. ABCDE). If omitted, header of input used.')
    p.add_argument('--match', type=int, default=4, help='Default match score')
    p.add_argument('--mismatch', type=int, default=-1, help='Default mismatch score')
    args = p.parse_args()

    inp = Path(args.input)
    out = Path(args.output) if args.output else None
    if not inp.exists():
        print('Input matrix not found:', inp, file=sys.stderr)
        sys.exit(2)
    # Diagnostic: print file info for debugging unreadable/absent s.mat
    try:
        st = inp.stat()
        print(f"[DEBUG] convert_smat_for_ssw: input exists={inp.exists()}, abs={inp.resolve()}, size={st.st_size}")
    except Exception:
        print(f"[DEBUG] convert_smat_for_ssw: input exists={inp.exists()}, abs={inp if inp else 'N/A'}")
    header, rows = read_matrix(inp)
    # Diagnostic: report what was parsed
    try:
        print(f"[DEBUG] convert_smat_for_ssw: parsed header={header}")
        print(f"[DEBUG] convert_smat_for_ssw: parsed rows_count={len(rows)}")
        sample_keys = list(rows.keys())[:5]
        print(f"[DEBUG] convert_smat_for_ssw: sample rows keys={sample_keys}")
    except Exception:
        pass
    if args.letters:
        letters = list(args.letters)
    else:
        if header:
            letters = header
        else:
            # fallback: use row keys
            letters = sorted(rows.keys())
    # Ensure single-character letters and keep order
    letters = [l for l in letters if isinstance(l, str) and len(l) == 1]

    # Determine repository root (parent of scripts/)
    # Ensure single-character letters and keep order
    letters = [l for l in letters if isinstance(l, str) and len(l) == 1]

    # If pairwise-dir supplied, place the converted file inside it with a name depending on input origin
    if args.pairwise_dir:
        pair_dir = Path(args.pairwise_dir)
        pair_dir.mkdir(parents=True, exist_ok=True)
        inp_lower = str(inp).lower()
        # Heuristics to determine new vs orig
        if 'training_3di_gpu_new' in inp_lower or '/new/' in inp_lower or 'new' in inp_lower:
            outname = 'converted_8f_s.mat'
        elif 'training_3di_gpu_original' in inp_lower or '/original/' in inp_lower or 'orig' in inp_lower:
            outname = 'converted_10f_s.mat'
        else:
            # fallback: if input filename contains '8' or '8f' assume new
            if '8f' in inp_lower or '8_f' in inp_lower or '8' in inp_lower:
                outname = 'converted_8f_s.mat'
            elif '10f' in inp_lower or '10_f' in inp_lower or '10' in inp_lower:
                outname = 'converted_10f_s.mat'
            else:
                outname = 'converted_s.mat'
        out = pair_dir / outname
    else:
        repo_root = Path(__file__).resolve().parents[1]
        if out is None:
            out = repo_root / 'converted_s.mat'
        else:
            # if user supplied an absolute path outside the repo, rewrite to repo root with same basename
            try:
                out = out.resolve()
                if repo_root not in out.parents:
                    out = repo_root / out.name
            except Exception:
                out = repo_root / (out.name if out is not None else 'converted_s.mat')

    write_matrix(out, letters, rows, default_match=args.match, default_mismatch=args.mismatch)
    print('Wrote converted matrix to', out)


if __name__ == '__main__':
    main()
