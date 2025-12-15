#!/usr/bin/env python3
import os
import sys
import subprocess
import shutil
import re
from pathlib import Path

def _resolve_foldseek_bin() -> str:
    # Prefer env override
    cand = os.environ.get('FOLDSEEK_BIN')
    if cand and Path(cand).exists():
        return cand
    # Try PATH
    from shutil import which
    w = which('foldseek')
    if w:
        return w
    # Fallback to common conda path used in this workspace
    fallback = '/home/jugipalace/miniconda3/envs/torchclean/bin/foldseek'
    if Path(fallback).exists():
        return fallback
    return 'foldseek'

FOLDSEEK_BIN = _resolve_foldseek_bin()

def cleanse_createdb_file(path: Path) -> str:
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
    lines = [l.strip() for l in content.splitlines() if l.strip()]
    if not lines:
        return ''
    if lines[0].startswith('>'):
        return ''.join(lines[1:])
    return ''.join(lines)


def run_createdb_get_3di(pdb_path: Path, work_dir: Path) -> str:
    base = pdb_path.stem
    work_dir.mkdir(parents=True, exist_ok=True)
    pdb_tmp = work_dir / pdb_path.name
    shutil.copy2(str(pdb_path), str(pdb_tmp))

    prefix = f"{base}_fsdb"
    cmd = [FOLDSEEK_BIN, 'createdb', str(pdb_tmp), prefix]
    res = subprocess.run(cmd, cwd=str(work_dir), capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"foldseek createdb failed ({res.returncode}): {res.stderr.strip()}")

    # locate *_ss file
    candidates = [work_dir / f"{prefix}_ss", work_dir / f"{prefix}_fsdb_ss"]
    ss_file = None
    for c in candidates:
        if c.exists():
            ss_file = c
            break
    if ss_file is None:
        # fallback: first file ending with _ss
        for f in work_dir.iterdir():
            if f.name.endswith('_ss'):
                ss_file = f
                break
    if ss_file is None or not ss_file.exists():
        raise FileNotFoundError('could not find createdb _ss file')

    text = cleanse_createdb_file(ss_file)
    seq = extract_sequence_from_file_content(text)
    if not seq:
        raise ValueError('empty 3Di sequence extracted')
    return seq


def main():
    if len(sys.argv) != 3:
        print('usage: extract_3di_from_createdb.py <pdb_path> <out_fasta>', file=sys.stderr)
        sys.exit(2)
    pdb_path = Path(sys.argv[1]).resolve()
    out_fa = Path(sys.argv[2]).resolve()
    work_dir = out_fa.parent / '.createdb_tmp' / pdb_path.stem

    try:
        seq = run_createdb_get_3di(pdb_path, work_dir)
    finally:
        # best-effort cleanup
        try:
            shutil.rmtree(work_dir)
        except Exception:
            pass

    out_fa.write_text(f'>{pdb_path.stem}\n{seq}\n')

if __name__ == '__main__':
    main()
