#!/usr/bin/env python3
"""
Utility functions for PDZ pipeline (TM-align, SSW, and 3Di sequence handling).

These functions are shared across:
- pairwise_3di_pipeline.py (4_trimmed_ssw_comparison)
- tmalign_3di_match_pipeline.py (5_tm_ssw_coverage_analysis)
"""

import re
import subprocess
from pathlib import Path
from typing import Optional, Tuple, List
from Bio.PDB import PDBParser


def run_cmd(cmd: List[str], cwd: Optional[str] = None, capture: bool = False) -> Tuple[int, str, str]:
    """Run a command and return (returncode, stdout, stderr)."""
    print('> ' + ' '.join(cmd))
    try:
        proc = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE if capture else None,
                              stderr=subprocess.PIPE if capture else None, check=False, text=True)
        stdout = proc.stdout if capture else ''
        stderr = proc.stderr if capture else ''
        return proc.returncode, stdout, stderr
    except FileNotFoundError as e:
        return 127, '', str(e)


def write_fasta_from_string(target_path: Path, header: str, seq: str):
    """Write a FASTA sequence to file."""
    target_path.write_text(f">{header}\n{seq}\n")


def read_fasta_sequence(path: Path) -> str:
    """Read first FASTA sequence from file and return sequence string."""
    try:
        text = path.read_text()
    except Exception:
        return ''
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if not lines:
        return ''
    if lines[0].startswith('>'):
        return ''.join(lines[1:])
    return ''.join(lines)


def make_x2_fasta(fasta_path: Path, tmpdir: Path) -> Path:
    """Create a doubled-sequence FASTA (seq concatenated twice) for SSW query doubling.

    Returns path to new FASTA file under tmpdir.
    """
    seq = read_fasta_sequence(fasta_path)
    doubled = seq + seq
    out = tmpdir / (fasta_path.stem + '_x2.fasta')
    write_fasta_from_string(out, fasta_path.stem + '_x2', doubled)
    return out


def run_ssw_and_parse(ssw_bin: str, matrix: Optional[str], target_fasta: Path, query_fasta: Path,
                      gap_open: int = 11, gap_ext: int = 1, min_score: int = 0,
                      sam: bool = False, cwd: Optional[str] = None) -> Tuple[int, str]:
    """Run ssw_test and return (returncode, stdout).

    If `sam` is True, include `-s` to request SAM output. If `matrix` is None,
    ssw_test will be invoked without `-a` to use the default BLOSUM.
    The command is executed with working directory `cwd` if provided.
    """
    cmd = [ssw_bin, '-o', str(gap_open), '-e', str(gap_ext)]
    if matrix:
        cmd += ['-a', matrix]
    if sam:
        cmd += ['-s']
    cmd += ['-p', '-f', str(min_score), str(target_fasta), str(query_fasta)]
    rc, out, err = run_cmd(cmd, cwd=cwd, capture=True)
    if rc != 0:
        raise RuntimeError(f'ssw_test failed (rc={rc}): {err}')
    return rc, out


def get_residue_number_list(pdb_path: Path) -> List[int]:
    """Get list of residue numbers from PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('s', str(pdb_path))
    model = structure[0]
    chain = list(model.get_chains())[0]
    residues = [r for r in chain.get_residues() if len(r.id[0].strip()) == 0]
    return [r.id[1] for r in residues]


def resnum_range_to_indices(resnums: List[int], start: int, end: int) -> Tuple[int, int]:
    """Map residue numbers inclusive range [start,end] to 0-based indices over resnums list.
    
    Returns (i0, i1) as inclusive indices. If not found, raise ValueError.
    """
    try:
        i0 = resnums.index(start)
        i1 = resnums.index(end)
        return i0, i1
    except ValueError:
        # if exact residue numbers not in list, try to find first residue >= start and last <= end
        i0 = next((i for i, r in enumerate(resnums) if r >= start), None)
        i1 = next((i for i, r in enumerate(resnums) if r > end), None)
        if i0 is None:
            raise ValueError('Start residue number not found or out of range')
        if i1 is None:
            i1 = len(resnums) - 1
        else:
            i1 = i1 - 1
        if i0 > i1:
            raise ValueError('Computed indices invalid')
        return i0, i1


def parse_tmalign_ranges(tmalign_stdout: str) -> Optional[Tuple[Tuple[int, int], Tuple[int, int]]]:
    """Parse TM-align textual alignment block and return 0-based inclusive index ranges
    for the two aligned chains.

    This function searches for the alignment display (blocks of three lines: query, match line, subject)
    and concatenates them if wrapped. It then finds the first and last columns where both sequences
    have non-gap characters and returns their 0-based inclusive indices as ((q_start,q_end),(s_start,s_end)).

    Returns None if no suitable alignment block is found.
    """
    lines = [l.rstrip() for l in tmalign_stdout.splitlines()]

    q_aln_parts = []
    s_aln_parts = []

    i = 0
    while i + 2 < len(lines):
        top = lines[i].strip()
        mid = lines[i + 1].strip()
        bot = lines[i + 2].strip()

        # mid line typically contains ":" and/or "." and spaces indicating alignment quality
        if re.match(r'^[\.:\s]+$', mid) and re.search(r'[A-Za-z\-]', top) and re.search(r'[A-Za-z\-]', bot):
            # top and bot are alignment sequence fragments
            q_aln_parts.append(re.sub(r'\s+', '', top))
            s_aln_parts.append(re.sub(r'\s+', '', bot))
            i += 3
            # consume any blank line after a wrapped block
            while i < len(lines) and lines[i].strip() == '':
                i += 1
            continue
        i += 1

    if not q_aln_parts or not s_aln_parts:
        return None

    q_aln = ''.join(q_aln_parts)
    s_aln = ''.join(s_aln_parts)

    if len(q_aln) != len(s_aln):
        # alignment fragments should be same length
        return None

    # walk alignment columns to find first/last column where both are non-gaps
    q_pos = -1
    s_pos = -1
    q_start = None
    q_end = None
    s_start = None
    s_end = None

    for col_idx in range(len(q_aln)):
        qc = q_aln[col_idx]
        sc = s_aln[col_idx]
        if qc != '-':
            q_pos += 1
        if sc != '-':
            s_pos += 1

        if qc != '-' and sc != '-':
            if q_start is None:
                q_start = q_pos
                s_start = s_pos
            q_end = q_pos
            s_end = s_pos

    if q_start is None:
        return None

    return (q_start, q_end), (s_start, s_end)


def extract_sequence_from_file_content(content: str) -> str:
    """If content is FASTA-like, extract first sequence; otherwise return content stripped.
    
    This is used for pdzdb_ss which contains 3di sequence (maybe with null param at end)
    """
    lines = [l.strip() for l in content.splitlines() if l.strip()]
    if not lines:
        return ''
    if lines[0].startswith('>'):
        # join subsequent lines until next header
        seq = ''.join(lines[1:])
        return seq
    # otherwise assume whole file is sequence
    return ''.join(lines)


def find_createdb_outputs(tmpdir: Path, prefix: str) -> Tuple[Optional[Path], Optional[Path]]:
    """Find createdb outputs that contain fasta and ss (3di) info.

    Returns (pdzdb_path, pdzdb_ss_path) or (None,None) if not found.
    """
    pdz = None
    pdz_ss = None
    for p in tmpdir.iterdir():
        name = p.name
        if prefix in name and 'pdz' in name and 'ss' in name:
            pdz_ss = p
        elif prefix in name and 'pdz' in name:
            pdz = p
        # fallback: user said produced files end with _pdzdb and _pdzdb_ss
        if prefix + '_pdzdb' in name:
            pdz = p
        if prefix + '_pdzdb_ss' in name or name.endswith('_pdzdb_ss'):
            pdz_ss = p
    return pdz, pdz_ss


def cleanse_createdb_file(path: Path) -> str:
    """Remove trailing null parameter markers and null bytes.

    The user's instruction: remove trailing null parameter. We apply conservative cleaning:
    - strip NULL bytes
    - strip trailing whitespace
    - remove trailing tokens like '\t0', '\tNULL', ' 0' at line ends
    """
    raw = path.read_bytes()
    try:
        text = raw.decode('utf-8')
    except UnicodeDecodeError:
        # try latin-1 fallback
        text = raw.decode('latin-1')
    lines = []
    for ln in text.splitlines():
        s = ln.rstrip('\x00').rstrip()
        # remove trailing tab followed by 0 or NULL
        s = re.sub(r"(\t|\s)+(0|NULL)$", '', s)
        lines.append(s)
    return '\n'.join(lines) + '\n'


def normalize_matrix_file(src: Path, dst: Path):
    """Create a normalized copy of a score matrix suitable for ssw_test.

    Normalization strategy:
    - Strip leading/trailing whitespace on each line
    - Ensure the header row (letters) is a single whitespace-separated row
    - Keep numeric spacing normalized to single spaces
    """
    try:
        text = src.read_text()
    except Exception as e:
        # Provide diagnostics for why reading the provided matrix failed
        try:
            print(f"[DEBUG] normalize_matrix_file: failed to read source matrix {src}: exists={src.exists()}, abs={src.resolve() if src.exists() else 'N/A'}")
            try:
                st = src.stat()
                print(f"[DEBUG] normalize_matrix_file: file size={st.st_size}")
            except Exception:
                pass
        except Exception:
            pass
        print(f"[DEBUG] normalize_matrix_file: exception: {e}")
        raise

    lines = [l.strip() for l in text.splitlines() if l.strip()]
    norm_lines = []
    for ln in lines:
        # collapse multiple spaces/tabs into single spaces
        parts = ln.split()
        norm_lines.append(' '.join(parts))
    dst.write_text('\n'.join(norm_lines) + '\n')


def read_matrix_alphabet(matrix_path: Path):
    """Attempt to read the scoring matrix header row and return the set of alphabet symbols.

    The function looks for the first non-empty line and treats whitespace-separated
    tokens as header symbols if they are alphabetic-like (single character tokens).
    Returns a set of single-character uppercase symbols. If parsing fails, returns None.
    """
    try:
        text = matrix_path.read_text()
    except Exception:
        return None
    for ln in text.splitlines():
        ln = ln.strip()
        if not ln:
            continue
        parts = ln.split()
        letters = [p for p in parts if p.isalpha() and len(p) == 1]
        if letters:
            return set([p.upper() for p in letters])
        # fallback: if line looks like a continuous string of letters
        contiguous = ''.join(parts)
        if contiguous.isalpha() and 1 < len(contiguous) <= 64:
            return set([c.upper() for c in contiguous])
    return None


def parse_ssw_results(text: str) -> List[str]:
    """Very small parser: return non-empty lines as 'entries' to display."""
    return [l for l in (ln.strip() for ln in text.splitlines()) if l]
