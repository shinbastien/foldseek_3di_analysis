#!/usr/bin/env python3
"""
Pairwise 3Di pipeline

This script implements the 13-step pipeline described by the user.

Usage:
    python pairwise_3di_pipeline.py pdb1.pdb pdb2.pdb \
        --model-dir encoders_and_tools/training_3di_gpu_8f \
        [--outdir ./tmp_run] [--foldseek foldseek] [--tmalign TMalign] [--ssw_new tmp/ssw_test] [--ssw_orig tmp/ssw_test] [--yes]

Behavior summary (defaults):
 - foldseek binary: 'foldseek' (from PATH)
 - TMalign binary: 'TMalign' (from PATH)
 - ssw_test for new: encoders_and_tools/training_3di_gpu_8f/ssw_test
 - ssw_test for orig: encoders_and_tools/training_3di_gpu_10f/ssw_test (falls back to 8f)
 - model_dir default: encoders_and_tools/training_3di_gpu_8f
 - Non-interactive: use --yes

The script prints explanations before each step and asks for confirmation (press Enter) unless --yes is provided.

"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time
import re
from pathlib import Path
from typing import Optional, Tuple, List
import logging

# Add project root to path so we can import pdb_to_3di helper
ROOT = os.path.abspath(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Early logging: create a small run_logs dir and configure a logger immediately so
# messages printed before tmpdir creation are captured and visible on stderr.
from pathlib import Path as _Path_for_logger
_run_logs_dir = _Path_for_logger(ROOT) / 'run_logs'
try:
    _run_logs_dir.mkdir(parents=True, exist_ok=True)
except Exception:
    pass
_ts0 = time.strftime('%Y%m%d%H%M%S')
_start_log = _run_logs_dir / f'pairwise_3di_start_{_ts0}.log'
logger = logging.getLogger('pairwise_3di')
logger.setLevel(logging.INFO)
if not logger.handlers:
    try:
        fh = logging.FileHandler(str(_start_log))
        fh.setLevel(logging.INFO)
        sh = logging.StreamHandler(sys.stderr)
        sh.setLevel(logging.INFO)
        fmt = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
        fh.setFormatter(fmt)
        sh.setFormatter(fmt)
        logger.addHandler(fh)
        logger.addHandler(sh)
    except Exception:
        # best-effort: if logger setup fails, continue without file handler
        pass
logger.info(f"Starting pairwise_3di_pipeline; project ROOT={ROOT}; start_log={_start_log}")

import pdb_to_3di as p2  # reuse functions like validate_model_files/pdb_to_3di
from Bio.PDB import PDBParser

# Import utility functions from shared module
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))  # Add project root
from utils.pdz_pipeline_utils import (
    run_cmd, find_createdb_outputs, cleanse_createdb_file,
    write_fasta_from_string, extract_sequence_from_file_content,
    get_residue_number_list, resnum_range_to_indices,
    parse_tmalign_ranges, run_ssw_and_parse, read_fasta_sequence,
    make_x2_fasta, parse_ssw_results, normalize_matrix_file,
    read_matrix_alphabet
)

def confirm(message: str, auto_yes: bool) -> bool:
    # Use logger when available so messages appear early (before tmpdir/log redirection).
    try:
        logger.info(message)
    except Exception:
        print('\n' + message)

    if auto_yes:
        try:
            logger.info('[Auto-confirmed with --yes]')
        except Exception:
            print('[Auto-confirmed with --yes]')
        return True

    resp = input('Press [Enter] to continue, or type n to abort: ')
    if resp.strip().lower().startswith('n'):
        try:
            logger.info('Aborted by user.')
        except Exception:
            print('Aborted by user.')
        return False
    return True


def main():
    # Parse CLI args and determine which model variants to run
    logger.info('Reached CLI parsing point')
    parser = argparse.ArgumentParser(description='Pairwise 3Di pipeline')
    parser.add_argument('pdb1')
    parser.add_argument('pdb2')
    parser.add_argument('--model-dir', help='Single model tmp dir (override). If provided, run only this model; otherwise auto-detect 8f/9f/10f', default=None)
    parser.add_argument('--outdir', help='Output directory', default=None)
    parser.add_argument('--foldseek', default='foldseek')
    parser.add_argument('--tmalign', default='TMalign')
    parser.add_argument('--ssw_new', default='tmp/ssw_test')
    parser.add_argument('--ssw_orig', default='tmp/ssw_test')
    parser.add_argument('--ssw', default=None, help='Path to a single ssw_test binary to use for all variants (overrides per-variant binaries). Can also be set via env SSW_TEST')
    parser.add_argument('--gap-open', type=int, default=8)
    parser.add_argument('--gap-extend', type=int, default=2)
    parser.add_argument('--yes', action='store_true', dest='yes')
    args = parser.parse_args()
    logger.info(f'Parsed CLI args: pdb1={args.pdb1} pdb2={args.pdb2} model_dir={args.model_dir} yes={args.yes}')

    # expose frequently-used variables expected by the rest of the script
    pdb1 = Path(args.pdb1)
    pdb2 = Path(args.pdb2)
    auto_yes = bool(args.yes)

    # Resolve global ssw_test binary if provided via CLI or environment
    ssw_global = None
    if args.ssw:
        ssw_global = Path(args.ssw).resolve()
    elif os.environ.get('SSW_TEST'):
        ssw_global = Path(os.environ.get('SSW_TEST')).resolve()
    if ssw_global:
        if not ssw_global.exists():
            print(f'Error: provided global ssw_test binary not found: {ssw_global}')
            sys.exit(1)
        if not os.access(str(ssw_global), os.X_OK):
            print(f'Error: provided global ssw_test is not executable: {ssw_global}')
            sys.exit(1)
        print(f'Using global ssw_test: {ssw_global}')

    # Determine model variants to run. Option B semantics: if --model-dir is provided,
    # treat it as an override and run only that model; otherwise auto-detect the three
    # standard variants under repository roots.
    variant_model_dirs = []  # list of (suffix, path-to-tmp)
    if args.model_dir:
        provided = Path(args.model_dir)
        # normalize to absolute
        provided_abs = provided.resolve()
        # try to detect a suffix like '8f','9f','10f' in the path
        m = re.search(r'(8f|9f|10f)', str(provided_abs).lower())
        suffix = m.group(1) if m else provided_abs.name
        variant_model_dirs.append((suffix, provided_abs))
    else:
        for v in ['8f', '9f', '10f']:
            cand = Path(f'encoders_and_tools/training_3di_gpu_{v}')
            if cand.exists():
                variant_model_dirs.append((v, cand.resolve()))
        if not variant_model_dirs:
            # fallback to the historical default (8f)
            variant_model_dirs.append(('8f', (Path('encoders_and_tools/training_3di_gpu_8f')).resolve()))

    # set new_suffix to the first variant (used by some downstream file naming logic)
    new_suffix = variant_model_dirs[0][0]

    # compute a single timestamp used for both tmp workspace and outdir
    ts = time.strftime('%Y%m%d%H%M%S')

    # output dir (final results). Default: project root `ROOT/pairwise_3di_{ts}`
    if args.outdir:
        outdir = Path(args.outdir).resolve()
        outdir.mkdir(parents=True, exist_ok=True)
    else:
        outdir = Path(ROOT) / f"pairwise_3di_{ts}"
        outdir.mkdir(parents=True, exist_ok=True)

    logger.info(f'Outputs will be written to {outdir}')

    # Step 1: copy pdbs to tmp dir
    if not confirm('Step 1: copy input PDBs to temporary workspace', auto_yes):
        sys.exit(1)

    # Create timestamped temporary workspace under project tmp directory (not /tmp)
    project_tmp_root = Path(ROOT) / 'tmp'
    project_tmp_root.mkdir(parents=True, exist_ok=True)
    ts = time.strftime('%Y%m%d%H%M%S')
    tmpdir = project_tmp_root / f"pairwise_3di_{ts}"
    # ensure unique by appending a counter if exists
    i = 0
    base_tmpdir = tmpdir
    while tmpdir.exists():
        i += 1
        tmpdir = Path(f"{str(base_tmpdir)}_{i}")
    tmpdir.mkdir(parents=True, exist_ok=False)
    print('Using temporary directory (project tmp):', tmpdir)
    # create a run-local log file and tee stdout/stderr into it so every run's output is preserved
    log_path = tmpdir / 'log.txt'
    logf = open(log_path, 'a')
    class Tee:
        def __init__(self, a, b):
            self.a = a
            self.b = b
        def write(self, s):
            """Write to both streams and flush immediately so interactive terminals see output.

            Flushing helps avoid output being stuck in buffers when the script is run
            non-interactively or when stdout/stderr are redirected.
            """
            try:
                self.a.write(s)
            except Exception:
                pass
            try:
                self.b.write(s)
            except Exception:
                pass
            # best-effort flush to minimize lost/stalled output
            try:
                self.a.flush()
            except Exception:
                pass
            try:
                self.b.flush()
            except Exception:
                pass
        def flush(self):
            try:
                self.a.flush()
            except Exception:
                pass
            try:
                self.b.flush()
            except Exception:
                pass
    # replace stdout/stderr so all prints and subprocess outputs go into log.txt as well
    sys.stdout = Tee(sys.stdout, logf)
    sys.stderr = Tee(sys.stderr, logf)
    p1_tmp = tmpdir / pdb1.name
    p2_tmp = tmpdir / pdb2.name
    shutil.copy2(pdb1, p1_tmp)
    shutil.copy2(pdb2, p2_tmp)
    print('Copied:', p1_tmp, p2_tmp)

    # Step 2: foldseek createdb for each pdb
    if not confirm('Step 2: run foldseek createdb on each PDB', auto_yes):
        sys.exit(1)

    for ptmp in (p1_tmp, p2_tmp):
        prefix = ptmp.stem + '_pdzdb'
        # create a protein-specific subdirectory for createdb outputs
        p_sub = tmpdir / ptmp.stem
        # ensure the directory exists and run createdb from inside it
        p_sub.mkdir(parents=True, exist_ok=True)
        # use basename as out_base so foldseek writes files to the current working dir (p_sub)
        out_base_name = prefix
        cmd = [args.foldseek, 'createdb', str(ptmp), out_base_name]
        rc, out, err = run_cmd(cmd, cwd=str(p_sub), capture=True)
        if rc != 0:
            print(f'foldseek createdb failed for {ptmp} (rc={rc})')
            print('stderr:', err)
            print('stdout:', out)
            sys.exit(1)
        else:
            # verify that createdb produced files inside p_sub
            found = any(f.name.startswith(prefix) for f in p_sub.iterdir())
            if not found:
                print('Warning: foldseek reported success but no createdb files with prefix', prefix, 'found in', p_sub)
                print('Listing parent tmpdir contents for debugging:')
                for x in tmpdir.iterdir():
                    print(' ', x.name)
            print('foldseek createdb done for', ptmp, '->', p_sub)

    # Step 3: generate FASTA and 3di_original from createdb outputs
    if not confirm('Step 3: create FASTA and 3di_original from createdb outputs (strip null param)', auto_yes):
        sys.exit(1)

    def process_createdb_dir(p_sub: Path, prefix: str, orig_name: str):
        """Process createdb outputs in p_sub. Write two files and remove others:
        - {orig_name}.fasta  (cleaned FASTA)
        - {orig_name}_10f_3di.txt  (sequence only)
        """
        fasta_path = p_sub / f"{orig_name}.fasta"
        seq_path = p_sub / f"{orig_name}_10f_3di.txt"

        # Prefer files that exactly match the createdb base names (no startswith)
        pdz = None
        pdz_ss = None

        # decide where to search: prefer p_sub, else its parent tmpdir
        search_dir = p_sub if p_sub.exists() else p_sub.parent

        # exact expected filenames
        candidate_pdz = search_dir / prefix
        candidate_pdz_ss = search_dir / f"{prefix}_ss"
        candidate_pdz_ss_alt = search_dir / f"{prefix}_pdzdb_ss"

        if candidate_pdz.exists():
            pdz = candidate_pdz
        if candidate_pdz_ss.exists():
            pdz_ss = candidate_pdz_ss
        elif candidate_pdz_ss_alt.exists():
            pdz_ss = candidate_pdz_ss_alt

        # fallback: look for a file whose name exactly equals prefix or known ss variants
        if pdz is None:
            for f in search_dir.iterdir():
                if f.name == prefix:
                    pdz = f
                    break
        if pdz_ss is None:
            for f in search_dir.iterdir():
                if f.name in (f"{prefix}_ss", f"{prefix}_pdzdb_ss"):
                    pdz_ss = f
                    break

        # ensure protein folder exists so we place cleaned outputs there
        try:
            p_sub.mkdir(parents=True, exist_ok=True)
        except Exception:
            pass

        # Write cleaned FASTA only from the exact pdz file
        if pdz and pdz.exists():
            text = cleanse_createdb_file(pdz)
            seq_text = extract_sequence_from_file_content(text)
            # write as proper FASTA with header
            write_fasta_from_string(fasta_path, orig_name, seq_text)
            print('Wrote FASTA:', fasta_path)
        else:
            print('pdz file not found for', orig_name, 'in', search_dir)

        # Write 3Di original only from the exact pdz_ss file
        if pdz_ss and pdz_ss.exists():
            text = cleanse_createdb_file(pdz_ss)
            seq = extract_sequence_from_file_content(text)
            seq_path.write_text(seq + '\n')
            print('Wrote 10f_3di seq file:', seq_path)
        else:
            print('pdz_ss (3di) file not found for', orig_name, 'in', search_dir)

        # remove other files in protein folder, keep only fasta_path and seq_path
        for f in list(p_sub.iterdir()):
            try:
                if f.resolve() == fasta_path.resolve() or f.resolve() == seq_path.resolve():
                    continue
                if f.is_file():
                    f.unlink()
                else:
                    shutil.rmtree(f)
            except Exception as e:
                print('Warning: failed to remove', f, e)

        # If the raw createdb files were found in the parent search_dir (not p_sub), remove them to avoid duplicates
        for candidate in [pdz, pdz_ss]:
            try:
                if candidate and candidate.exists() and candidate.parent.resolve() != p_sub.resolve():
                    try:
                        candidate.unlink()
                    except Exception:
                        pass
            except Exception:
                pass

        return fasta_path, seq_path

    # call process_createdb_dir for each protein subdir (created above)
    p1_sub = tmpdir / p1_tmp.stem
    p2_sub = tmpdir / p2_tmp.stem
    prefix1 = p1_tmp.stem + '_pdzdb'
    prefix2 = p2_tmp.stem + '_pdzdb'
    fasta1, orig3di1 = process_createdb_dir(p1_sub, prefix1, p1_tmp.stem)
    fasta2, orig3di2 = process_createdb_dir(p2_sub, prefix2, p2_tmp.stem)

    # Step 4: use pdb_to_3di.py to create 3di_new for both PDBs for each detected variant
    if not confirm('Step 4: create 3di_new using pdb_to_3di.py for each model variant', auto_yes):
        sys.exit(1)

    # new3di_by_var: maps suffix -> (path1, path2)
    new3di_by_var = {}
    for suffix, model_tmp in variant_model_dirs:
        for ptmp in (p1_tmp, p2_tmp):
            p_sub = tmpdir / ptmp.stem
            p_sub.mkdir(parents=True, exist_ok=True)
            out_name = f"{ptmp.stem}_{suffix}_3di.txt"
            out_path = p_sub / out_name
            model_dir_abs = str(model_tmp)
            cmd = [sys.executable, str(Path(__file__).parent / 'pdb_to_3di.py'), str(ptmp), model_dir_abs, '-o', str(out_path), '--no-header']
            rc, out, err = run_cmd(cmd, cwd=str(p_sub), capture=True)
            if rc != 0:
                print(f'pdb_to_3di failed for variant {suffix} on {ptmp} (rc={rc})')
                print('stdout:', out)
                print('stderr:', err)
                # do not abort entire run; skip this variant
                break
            else:
                print(f'pdb_to_3di ({suffix}) output:', out_path)
        # after both produced, verify
        f1 = tmpdir / p1_tmp.stem / f"{p1_tmp.stem}_{suffix}_3di.txt"
        f2 = tmpdir / p2_tmp.stem / f"{p2_tmp.stem}_{suffix}_3di.txt"
        if f1.exists() and f2.exists():
            new3di_by_var[suffix] = (f1, f2)
        else:
            print(f'Warning: new3di files for variant {suffix} not found; skipping this variant')

    if not new3di_by_var:
        print('Error: no variant produced 3di outputs; aborting')
        sys.exit(1)

    # Step 5: run TM-align
    if not confirm('Step 5: run TM-align to align the two PDBs', auto_yes):
        sys.exit(1)

    tmalign_out = tmpdir / 'TM_align_result.txt'
    cmd = [args.tmalign, str(p1_tmp), str(p2_tmp)]
    rc, out, err = run_cmd(cmd, capture=True)
    (tmalign_out).write_text(out + '\n' + err)
    print('TM-align output written to', tmalign_out)

    # Step 6/7: parse TM-align output for aligned ranges
    if not confirm('Step 6: parse TM-align result to get aligned residue ranges', auto_yes):
        sys.exit(1)

    parsed = parse_tmalign_ranges(out + '\n' + err)
    if parsed is None:
        print('Failed to parse TM-align alignment block. The script will continue, but trimming will be skipped.')
        parsed_ranges = None
    else:
        # parsed now contains 0-based inclusive index ranges for aligned region: ((q0,q1),(s0,s1))
        parsed_ranges = parsed
        print('Parsed TM-align index ranges (0-based inclusive):', parsed_ranges)

    # Step 8: Trim sequences based on ranges
    if not confirm('Step 8: trim original and new 3di sequences according to TM-align ranges (if parsed)', auto_yes):
        sys.exit(1)

    def read_sequence_from_3di_file(path: Path) -> str:
        text = path.read_text()
        lines = [l.strip() for l in text.splitlines() if l.strip()]
        if lines and lines[0].startswith('>'):
            return ''.join(lines[1:])
        return ''.join(lines)

    # produce trimmed files named per the variant: e.g. *_8f_trim.txt, *_9f_trim.txt and original remains 10f
    trimmed_txt_new_1 = tmpdir / (p1_tmp.stem + f'_{new_suffix}_trim.txt')
    trimmed_txt_new_2 = tmpdir / (p2_tmp.stem + f'_{new_suffix}_trim.txt')
    trimmed_txt_10f_1 = tmpdir / (p1_tmp.stem + '_10f_trim.txt')
    trimmed_txt_10f_2 = tmpdir / (p2_tmp.stem + '_10f_trim.txt')

    fasta_new_1 = tmpdir / (p1_tmp.stem + f'_{new_suffix}_trim.fasta')
    fasta_new_2 = tmpdir / (p2_tmp.stem + f'_{new_suffix}_trim.fasta')
    fasta_10f_1 = tmpdir / (p1_tmp.stem + '_10f_trim.fasta')
    fasta_10f_2 = tmpdir / (p2_tmp.stem + '_10f_trim.fasta')

    # Initialize fasta-based trimmed paths to None; later we'll point trimmed_new*/trimmed_orig* to FASTA paths
    trimmed_new1 = None
    trimmed_new2 = None
    trimmed_orig1 = None
    trimmed_orig2 = None

    # produce trimmed FASTAs per variant based on parsed TM-align ranges
    trimmed_fastas_by_var = {}  # suffix -> (fasta_a, fasta_b)
    if parsed_ranges:
        (r1, r2) = parsed_ranges
        i0_1, i1_1 = r1
        i0_2, i1_2 = r2

        for suffix, (f1, f2) in new3di_by_var.items():
            seq_new1 = read_sequence_from_3di_file(f1)
            seq_new2 = read_sequence_from_3di_file(f2)
            # validate indices
            if i0_1 < 0 or i1_1 >= len(seq_new1) or i0_2 < 0 or i1_2 >= len(seq_new2):
                print(f'Parsed TM-align indices out of range for variant {suffix}; skipping trimming for this variant')
                continue
            sub_a = seq_new1[i0_1:i1_1+1]
            sub_b = seq_new2[i0_2:i1_2+1]
            fasta_a = tmpdir / (p1_tmp.stem + f'_{suffix}_trim.fasta')
            fasta_b = tmpdir / (p2_tmp.stem + f'_{suffix}_trim.fasta')
            write_fasta_from_string(fasta_a, p1_tmp.stem + f'_{suffix}_trim', sub_a)
            # query doubled (x2) for ssw usage where needed
            write_fasta_from_string(fasta_b, p2_tmp.stem + f'_{suffix}_trim', sub_b + sub_b)
            trimmed_fastas_by_var[suffix] = (fasta_a, fasta_b)
            print('Wrote trimmed FASTAs for variant', suffix, fasta_a, fasta_b)
    else:
        print('Skipping trimming because no parsed TM-align ranges are available.')

    # Step 9/10: run ssw_test for each variant using the variant's model tmp dir
    if not confirm('Step 9/10: run ssw_test for each model variant (8f,9f,10f) using their sub_score.mat', auto_yes):
        sys.exit(1)

    align_dir = outdir / 'alignment_results'
    align_dir.mkdir(parents=True, exist_ok=True)

    for suffix, model_tmp in variant_model_dirs:
        print(f'--- ssw for variant {suffix} ---')
        if suffix not in trimmed_fastas_by_var:
            print('No trimmed FASTAs for variant', suffix, '; skipping ssw for this variant')
            continue
        fasta_a, fasta_b = trimmed_fastas_by_var[suffix]
        # model_root is parent of model_tmp (tmp dir)
        model_root = Path(model_tmp).resolve().parent
        # Do NOT change cwd when invoking ssw; call the binary with absolute paths instead.
        ssw_cwd = None
        # choose ssw binary: global override (args.ssw or env SSW_TEST) or per-variant binary
        if 'ssw_global' in locals() and ssw_global:
            ssw_exec = str(ssw_global)
        else:
            ssw_exec = str(Path(model_tmp) / 'ssw_test')
            if not Path(ssw_exec).exists():
                logger.warning(f'ssw_test not found for variant {suffix} at {ssw_exec}; invocation may fail')

        # matrix path inside model tmp (absolute)
        # Try multiple locations: model_tmp, then encoders_and_tools variant dir
        matrix_path = None
        candidates = [
            Path(model_tmp) / 'sub_score.mat',
            Path(model_tmp) / 's.mat',
            Path(model_tmp) / f's_{suffix}.mat',
        ]
        for cand in candidates:
            if cand.exists():
                matrix_path = cand
                break
        if not matrix_path:
            print(f'Error: scoring matrix not found for variant {suffix} in {model_tmp}.')
            sys.exit(1)

        # ensure we copy the matrix into outdir for record
        try:
            shutil.copy2(matrix_path, outdir / matrix_path.name)
        except Exception:
            pass
        matrix_abs = str(matrix_path.resolve())

        # Validate that the matrix alphabet covers the sequences we'll align
        mat_alpha = read_matrix_alphabet(matrix_path)
        if mat_alpha is None:
            # try normalizing and re-read
            try:
                norm_mat = tmpdir / (matrix_path.stem + '_norm.mat')
                normalize_matrix_file(matrix_path, norm_mat)
                mat_alpha = read_matrix_alphabet(norm_mat)
                if mat_alpha:
                    matrix_abs = str(norm_mat.resolve())
                    try:
                        shutil.copy2(norm_mat, outdir / norm_mat.name)
                    except Exception:
                        pass
            except Exception:
                pass
        if mat_alpha is None:
            print(f'Error: could not parse alphabet from scoring matrix for variant {suffix}: {matrix_path}. Aborting.')
            sys.exit(1)

        seqA = set(read_fasta_sequence(fasta_a).upper())
        seqB = set(read_fasta_sequence(fasta_b).upper())
        seq_chars = set([c for c in seqA.union(seqB) if c and not c.isspace()])
        missing = sorted([c for c in seq_chars if c not in mat_alpha])
        if missing:
            print(f'Error: sequence contains symbols not present in scoring matrix for variant {suffix}: missing={missing}. Matrix alphabet={sorted(mat_alpha)}')
            sys.exit(1)

        # prepare output paths per-variant
        fasta1_name = fasta_a.stem
        fasta2_name = fasta_b.stem
        out_blast = align_dir / f"{fasta1_name}_{fasta2_name}_{suffix}_blast_output.txt"
        out_ssw = align_dir / f"{fasta1_name}_{fasta2_name}_{suffix}_ssw_output.txt"

        # Always invoke ssw_test with a fixed set of options that produce SAM/CIGAR and use the variant matrix
        cmd_unified = [ssw_exec, '-p', '-c', '-o', str(args.gap_open), '-e', str(args.gap_extend), '-a', matrix_abs, str(fasta_a), str(fasta_b)]
        rc_u, out_u, err_u = run_cmd(cmd_unified, cwd=ssw_cwd, capture=True)

        def try_convert_ssw_to_blast(ssw_out_path: Path):
            sam2blast = Path(__file__).parent / 'scripts' / 'sam_to_blastlike.py'
            if sam2blast.exists():
                try:
                    cmd_conv = [sys.executable, str(sam2blast), str(fasta_a), str(ssw_out_path)]
                    rc_c, out_c, err_c = run_cmd(cmd_conv, cwd=ROOT, capture=True)
                    if rc_c == 0:
                        out_blast.write_text(out_c)
                        logger.info(f'Wrote BLAST-like (from SAM) for {suffix} {out_blast}')
                        return True
                    else:
                        logger.warning(f'sam_to_blastlike failed for {suffix} rc={rc_c} err={err_c}')
                except Exception:
                    logger.exception('Failed to invoke sam_to_blastlike')
            else:
                logger.debug('sam_to_blastlike.py not found; skipping conversion')
            return False

        if rc_u == 0:
            out_ssw.write_text(out_u)
            logger.info(f'Wrote ssw (-c) output for {suffix} {out_ssw}')
            try_convert_ssw_to_blast(out_ssw)
        else:
            logger.warning(f'ssw (-c) run failed for {suffix} rc={rc_u} stderr={err_u}')

            # Attempt a simple retry without the -a matrix override (use default matrix)
            if matrix_abs:
                cmd_retry = [ssw_exec, '-p', '-c', '-o', str(args.gap_open), '-e', str(args.gap_extend), str(fasta_a), str(fasta_b)]
                rc_r, out_r, err_r = run_cmd(cmd_retry, cwd=ssw_cwd, capture=True)
                if rc_r == 0:
                    out_ssw.write_text(out_r)
                    logger.info(f'Wrote ssw (-c) output for {suffix} {out_ssw} (fallback without -a)')
                    try_convert_ssw_to_blast(out_ssw)
                    continue

            # Try normalizing the provided matrix and retry with the normalized file
            try:
                if matrix_path.exists():
                    norm_mat = tmpdir / (matrix_path.stem + '_norm.mat')
                    normalize_matrix_file(matrix_path, norm_mat)
                    cmd_norm = [ssw_exec, '-p', '-c', '-a', str(norm_mat), '-o', str(args.gap_open), '-e', str(args.gap_extend), str(fasta_a), str(fasta_b)]
                    rc_n, out_n, err_n = run_cmd(cmd_norm, cwd=ssw_cwd, capture=True)
                    if rc_n == 0:
                        out_ssw.write_text(out_n)
                        logger.info(f'Wrote ssw (-c) output for {suffix} {out_ssw} (normalized matrix)')
                        try_convert_ssw_to_blast(out_ssw)
                        continue
                    else:
                        logger.error(f'ssw (-c) run failed for {suffix} even with normalized matrix (rc={rc_n}) stderr={err_n}')
            except Exception:
                logger.exception('Matrix normalization attempt failed')

            # As a last resort, try the converter script to produce a converted matrix and retry
            conv_script = Path(__file__).parent / 'scripts' / 'convert_smat_for_ssw.py'
            try:
                if conv_script.exists() and matrix_path.exists():
                    cmd_conv = [sys.executable, str(conv_script), '-i', str(matrix_path), '--pairwise-dir', str(model_tmp)]
                    rc_conv, out_conv_text, err_conv = run_cmd(cmd_conv, cwd=ROOT, capture=True)
                    logger.info(f'Converter rc={rc_conv} out={out_conv_text} err={err_conv}')
                    if rc_conv == 0:
                        converted_candidates = list(Path(model_tmp).glob('converted*.mat'))
                        if converted_candidates:
                            conv_mat = converted_candidates[0]
                            rel_conv = os.path.join('tmp', conv_mat.name)
                            cmd_conv_retry = [ssw_exec, '-p', '-c', '-a', rel_conv, '-o', str(args.gap_open), '-e', str(args.gap_extend), str(fasta_a), str(fasta_b)]
                            rc4, out4, err4 = run_cmd(cmd_conv_retry, cwd=ssw_cwd, capture=True)
                            if rc4 == 0:
                                out_ssw.write_text(out4)
                                logger.info(f'Wrote ssw (-c) output for {suffix} {out_ssw} (converted matrix)')
                                try_convert_ssw_to_blast(out_ssw)
                                continue
                            else:
                                logger.error(f'ssw (-c) run failed for {suffix} even after convert_smat_for_ssw (rc={rc4}) stderr={err4}')
            except Exception:
                logger.exception('Failed running convert_smat_for_ssw')

    # Step 11: display per-variant ssw_test outputs and compare with TM-align ranges
    print('\nStep 11: ssw_test outputs (per-variant)')
    if parsed_ranges:
        print('TM-align parsed ranges:', parsed_ranges)
    else:
        print('No TM-align ranges parsed.')

    for suffix in trimmed_fastas_by_var.keys():
        # find files we wrote earlier
        pattern_blast = f"{p1_tmp.stem}_{p2_tmp.stem}_{suffix}_blast_output.txt"
        # but actual filenames use fasta stems; relax search in alignment_results
        print('\nVariant', suffix, 'alignment results:')
        found = False
        for f in (outdir / 'alignment_results').iterdir():
            if f.is_file() and f.name.endswith(f"_{suffix}_blast_output.txt"):
                found = True
                print('\n--- BLAST-like (first 80 chars) ---')
                txt = f.read_text()
                print('\n'.join(txt.splitlines()[:20]))
                # try parse
                try:
                    entries = parse_ssw_results(txt)
                    print('\nParsed top entries:')
                    for e in entries[:5]:
                        print(e)
                except Exception:
                    pass
        if not found:
            print(' No BLAST-like file found for variant', suffix)

    # Step 13: gather and copy intermediate outputs to outdir
    print('\nStep 13: collecting intermediate outputs to outdir')
    for p in tmpdir.iterdir():
        try:
            # copy relevant files (small)
            if p.is_file():
                shutil.copy2(p, outdir / p.name)
        except Exception:
            pass
    print('All temp files copied to', outdir)
    print('Temporary workspace retained at', tmpdir)
    print('Pipeline finished.')


if __name__ == '__main__':
    # Execute the pipeline via main()
    try:
        main()
    except Exception as e:
        # Ensure any uncaught exception is logged and re-raised so the caller sees it
        try:
            logger.exception('Unhandled exception in main: %s', e)
        except Exception:
            print('Unhandled exception in main:', e)
        raise
