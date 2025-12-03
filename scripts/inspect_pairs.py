#!/usr/bin/env python3
"""Standalone inspect_pairs utility that reuses the pipeline's parsing and
index-correction logic so results match the pipeline's comparisons.

Usage: scripts/inspect_pairs.py [--rundir <run_dir>] [--verbose]
If no rundir is given, the most recent /tmp/tmalign_3di_* run is used.
"""

from pathlib import Path
import argparse
import logging
import sys
from typing import List, Tuple

# import core helpers from the main pipeline to ensure identical behavior
from tmalign_3di_match_pipeline import (
    parse_ssw_sam_pairs,
    unwrap_bprime_to_b,
    slice_to_b,
    rotate_j,
    choose_cutpoint_from_tm,
    read_fasta_sequence,
)

logger = logging.getLogger('inspect_pairs')


def read_tsv_pairs(p: Path) -> List[Tuple[int,int]]:
    pairs = []
    for ln in p.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        toks = ln.split()
        if len(toks) < 2:
            continue
        try:
            a = int(toks[0]); b = int(toks[1])
            pairs.append((a,b))
        except Exception:
            continue
    return pairs


def compare_pairs(tm_pairs: List[Tuple[int,int]], three_pairs: List[Tuple[int,int]], w: int = 1):
    tm_matched = [False] * len(tm_pairs)
    three_matched = [False] * len(three_pairs)
    matched = []
    for j_idx, (i3, j3) in enumerate(three_pairs):
        for t_idx, (it, jt) in enumerate(tm_pairs):
            if tm_matched[t_idx]:
                continue
            if abs(it - i3) <= w and abs(jt - j3) <= w:
                tm_matched[t_idx] = True
                three_matched[j_idx] = True
                matched.append((t_idx, j_idx, (it, jt), (i3, j3)))
                break
    TP = sum(1 for x in tm_matched if x)
    matched_three = sum(1 for x in three_matched if x)
    FN = max(0, len(tm_pairs) - TP)
    FP = max(0, len(three_pairs) - matched_three)
    return {'TP': TP, 'FP': FP, 'FN': FN, 'matched': matched, 'n_tm': len(tm_pairs), 'n_3di': len(three_pairs)}


def detect_and_apply_global_j_offset(pairs_3di_B: List[Tuple[int,int]], tm_pairs: List[Tuple[int,int]], lenT: int, logger) -> List[Tuple[int,int]]:
    try:
        tm_map = {i: j for i, j in tm_pairs}
        deltas = []
        for (ia, jb) in pairs_3di_B:
            if ia in tm_map:
                deltas.append(tm_map[ia] - jb)
        if deltas:
            from collections import Counter
            norm = []
            for d in deltas:
                if lenT:
                    dmod = ((d + (lenT//2)) % lenT) - (lenT//2)
                else:
                    dmod = d
                norm.append(dmod)
            cnt = Counter(norm)
            mode_delta, cnt_mode = cnt.most_common(1)[0]
            if cnt_mode >= max(2, int(0.5 * len(norm))):
                logger.info(f'Detected global SSW->TM j-offset {mode_delta}; applying correction')
                def shift_j(jb, delta):
                    return ((jb + delta - 1) % lenT) + 1 if lenT else jb
                return [(ia, shift_j(jb, mode_delta)) for (ia, jb) in pairs_3di_B]
    except Exception:
        logger.debug('SSW->TM offset detection failed', exc_info=True)
    return pairs_3di_B


def find_rotation_offset(original: str, printed: str) -> int:
    """Return 0-based offset where printed appears within doubled original or -1.

    Matches the pipeline's find_rotation_offset implementation.
    """
    if not printed or not original:
        return -1
    doubled = original + original
    pos = doubled.find(printed)
    if pos != -1 and pos < len(original):
        return pos
    return -1


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--rundir', '-r', help='Path to tmalign tmp run dir (defaults to latest under tmp/)', default=None)
    p.add_argument('--w', type=int, default=1, help='pairing tolerance window (default 1)')
    p.add_argument('--verbose', '-v', action='store_true')
    args = p.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

    ROOT = Path.cwd()
    base = ROOT / 'tmp'
    if args.rundir:
        run = Path(args.rundir)
    else:
        runs = sorted([d for d in base.glob('tmalign_3di_*') if d.is_dir()])
        if not runs:
            logger.error('No runs found under %s', base)
            sys.exit(1)
        run = runs[-1]

    logger.info('Using run dir: %s', run)

    tm_pairs_path = run / 'tm_pairs.tsv'
    if not tm_pairs_path.exists():
        logger.error('tm_pairs.tsv not found at %s', tm_pairs_path)
        sys.exit(1)

    tm_pairs = read_tsv_pairs(tm_pairs_path)
    logger.info('TM pairs count= %d', len(tm_pairs))
    logger.debug('TM pairs (first 50): %s', tm_pairs[:50])

    # collect fasta lengths to assist with query_dup detection
    fasta_lens = {}
    for f in run.rglob('*.fasta'):
        # key by stem
        try:
            seq = ''.join([ln.strip() for ln in f.read_text().splitlines() if ln and not ln.startswith('>')])
            fasta_lens[f.stem] = len(seq)
        except Exception:
            continue

    # set up an optional file log in the run directory so the full process can be inspected later
    try:
        fh = logging.FileHandler(run / 'inspect_pairs_verbose.log')
        fh.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))
        logging.getLogger().addHandler(fh)
        logger.info('Detailed inspect log will be written to %s', str(run / 'inspect_pairs_verbose.log'))
    except Exception:
        logger.debug('Could not create run logfile', exc_info=True)

    # attempt to read TM-align printed snippet so we can compute rotation offsets
    tm_snip = ''
    try:
        snipf = run / 'TM_align_cp_result_snippet.txt'
        if snipf.exists():
            tm_snip = snipf.read_text()
            # attempt to extract two printed sequences (top and bottom) from the snippet
            printed_lines = [ln.strip() for ln in tm_snip.splitlines() if ln.strip()]
            # heuristic: find the first block of three lines where middle is marker line (contains ':' or '.' or '*')
            pt_top = ''
            pt_bot = ''
            for i in range(len(printed_lines)-2):
                a = printed_lines[i]
                b = printed_lines[i+1]
                c = printed_lines[i+2]
                # marker line typically contains ':' or '.' characters and spaces
                if (set(b) & set(':.*')) and any(ch.isalpha() or ch=='-' or ch=='*' for ch in a) and any(ch.isalpha() or ch=='-' or ch=='*' for ch in c):
                    pt_top = a
                    pt_bot = c
                    break
            if pt_top and pt_bot:
                logger.debug('Extracted printed sequences from TM snippet (len top=%d, bot=%d)', len(pt_top), len(pt_bot))
            else:
                logger.debug('Could not heuristically extract printed sequences from TM snippet')
        else:
            logger.debug('TM_align snippet not found at %s', snipf)
    except Exception:
        logger.debug('Error while reading TM_align snippet', exc_info=True)

    # try to parse tm_pairs header to get pdb stems (tm_pairs.tsv contains #pdb1 and #pdb2 lines)
    p1_stem = None
    p2_stem = None
    try:
        for ln in tm_pairs_path.read_text().splitlines():
            if ln.startswith('#pdb1'):
                p1_stem = Path(ln.split('\t',1)[1].strip()).stem
            elif ln.startswith('#pdb2'):
                p2_stem = Path(ln.split('\t',1)[1].strip()).stem
            if p1_stem and p2_stem:
                break
    except Exception:
        logger.debug('Failed to parse pdb stems from tm_pairs header', exc_info=True)

    seqA = ''
    seqB = ''
    try:
        if p1_stem:
            # find fasta under run dir that contains p1_stem
            cand = list(run.rglob(f"*{p1_stem}*.fasta"))
            if cand:
                seqA = read_fasta_sequence(cand[0])
        if p2_stem:
            cand = list(run.rglob(f"*{p2_stem}*.fasta"))
            if cand:
                seqB = read_fasta_sequence(cand[0])
    except Exception:
        logger.debug('Failed to read fasta sequences for printed-offset detection', exc_info=True)

    # compute rotation offsets using pipeline's logic if we have printed sequences and originals
    offA = find_rotation_offset(seqA, pt_top) if seqA and 'pt_top' in locals() and pt_top else -1
    offB = find_rotation_offset(seqB, pt_bot) if seqB and 'pt_bot' in locals() and pt_bot else -1
    if offA >= 0 or offB >= 0:
        logger.info(f'Detected rotation offsets from TM snippet: A={offA}, B={offB} (0-based)')

    # find SSW output files (prefer raw ssw outputs)
    # choose cutpoint s using pipeline logic (choose_cutpoint_from_tm)
    lenT = max((j for _, j in tm_pairs), default=0)
    try:
        s = choose_cutpoint_from_tm(tm_pairs, 'minj', lenT)
    except Exception:
        s = min((j for _, j in tm_pairs), default=1)
    logger.info(f'Chosen cutpoint s={s} len_target={lenT}')

    ssw_files = list(run.glob('*_ssw_output.txt')) + list((run / 'alignment_results').glob('*_ssw_output.txt'))
    ssw_files = [f for f in ssw_files if f.is_file()]
    if not ssw_files:
        logger.error('No *_ssw_output.txt files found in run dir')
        sys.exit(1)

    logger.info('Found %d ssw output files', len(ssw_files))

    for f in ssw_files:
        logger.info('Parsing SSW file: %s', f)
        sam_text = f.read_text()

        # try to infer query_duped_len: if a fasta stem with '_x2' appears in filename, use half of its length
        q_duped_len = None
        for stem, L in fasta_lens.items():
            if stem in f.name and stem.endswith('_x2'):
                q_duped_len = L // 2
                logger.debug('Inferred query_duped_len=%d from fasta %s', q_duped_len, stem)
                break
        # fallback: derive from first SAM seq field
        if q_duped_len is None:
            for ln in sam_text.splitlines():
                if ln.strip() and not ln.startswith('@'):
                    flds = ln.split('\t')
                    if len(flds) >= 10:
                        seq = flds[9]
                        if seq:
                            q_duped_len = len(seq) // 2
                    break

        pairs_3di = parse_ssw_sam_pairs(sam_text, query_duped_len=q_duped_len, debug=args.verbose)
        logger.info('  %s: parsed %d 3di aligned pairs (after mod mapping)', f.name, len(pairs_3di))

    # Map target indices back to B coordinates using unwrap or slice heuristics if needed
    # lenT and cutpoint s were computed earlier from TM pairs
        pairs_3di_B = []
        max_jp = max((j for _, j in pairs_3di), default=0)
        for (ia, jp) in pairs_3di:
            # decide mapping: if jp > lenT and lenT>0, assume using B' coordinates
            if max_jp > lenT and lenT>0:
                jb = unwrap_bprime_to_b(jp, lenT)
            else:
                jb = jp
            pairs_3di_B.append((ia, jb))

        # detect and apply global j offset similar to pipeline
        pairs_3di_B = detect_and_apply_global_j_offset(pairs_3di_B, tm_pairs, lenT, logger)

    # cutpoint s is pre-computed; reuse it for consistency with the pipeline
        pairs_3di_rot = [(ia, rotate_j(jb, s, lenT)) for (ia, jb) in pairs_3di_B]
        pairs_3di_rot_sorted = sorted(pairs_3di_rot)
        tm_rot_sorted = sorted([(ia, rotate_j(jb, s, lenT)) for (ia, jb) in tm_pairs])

        metrics = compare_pairs(tm_rot_sorted, pairs_3di_rot_sorted, w=args.w)

        print('\nFile:', f.name)
        print('  TM count=%d, SSW count=%d' % (metrics['n_tm'], metrics['n_3di']))
        print('  TP=%d, FP=%d, FN=%d' % (metrics['TP'], metrics['FP'], metrics['FN']))
        if metrics['matched']:
            print('  sample matches:')
            for m in metrics['matched'][:10]:
                print('   ', m)
        else:
            print('  No matched pairs.')
        print('  first TM pairs:', tm_pairs[:10])
        print('  first SSW pairs (raw parsed):', pairs_3di[:10])
        print('  first SSW pairs (mapped B):', pairs_3di_B[:10])
        print('  first SSW pairs (rotated):', pairs_3di_rot_sorted[:10])


if __name__ == '__main__':
    main()
