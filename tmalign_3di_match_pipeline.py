#!/usr/bin/env python3
"""
tmalign_3di_match_pipeline.py

New pipeline that:
 - Runs TMalign with -cp (circular permutation option) and writes raw output to TM_align_cp_result.txt
 - Does NOT parse TM-align output for trimming
 - Locates 3Di SSW output files (SAM-like) produced by the pairwise pipeline
 - Extracts aligned residue pairs from the SSW SAM CIGAR (handles '=', 'X', 'M', 'I', 'D')
 - For doubled-query sequences (SSW query doubled), maps query indices > lenB back via modulo
 - Compares TM-align residue pairs (reference) to 3Di pairs using tolerance window w
 - Computes precision, recall (sensitivity), F1 and Jaccard indices and writes a CSV summary

Usage examples:
  python tmalign_3di_match_pipeline.py /abs/path/A.pdb /abs/path/B.pdb \
      --alignment-results-dir /abs/path/pairwise_3di_<ts>/alignment_results --w 2 --tmalign /usr/local/bin/TMalign

        row = {
 - This script expects TMalign to be accessible at the path given by --tmalign (or in PATH).
 - If --skip-tmalign is provided, the script will attempt to load an existing TM_align_cp_result.txt from the temp dir or from --tmalign-output-file.
 - The mapping extraction from TMalign output is heuristic: it extracts lines that look like a residue-pair list (two integers per line). Inspect the saved TM_align_cp_result.txt if needed.
"""

from pathlib import Path
import argparse
import csv
import subprocess
import sys
import re
import os
import shutil
import time
import logging
from typing import List, Tuple, Optional

# basic project/root setup and early logging
ROOT = os.path.abspath(os.path.dirname(__file__))
_run_logs_dir = Path(ROOT) / 'run_logs'
try:
    _run_logs_dir.mkdir(parents=True, exist_ok=True)
except Exception:
    pass
_ts0 = time.strftime('%Y%m%d%H%M%S')
_start_log = _run_logs_dir / f'tmalign_3di_match_start_{_ts0}.log'
logger = logging.getLogger('tmalign_3di_match')
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
        pass


def write_fasta_from_string(target_path: Path, header: str, seq: str):
    target_path.write_text(f">{header}\n{seq}\n")


def write_3di_fasta(out_path: Path, header: str, token_seq: str):
    out_path.write_text(f">{header}\n{token_seq}\n")


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


def read_fasta_sequence(path: Path) -> str:
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



def run_cmd(cmd: List[str], capture: bool = True, cwd: Optional[str] = None) -> Tuple[int, str, str]:
    print('> ' + ' '.join(cmd))
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE if capture else None,
                              stderr=subprocess.PIPE if capture else None, check=False, text=True, cwd=cwd)
        out = proc.stdout if capture else ''
        err = proc.stderr if capture else ''
        return proc.returncode, out, err
    except FileNotFoundError as e:
        return 127, '', str(e)


# ------------------ TM-align parsing ------------------
def parse_tmalign_cp_pairs(text: str) -> List[Tuple[int,int]]:
    """Heuristic parse of TMalign -cp output to extract residue pair mappings.

    We look for lines containing two integers (optionally more columns) and return
    the first two integers as a pair. We skip lines that are clearly header/summary
    (containing letters like 'TM-score') by requiring the line have at least two integer tokens.
    """
    pairs = []
    # Heuristics: look for patterns like '123 456', '123->456', '123:456', or '(123,456)'
    for ln in text.splitlines():
        ln = ln.strip()
        if not ln:
            continue
        # prefer explicit arrow/colon patterns first
        m = re.findall(r'(\d+)\s*[-:>]+\s*(\d+)', ln)
        if m:
            for a,b in m:
                try:
                    pairs.append((int(a), int(b)))
                except Exception:
                    pass
            continue
        # look for parenthesized pairs
        m2 = re.findall(r'\(\s*(\d+)\s*,\s*(\d+)\s*\)', ln)
        if m2:
            for a,b in m2:
                try:
                    pairs.append((int(a), int(b)))
                except Exception:
                    pass
            continue
        # fallback: if line contains only digits, spaces, commas, parentheses or separators, extract first two ints
        if re.search(r'[A-Za-z]', ln):
            # skip typical header lines containing text
            continue
        toks = re.findall(r'\d+', ln)
        if len(toks) >= 2:
            try:
                a = int(toks[0]); b = int(toks[1])
                pairs.append((a,b))
            except Exception:
                pass
    return pairs


def parse_tmalign_visual_pairs(text: str, debug: bool = False) -> Tuple[List[Tuple[int,int]], List[Tuple[int,int,str]], str, str]:
    """Parse TM-align visual 3-line alignment blocks (top, mid, bot) to recover (i,j) pairs.

    Behaviour:
    - Locate visual mid-lines (those containing ':' or '.') and choose the last such region.
      Start scanning from 4 lines above the last mid-line (or from 0 if fewer lines exist).
    - Scan non-overlapping 3-line visual blocks from that start position to the end.
    - Maintain global residue counters (i for top, j for bot) across blocks. For each column
      where both top and bot present (not '-' or space) and mid is ':' or '.', emit a pair.
    - Collect concatenated printed top/bot sequences (including any '*' markers). Return them
      so the caller can detect circular permutation markers and map printed indices back to
      original fasta indices.

    Returns (pairs, pairs_labeled, printed_top_seq, printed_bot_seq)
    """
    lines = [ln.rstrip('\n') for ln in text.splitlines()]
    pairs: List[Tuple[int,int]] = []
    pairs_labeled: List[Tuple[int,int,str]] = []
    top_concat: List[str] = []
    bot_concat: List[str] = []

    L = len(lines)
    if L < 3:
        return pairs, pairs_labeled, '', ''

    # find all mid-line indices that look like visual alignment mid-lines
    mid_idxs = []
    for ii in range(0, L - 2):
        t = lines[ii].rstrip()
        m = lines[ii + 1].rstrip()
        b = lines[ii + 2].rstrip()
        if (re.search(r'[A-Za-z\*\-]', t) and re.search(r'[A-Za-z\*\-]', b)
                and re.search(r'[:\.]', m) and len(t) >= 5 and len(b) >= 5):
            mid_idxs.append(ii + 1)  # index of the middle line

    if not mid_idxs:
        return pairs, pairs_labeled, '', ''

    # start scanning from 4 lines above the last mid-line (user requested behaviour)
    last_mid = max(mid_idxs)
    start_idx = max(0, last_mid - 4)
    if debug:
        logger.info(f'Found {len(mid_idxs)} visual mid-line(s); starting parse at line {start_idx+1} (4 lines above last mid-line {last_mid+1})')

    idx = start_idx
    global_i = 0
    global_j = 0
    # track star positions (printed index among concatenated printed sequence, 1-based)
    top_star_pos: Optional[int] = None
    bot_star_pos: Optional[int] = None

    while idx < L - 2:
        top = lines[idx].rstrip()
        mid = lines[idx + 1].rstrip()
        bot = lines[idx + 2].rstrip()
        if not (re.search(r'[A-Za-z\*\-]', top) and re.search(r'[A-Za-z\*\-]', bot)
                and re.search(r'[:\.]', mid) and len(top) >= 5 and len(bot) >= 5):
            idx += 1
            continue

        # Normalize widths and scan columns
        width = max(len(top), len(mid), len(bot))
        top_s = top.ljust(width)
        mid_s = mid.ljust(width)
        bot_s = bot.ljust(width)

        if debug:
            logger.info('Found visual block at lines %d-%d' % (idx + 1, idx + 3))
            logger.info('TOP : %s' % top_s)
            logger.info('MID : %s' % mid_s)
            logger.info('BOT : %s' % bot_s)
            logger.info('Starting global_i=%d global_j=%d' % (global_i, global_j))

        for col_idx, (c1, cm, c2) in enumerate(zip(top_s, mid_s, bot_s), start=1):
            inc1 = (c1 != '-' and c1 != ' ')
            inc2 = (c2 != '-' and c2 != ' ')
            before_i = global_i
            before_j = global_j
            if inc1:
                global_i += 1
                # record star position if present at this printed residue
                if c1 == '*':
                    top_star_pos = len([c for c in top_concat if c and c != ' ' and c != '-']) + 1
                top_concat.append(c1)
            if inc2:
                global_j += 1
                if c2 == '*':
                    bot_star_pos = len([c for c in bot_concat if c and c != ' ' and c != '-']) + 1
                bot_concat.append(c2)
            if debug:
                logger.info(f' COL {col_idx:03d}: top="{c1}" mid="{cm}" bot="{c2}" -> inc1={inc1} inc2={inc2} i:{before_i}->{global_i} j:{before_j}->{global_j}')
            if inc1 and inc2 and cm in (':', '.'):
                label = 'close' if cm == ':' else 'aligned'
                pairs.append((global_i, global_j))
                pairs_labeled.append((global_i, global_j, label))
                if debug:
                    logger.info(f'  -> ADDED pair ({global_i},{global_j}) label={label}')

        # advance past this 3-line block
        idx += 3

    # deduplicate while preserving order
    seen = set()
    uniq_pairs: List[Tuple[int,int]] = []
    uniq_pairs_labeled: List[Tuple[int,int,str]] = []
    for ii, p in enumerate(pairs):
        if p in seen:
            continue
        seen.add(p)
        uniq_pairs.append(p)
        uniq_pairs_labeled.append(pairs_labeled[ii])

    top_seq = ''.join([c for c in top_concat if c and c != ' ' and c != '-'])
    bot_seq = ''.join([c for c in bot_concat if c and c != ' ' and c != '-'])
    if debug:
        logger.info(f'Parsed {len(uniq_pairs)} unique pairs; top_seq len={len(top_seq)} bot_seq len={len(bot_seq)} top_star={top_star_pos} bot_star={bot_star_pos}')
    return uniq_pairs, uniq_pairs_labeled, top_seq, bot_seq


# ------------------ SSW SAM parsing -> aligned pairs ------------------

def parse_ssw_sam_pairs(sam_text: str, query_duped_len: Optional[int] = None, debug: bool = False) -> List[Tuple[int,int]]:
    """Parse a SAM-like SSW output and return list of aligned residue pairs (1-based).

    IMPORTANT: SSW SAM lines use reference (POS) then query coordinates in the alignment
    building logic. For compatibility with the rest of this script (which uses pairs as
    (query_index, target_index) == (i, j) from TM-align), this function returns pairs
    in (query, ref) order — i.e. (i, j).

    If debug=True, emit step-by-step logs showing how ref_pos and query_pos change while
    walking the CIGAR string and when pairs are appended.

    For each alignment line (skip headers starting with @), parse fields and build
    the alignment by walking the CIGAR string. We include '=' and 'X' and 'M' as aligned
    columns. Insertions (I) advance query position only. Deletions (D) advance reference position only.

    If query_duped_len is provided (the original non-doubled length), any query index > query_duped_len
    will be mapped back via ((idx-1) % query_duped_len) + 1.
    """
    pairs: List[Tuple[int,int]] = []
    for ln in sam_text.splitlines():
        ln = ln.strip()
        if not ln or ln.startswith('@'):
            continue
        fields = ln.split('\t')
        if len(fields) < 11:
            continue
        # SAM columns: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [OPT]
        try:
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])  # 1-based leftmost ref position
            cigar = fields[5]
            seq = fields[9]
        except Exception:
            continue

        if debug:
            logger.info(f"Parsing SAM entry qname={qname} flag={flag} rname={rname} pos={pos} cigar={cigar} seq_len={len(seq)} query_duped_len={query_duped_len}")

        # parse CIGAR
        ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        ref_pos = pos  # 1-based
        query_pos = 1  # assume query coordinates start at 1
        for num_s, op in ops:
            n = int(num_s)
            if debug:
                logger.info(f" CIGAR op {num_s}{op}: before ref_pos={ref_pos} query_pos={query_pos}")
            if op in ('=', 'X', 'M'):
                # aligned columns: consume both
                for k in range(n):
                    ri = ref_pos + k
                    qi = query_pos + k
                    if query_duped_len:
                        # map back into original range
                        qi_mapped = ((qi - 1) % query_duped_len) + 1
                    else:
                        qi_mapped = qi
                    # append as (query_index, ref_index) to match TM-align (i,j) ordering
                    pairs.append((qi_mapped, ri))
                    if debug:
                        logger.info(f"  -> aligned col k={k}: ref={ri} query(raw)={qi} query(mapped)={qi_mapped}")
                ref_pos += n
                query_pos += n
            elif op == 'I':
                # insertion to reference: consumes query only
                if debug:
                    logger.info(f"  -> insertion: advance query_pos by {n} ({query_pos}->{query_pos + n})")
                query_pos += n
            elif op == 'D':
                # deletion from reference: consumes ref only
                if debug:
                    logger.info(f"  -> deletion: advance ref_pos by {n} ({ref_pos}->{ref_pos + n})")
                ref_pos += n
            else:
                # S,H,P,N - clipping/padding/skip: do not count in aligned pairs
                if op in ('S', 'H', 'P', 'N'):
                    if op == 'S' or op == 'H':
                        # clipped sequences consume query (S) or not (H) depending; conservative: advance query for S
                        if op == 'S':
                            if debug:
                                logger.info(f"  -> soft-clip: advance query_pos by {n} ({query_pos}->{query_pos + n})")
                            query_pos += n
                    continue
        # only process first alignment line per SAM entry (typical for SSW output)
        # break after parsing first alignment
        break
    return pairs


# ------------------ comparison / metrics ------------------

# ------------------ CP / rotation helpers ------------------
def pos_mod(a: int, m: int) -> int:
    r = a % m
    return r + m if r < 0 else r


def rotate_j(j: int, s: int, lenB: int) -> int:
    # 1-based rotation to put cutpoint s at position 1
    return pos_mod(j - s, lenB) + 1


def unwrap_bprime_to_b(j_prime: int, lenB: int) -> int:
    # map B' index back to original B (1..lenB)
    return pos_mod(j_prime - 1, lenB) + 1


def slice_to_b(j_slice: int, t: int, lenB: int) -> int:
    # slice B[t..] + B[1..t-1] back to original B
    return pos_mod(t + j_slice - 2, lenB) + 1


def choose_cutpoint_from_tm(pairs: List[Tuple[int,int]], mode: str, lenB: int) -> int:
    # mode: "minj" | "lcr" | int-as-string
    if not pairs:
        return 1
    if mode.lower() == "minj":
        return min(j for _, j in pairs)
    if mode.lower() == "lcr":
        # longest monotonic run in TM j (Δj>=0 with at most one wrap)
        best_s, best_len = 1, 0
        js = [j for _, j in pairs]
        # detect start of the longest non-decreasing run (split at wrap)
        start = 0
        for k in range(1, len(js)+1):
            if k == len(js) or js[k] < js[k-1]:
                # run [start..k-1]
                run_start_j = js[start]
                run_len = k - start
                if run_len > best_len:
                    best_len = run_len
                    best_s = run_start_j
                start = k
        return best_s
    # integer string
    try:
        s = int(mode)
        s = 1 if s < 1 else (lenB if s > lenB else s)
        return s
    except Exception:
        return 1


def infer_ranges_from_pairs(pairs: List[Tuple[int,int]]) -> Tuple[Tuple[int,int], Tuple[int,int]]:
    # Returns (target_range_minmax, query_range_minmax) in original coordinates (no rotation)
    if not pairs:
        return (0,0), (0,0)
    is_, js_ = zip(*pairs)
    return (min(js_), max(js_)), (min(is_), max(is_))


def circular_overlap_len(a: Tuple[int,int], bs: List[Tuple[int,int]], lenB: int) -> int:
    # a = [l..r] with l<=r (assume caller normalized)
    # bs can include segments that already account for wrapping (caller splits TM wrap into two)
    l,r = a
    total = 0
    for bl, br in bs:
        lo = max(l, bl)
        hi = min(r, br)
        if hi >= lo:
            total += (hi - lo + 1)
    return total


def split_wrap_interval(seg: Tuple[int,int], lenB: int) -> List[Tuple[int,int]]:
    l, r = seg
    if l <= r:
        return [(l, r)]
    # wrapped: [l..lenB] ∪ [1..r]
    return [(l, lenB), (1, r)]


def range_f1(precision: float, recall: float) -> float:
    return (2*precision*recall/(precision+recall)) if (precision+recall)>0 else 0.0


def lctr_longest_consistent_run(pairs_sorted: List[Tuple[int,int]]) -> int:
    # longest run with Δi=Δj=+1
    if not pairs_sorted:
        return 0
    best=1; cur=1
    for k in range(1, len(pairs_sorted)):
        if pairs_sorted[k][0]==pairs_sorted[k-1][0]+1 and pairs_sorted[k][1]==pairs_sorted[k-1][1]+1:
            cur+=1
            best=max(best,cur)
        else:
            cur=1
    return best


def block_recall_tau(tm_pairs: List[Tuple[int,int]], pred_pairs: List[Tuple[int,int]], tau: int) -> float:
    # fraction of TM blocks (Δi=Δj=+1, length≥tau) that are reproduced by any contiguous run in pred of length≥tau
    # Build TM blocks
    blocks=[]
    if tm_pairs:
        s=0
        for k in range(1,len(tm_pairs)+1):
            if k==len(tm_pairs) or not (tm_pairs[k][0]==tm_pairs[k-1][0]+1 and tm_pairs[k][1]==tm_pairs[k-1][1]+1):
                if k-s>=tau:
                    blocks.append((tm_pairs[s][0], tm_pairs[s][1], k-s))
                s=k
    if not blocks:
        return 1.0  # nothing to recall
    # Build lengths of all contiguous runs in pred
    runs=[]
    if pred_pairs:
        s=0
        for k in range(1,len(pred_pairs)+1):
            if k==len(pred_pairs) or not (pred_pairs[k][0]==pred_pairs[k-1][0]+1 and pred_pairs[k][1]==pred_pairs[k-1][1]+1):
                runs.append((pred_pairs[s][0], pred_pairs[s][1], k-s))
                s=k
    # Hit if any run has length≥tau and overlaps a TM block by ≥tau in both i and j steps (monotonic)
    hit=0
    for (ti,tj,tl) in blocks:
        ok=False
        for (pi,pj,pl) in runs:
            if pl>=tau:
                # check overlap in aligned grid by length
                overlap = min(ti+tl-1, pi+pl-1) - max(ti, pi) + 1
                if overlap is not None and overlap>=tau:
                    ok=True; break
        if ok: hit+=1
    return hit/len(blocks)


    # print_pretty_table removed (unused) to reduce code size per user request


def format_pretty_table(rows, columns):
    """Return the pretty table as a string (same formatting as print_pretty_table)."""
    if not rows:
        return ''
    widths = {}
    for c in columns:
        w = len(c)
        for r in rows:
            v = r.get(c, '')
            s = str(v)
            if '\n' in s:
                s = s.split('\n', 1)[0]
            if len(s) > w:
                w = len(s)
        widths[c] = w + 2

    hdr = ''.join(c.center(widths[c]) for c in columns)
    sep = ''.join('-' * widths[c] for c in columns)
    out_lines = [hdr, sep]
    for r in rows:
        parts = []
        for c in columns:
            v = r.get(c, '')
            s = '' if v is None else str(v)
            if '\n' in s:
                s = s.split('\n', 1)[0]
            if c == 'variant':
                parts.append(s.ljust(widths[c]))
            else:
                parts.append(s.rjust(widths[c]))
        out_lines.append(''.join(parts))
    return '\n'.join(out_lines)


def render_three_row_coverage(indices: List[int], tm_line_chars: List[str], ssw_line_chars: List[str], label_tm: str = 'TM', label_ssw: str = 'SSW') -> str:
    """Render a compact 3-row coverage view for the given index range.

    - indices: list of integers (1-based positions)
    - tm_line_chars: list of single-character strings for TM line (use ' ' for empty)
    - ssw_line_chars: list of single-character strings for SSW line (use ' ' for empty)

    The output is three lines: a numeric index ruler (fixed-width columns),
    a TM line prefixed with label_tm, and an SSW line prefixed with label_ssw.
    Columns are aligned using the width of the largest index.
    """
    if not indices:
        return ''
    max_idx = indices[-1]
    w = max(2, len(str(max_idx)))
    idx_fmt = f"{{:>{w}}}"
    idx_row = 'Idx: ' + ' '.join(idx_fmt.format(i) for i in indices)
    tm_row = label_tm.ljust(4) + ': ' + ' '.join((ch if ch is not None else ' ') .rjust(w) for ch in tm_line_chars)
    ssw_row = label_ssw.ljust(4) + ': ' + ' '.join((ch if ch is not None else ' ') .rjust(w) for ch in ssw_line_chars)
    return '\n'.join([idx_row, tm_row, ssw_row])


# Legacy rendering helpers removed: we produce compact 3-row coverage outputs instead


# render_residuewise_tracks removed — we now emit compact 3-row TSV-like files (Index, TM, SSW)


# _ruler_line removed (legacy)


# _format_range_label removed (legacy)


# render_union_block removed (legacy)


# render_union_tracks removed (legacy)

# ------------------ comparison / metrics ------------------

def compare_pairs(tm_pairs: List[Tuple[int,int]], three_pairs: List[Tuple[int,int]], w: int = 0) -> dict:
    """Compare tm_pairs (ground truth) to three_pairs (predictions) using window w.

    Returns dict with TP, FP, FN, precision, recall, f1, jaccard and counts.

    Matching strategy:
    - For each 3di pair, find any TM pair within |i - i'| <= w and |j - j'| <= w. If found,
      consider the TM pair 'matched' (mark it matched) and mark the 3di pair as matched.
    - TP = number of matched TM pairs
    - FN = total TM pairs - TP
    - FP = total 3di pairs - number of matched 3di pairs
    """
    # Robust matching by sorting on reference (j) coordinate first.
    # Rationale: different tools may produce pairs in different orders or with a
    # global shift; sorting by reference (target) index and using a sliding
    # window reduces spurious misses caused by ordering differences.
    tm_matched = [False] * len(tm_pairs)
    three_matched = [False] * len(three_pairs)
    matched_pairs = []  # list of (tm_idx, three_idx)

    # build indexed lists and sort by reference (j)
    tm_indexed = [(t_idx, it, jt) for t_idx, (it, jt) in enumerate(tm_pairs)]
    three_indexed = [(s_idx, i3, j3) for s_idx, (i3, j3) in enumerate(three_pairs)]
    tm_indexed.sort(key=lambda x: x[2])
    three_indexed.sort(key=lambda x: x[2])

    # sliding window over sorted tm_indexed by j coordinate
    left = 0
    n_tm = len(tm_indexed)
    for s_idx, i3, j3 in three_indexed:
        # advance left to first tm with jt >= j3 - w
        while left < n_tm and tm_indexed[left][2] < j3 - w:
            left += 1
        # scan candidates until jt > j3 + w
        k = left
        while k < n_tm and tm_indexed[k][2] <= j3 + w:
            t_idx, it, jt = tm_indexed[k]
            if not tm_matched[t_idx] and abs(it - i3) <= w and abs(jt - j3) <= w:
                tm_matched[t_idx] = True
                three_matched[s_idx] = True
                matched_pairs.append((t_idx, s_idx))
                break
            k += 1

    TP = sum(1 for x in tm_matched if x)
    matched_three = sum(1 for x in three_matched if x)
    FN = max(0, len(tm_pairs) - TP)
    FP = max(0, len(three_pairs) - matched_three)

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0
    jaccard = TP / (TP + FP + FN) if (TP + FP + FN) > 0 else 0.0

    return {
        'TP': TP, 'FP': FP, 'FN': FN,
        'precision': precision, 'recall': recall, 'f1': f1, 'jaccard': jaccard,
        'n_tm': len(tm_pairs), 'n_3di': len(three_pairs), 'matched_3di': matched_three
        , 'matched_pairs': matched_pairs, 'tm_matched_mask': tm_matched, 'three_matched_mask': three_matched
    }





# ------------------ main pipeline ------------------

def main():
    parser = argparse.ArgumentParser(description='Run TMalign -cp and compare mapping to 3di SSW outputs')
    parser.add_argument('pdb1')
    parser.add_argument('pdb2')
    parser.add_argument('--tmalign', default='TMalign', help='Path to TMalign binary')
    # NOTE: --outdir removed (we persist TM-align copies into run_logs). If you have existing SSW outputs,
    # provide their directory via --alignment-results-dir and they will be copied into the temporary run dir.
    parser.add_argument('--alignment-results-dir', help='Directory containing 3di alignment *_ssw_output.txt files (alignment_results). If provided, these files will be copied into the temporary run alignment_results for reuse')
    parser.add_argument('--skip-tmalign', action='store_true', help='Do not run TMalign; instead load existing TM_align_cp_result.txt from --tmalign-output-file or run_logs/tmp')
    parser.add_argument('--tmalign-output-file', help='If provided, use this file instead of running TMalign')
    parser.add_argument('--model-dir', help='Single model tmp dir (override). If provided, run only this model', default=None)
    
    # (removed) --embed-outdir-in-tmp: no longer supported; outdir is used for summaries only
    parser.add_argument('--yes', action='store_true', dest='yes', help='Auto-confirm prompts')
    parser.add_argument('--gap-open', type=int, default=8, help='SSW gap open')
    parser.add_argument('--gap-extend', type=int, default=2, help='SSW gap extend')
    parser.add_argument('--debug-visual-parse', action='store_true', help='Emit detailed logs while parsing TM-align visual blocks')
    parser.add_argument('--debug-ssw-parse', action='store_true', help='Emit detailed logs while parsing SSW SAM/CIGAR')
    # CP / indexing options
    parser.add_argument('--len-target', type=int, default=None, help='Length of target (B). If omitted, inferred from FASTA for target')
    parser.add_argument('--len-query', type=int, default=None, help='Length of query (A). If omitted, inferred from FASTA for query')
    parser.add_argument('--cutpoint', type=str, default='minj', help='Cutpoint choice: "minj", "lcr", or integer (1-based)')
    parser.add_argument('--w-pair', type=int, default=1, help='Tolerance window for pair matching (use with compare_pairs)')
    parser.add_argument('--tau', type=int, default=8, help='tau for BlockRecall@tau')
    parser.add_argument('--using-bprime', action='store_true', help="Indicates SSW target used B' = B||B and needs unwrapping")
    parser.add_argument('--window-start', type=int, default=None, help="If --using-bprime, the B' window start t (1-based) used for SSW target")
    parser.add_argument('--target-slice', type=int, default=None, help='If SSW used circular slice B[t..]+B[..t-1], provide t (1-based)')
    
    # legacy track-writing options removed; we now always write compact per-variant coverage files into tmp alignment_results
    # allow overriding global ssw location and forcing a single variant
    parser.add_argument('--ssw-root', default=str(Path(ROOT) / 'ssw'), help='Directory containing global ssw_test executable and matrices (default: <root>/ssw)')
    parser.add_argument('--variant', choices=['8f', '9f', '10f'], help='If set, run only the specified model variant (8f/9f/10f)')
    parser.add_argument('--run-dir', help='Explicit run directory to use for temporary outputs (overrides auto tmpdir)')
    args = parser.parse_args()

    # Validate exclusive options: target-slice and using-bprime/window-start
    if args.target_slice is not None and args.using_bprime:
        logger.error('Provide exactly one of --target-slice or --using-bprime (with optional --window-start) when index correction is needed')
        sys.exit(1)

    p1 = Path(args.pdb1)
    p2 = Path(args.pdb2)
    # We no longer accept --outdir; persistent TM-align copies are written to project run_logs
    tmalign_out = _run_logs_dir / f'TM_align_cp_result_{p1.name}_{p2.name}.txt'

    # additional CLI options defaults
    # allow model-dir override, foldseek path, yes flag, gap params
    # These may be provided as environment variables or defaults were already used by caller
    # Auto-detect variants unless --model-dir provided
    model_variant_dirs = []
    if getattr(args, 'model_dir', None):
        provided = Path(args.model_dir)
        provided_abs = provided.resolve()
        m = re.search(r'(8f|9f|10f)', str(provided_abs).lower())
        suffix = m.group(1) if m else provided_abs.name
        model_variant_dirs.append((suffix, provided_abs))
    else:
        variants_to_check = ['8f', '9f', '10f']
        if getattr(args, 'variant', None):
            variants_to_check = [args.variant]
        for v in variants_to_check:
            cand = Path(f'training_3di_gpu_{v}') / 'tmp'
            if cand.exists():
                model_variant_dirs.append((v, cand.resolve()))
        if not model_variant_dirs:
            # fallback to 8f tmp even if missing (preserves previous behavior)
            model_variant_dirs.append(('8f', (Path('training_3di_gpu_8f') / 'tmp').resolve()))

    # create project-local tmp workspace and copy pdbs
    project_tmp_root = Path(ROOT) / 'tmp'
    project_tmp_root.mkdir(parents=True, exist_ok=True)
    # Use stable run directory name based on input stems: {proA}_vs_{proB}
    run_name = f"{p1.stem}_vs_{p2.stem}"
    # Allow caller to supply an explicit run directory to control where outputs go
    if getattr(args, 'run_dir', None):
        tmpdir = Path(args.run_dir)
        # if relative, make absolute inside project tmp root
        if not tmpdir.is_absolute():
            tmpdir = project_tmp_root / tmpdir
        tmpdir.parent.mkdir(parents=True, exist_ok=True)
        tmpdir.mkdir(parents=True, exist_ok=True)
    else:
        tmpdir = project_tmp_root / run_name
        if tmpdir.exists():
            # if exists, append timestamp to avoid clobbering
            tmpdir = project_tmp_root / f"{run_name}_{time.strftime('%Y%m%d%H%M%S')}"
        tmpdir.mkdir(parents=True, exist_ok=False)
    logger.info(f'Using run directory: {tmpdir}')
    # copy input pdbs into tmpdir
    p1_tmp = tmpdir / p1.name
    p2_tmp = tmpdir / p2.name
    shutil.copy2(p1, p1_tmp)
    shutil.copy2(p2, p2_tmp)

    # Note: --outdir was removed; persistent TM-align copies go to run_logs

    # Step 2: run foldseek createdb on each PDB (best-effort)
    # --foldseek CLI was removed; assume 'foldseek' is available in PATH
    foldseek_bin = 'foldseek'
    for ptmp in (p1_tmp, p2_tmp):
        prefix = ptmp.stem + '_pdzdb'
        p_sub = tmpdir / ptmp.stem
        p_sub.mkdir(parents=True, exist_ok=True)
        out_base_name = prefix
        rc, out, err = run_cmd([foldseek_bin, 'createdb', str(ptmp), out_base_name], cwd=str(p_sub))
        if rc != 0:
            logger.warning(f'foldseek createdb failed for {ptmp} (rc={rc}); continuing')
        else:
            logger.info(f'foldseek createdb done for {ptmp} -> {p_sub}')

    # Step 3: process createdb outputs -> write cleaned FASTA and 10f_3di files
    def process_createdb_dir(p_sub: Path, prefix: str, orig_name: str):
        fasta_path = p_sub / f"{orig_name}.fasta"
        tenf_txt = p_sub / f"{orig_name}_10f.txt"
        # produce variant-named 3di FASTAs: <orig>_<variant>_3di.fasta and _x2
        three_fa = p_sub / f"{orig_name}_10f_3di.fasta"
        three_fa_x2 = p_sub / f"{orig_name}_10f_3di_x2.fasta"
        search_dir = p_sub if p_sub.exists() else p_sub.parent
        candidate_pdz = search_dir / prefix
        candidate_pdz_ss = search_dir / f"{prefix}_ss"
        candidate_pdz_ss_alt = search_dir / f"{prefix}_pdzdb_ss"
        pdz = None
        pdz_ss = None
        if candidate_pdz.exists():
            pdz = candidate_pdz
        if candidate_pdz_ss.exists():
            pdz_ss = candidate_pdz_ss
        elif candidate_pdz_ss_alt.exists():
            pdz_ss = candidate_pdz_ss_alt
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
        try:
            p_sub.mkdir(parents=True, exist_ok=True)
        except Exception:
            pass
        if pdz and pdz.exists():
            text = cleanse_createdb_file(pdz)
            seq_text = extract_sequence_from_file_content(text)
            write_fasta_from_string(fasta_path, orig_name, seq_text)
            logger.info(f'Wrote FASTA: {fasta_path}')
        else:
            logger.warning(f'pdz file not found for {orig_name} in {search_dir}')

        # Handle *_ss token output produced by createdb: save as <stem>_10f.txt and generate 3di FASTA files
        if pdz_ss and pdz_ss.exists():
            text = cleanse_createdb_file(pdz_ss)
            token_seq = extract_sequence_from_file_content(text)
            # write single-line 10f token file
            # write single-line 10f token file temporarily, then remove it to keep only FASTAs
            try:
                tenf_txt.write_text(token_seq + '\n')
                logger.info(f'Wrote temporary 10f token file: {tenf_txt}')
            except Exception:
                logger.debug('Failed to write temporary tenf txt', exc_info=True)
            # write 10f FASTA (single) and x2 FASTA using *_10f_3di naming (no legacy _3di files)
            tenf_fa = p_sub / f"{orig_name}_10f_3di.fasta"
            tenf_fa_x2 = p_sub / f"{orig_name}_10f_3di_x2.fasta"
            write_3di_fasta(tenf_fa, orig_name + '_10f_3di', token_seq)
            write_3di_fasta(tenf_fa_x2, orig_name + '_10f_3di_x2', token_seq + token_seq)
            # also set three_fa/three_fa_x2 to these files for uniform downstream names
            three_fa = tenf_fa
            three_fa_x2 = tenf_fa_x2
            logger.info(f'Wrote 10f FASTA: {tenf_fa} and x2: {tenf_fa_x2}')
            # remove the temporary token txt to preserve only FASTA files per-protein
            try:
                if tenf_txt.exists():
                    tenf_txt.unlink()
                    logger.info(f'Removed temporary token file: {tenf_txt}')
            except Exception:
                logger.debug('Failed to remove temporary tenf txt', exc_info=True)
        else:
            logger.warning(f'pdz_ss (3di) file not found for {orig_name} in {search_dir}')
        # cleanup: remove createdb intermediate files in this dir (keep fasta_path, tenf_txt, tenf_fa, tenf_fa_x2 and legacy three_fa files)
        try:
            for f in list(search_dir.iterdir()):
                try:
                    # keep AA fasta and 10f token txt / 3di fasta files
                    # keep AA fasta, 10f token txt, and the produced variant 3di FASTA/x2
                    if f.resolve() == fasta_path.resolve() or f.resolve() == tenf_txt.resolve() or f.resolve() == three_fa.resolve() or f.resolve() == three_fa_x2.resolve():
                        continue
                    # remove files or directories produced by createdb
                    if f.is_file():
                        f.unlink()
                    else:
                        shutil.rmtree(f)
                except Exception:
                    logger.debug(f'Failed to remove {f} during cleanup', exc_info=True)
        except Exception:
            logger.debug('Createdb cleanup failed', exc_info=True)

        # Return aa fasta, 10f token txt, 10f fasta and 10f x2 fasta
        tenf_fa = p_sub / f"{orig_name}_10f_3di.fasta"
        tenf_fa_x2 = p_sub / f"{orig_name}_10f_3di_x2.fasta"
        return fasta_path, tenf_txt, tenf_fa, tenf_fa_x2

    p1_sub = tmpdir / p1_tmp.stem
    p2_sub = tmpdir / p2_tmp.stem
    prefix1 = p1_tmp.stem + '_pdzdb'
    prefix2 = p2_tmp.stem + '_pdzdb'
    # process_createdb_dir now returns: (aa_fasta, tenf_txt, tenf_fa, tenf_fa_x2)
    fasta1, tenf1, tenf_fa1, tenf_fa1_x2 = process_createdb_dir(p1_sub, prefix1, p1_tmp.stem)
    fasta2, tenf2, tenf_fa2, tenf_fa2_x2 = process_createdb_dir(p2_sub, prefix2, p2_tmp.stem)

    # Step 4: run pdb_to_3di.py for each variant to produce new3di files
    # Note: do NOT call pdb_to_3di for 10f — its tokens come from createdb output
    new3di_by_var = {}
    for suffix, model_tmp in model_variant_dirs:
        if suffix == '10f':
            # skip calling pdb_to_3di for 10f; we'll register 10f from createdb outputs below
            continue
        for ptmp in (p1_tmp, p2_tmp):
            p_sub = tmpdir / ptmp.stem
            p_sub.mkdir(parents=True, exist_ok=True)
            out_name = f"{ptmp.stem}_{suffix}_3di.txt"
            out_path = p_sub / out_name
            model_dir_abs = str(model_tmp)
            rc, out, err = run_cmd([sys.executable, str(Path(__file__).parent / 'pdb_to_3di.py'), str(ptmp), model_dir_abs, '-o', str(out_path), '--no-header'], cwd=str(p_sub))
            if rc != 0:
                logger.warning(f'pdb_to_3di failed for variant {suffix} on {ptmp} (rc={rc}); skipping this variant')
                break
        f1 = tmpdir / p1_tmp.stem / f"{p1_tmp.stem}_{suffix}_3di.txt"
        f2 = tmpdir / p2_tmp.stem / f"{p2_tmp.stem}_{suffix}_3di.txt"
        if f1.exists() and f2.exists():
            # produce persistent per-protein FASTA files named <stem>_<variant>_3di.fasta and _x2.fasta
            p1_sub = tmpdir / p1_tmp.stem
            p2_sub = tmpdir / p2_tmp.stem
            token1 = extract_sequence_from_file_content(f1.read_text())
            token2 = extract_sequence_from_file_content(f2.read_text())
            fa1 = p1_sub / f"{p1_tmp.stem}_{suffix}_3di.fasta"
            fa1_x2 = p1_sub / f"{p1_tmp.stem}_{suffix}_3di_x2.fasta"
            fa2 = p2_sub / f"{p2_tmp.stem}_{suffix}_3di.fasta"
            fa2_x2 = p2_sub / f"{p2_tmp.stem}_{suffix}_3di_x2.fasta"
            write_3di_fasta(fa1, fa1.stem, token1)
            write_3di_fasta(fa1_x2, fa1_x2.stem, token1 + token1)
            write_3di_fasta(fa2, fa2.stem, token2)
            write_3di_fasta(fa2_x2, fa2_x2.stem, token2 + token2)
            # remove the intermediate *_3di.txt files so only FASTAs remain per-protein
            try:
                if f1.exists():
                    f1.unlink()
                    logger.info(f'Removed intermediate token file: {f1}')
            except Exception:
                logger.debug(f'Failed to remove intermediate file {f1}', exc_info=True)
            try:
                if f2.exists():
                    f2.unlink()
                    logger.info(f'Removed intermediate token file: {f2}')
            except Exception:
                logger.debug(f'Failed to remove intermediate file {f2}', exc_info=True)
            new3di_by_var[suffix] = (fa1, fa2)
        else:
            logger.warning(f'new3di files for variant {suffix} not found; skipping')

    # Register 10f variant directly from createdb *_ss token outputs (tenf files)
    try:
        # register 10f using the generated 10f FASTA files (we removed temporary txt tokens earlier)
        if tenf_fa1.exists() and tenf_fa2.exists():
            new3di_by_var['10f'] = (tenf_fa1, tenf_fa2)
            # remember the pre-created 10f fasta and x2 files for use as direct SSW inputs
            three_10f_map = {'A': tenf_fa1, 'A_x2': tenf_fa1_x2, 'B': tenf_fa2, 'B_x2': tenf_fa2_x2}
            logger.info(f'Registered 10f variant from createdb outputs: {tenf_fa1}, {tenf_fa2}')
        else:
            logger.error(f'10f FASTA files not found: expected {tenf_fa1} and {tenf_fa2}; cannot register 10f variant')
    except Exception:
        logger.exception('Failed to register 10f variant from createdb outputs')

    if not new3di_by_var:
        logger.error('No variant produced 3di outputs; aborting')
        sys.exit(1)

    # Helper: search for an existing TM-align result file that matches the current pdb pair.
    def find_existing_tmalign(p1stem: str, p2stem: str) -> Optional[Path]:
        # Search order: run_logs, tmp/tmalign_3di_* directories
        # If a candidate file contains both pdb stems in its text, treat as a match.
        candidates = []
        try:
            if tmalign_out.exists():
                candidates.append(tmalign_out)
        except Exception:
            pass
        # run_logs under project root
        try:
            for f in sorted(_run_logs_dir.glob('TM_align_cp_result*.txt')):
                candidates.append(f)
        except Exception:
            pass
        # tmp/tmalign_3di_* run dirs
        try:
            for d in sorted((Path(ROOT) / 'tmp').glob('tmalign_3di_*')):
                p = d / 'TM_align_cp_result.txt'
                if p.exists():
                    candidates.append(p)
        except Exception:
            pass

        for cand in candidates:
            try:
                txt = cand.read_text()
                if p1stem in txt and p2stem in txt:
                    return cand
            except Exception:
                continue
        return None

    # 1) Run or load TM-align -cp output
    if args.tmalign_output_file:
        provided_path = Path(args.tmalign_output_file)
        if provided_path.exists():
            try:
                tmalign_text = provided_path.read_text()
                print('Loaded TM-align output from', args.tmalign_output_file)
                try:
                    tmp_p = tmpdir / provided_path.name
                    tmp_p.write_text(tmalign_text)
                except Exception:
                    logger.debug('Failed to write TM-align copy into tmpdir for provided file', exc_info=True)
            except Exception:
                logger.warning(f'Could not read provided --tmalign-output-file: {args.tmalign_output_file}; will search for existing outputs or run TMalign')
                tmalign_text = None
        else:
            logger.warning(f'Provided --tmalign-output-file not found: {args.tmalign_output_file}; will search for existing outputs or run TMalign')
            tmalign_text = None
    else:
        # search for an existing TM-align result matching this pair
        found = find_existing_tmalign(p1.name, p2.name)
        if found:
            tmalign_text = found.read_text()
            print('Loaded existing TM-align output from', str(found))
            try:
                tmp_t = tmpdir / found.name
                tmp_t.write_text(tmalign_text)
            except Exception:
                logger.debug('Failed to write existing TM-align into tmpdir', exc_info=True)
        else:
            if args.skip_tmalign:
                print('Error: --skip-tmalign specified but no matching TM_align_cp_result.txt found in run_logs/tmp')
                sys.exit(1)
                sys.exit(1)
            # run TMalign and write result into run_logs for future reuse
            cmd = [args.tmalign, str(p1), str(p2), '-cp']
            rc, out, err = run_cmd(cmd)
            if rc != 0:
                print('TMalign failed (rc=%d). stderr:' % rc)
                print(err)
                sys.exit(1)
            tmalign_text = out + '\n' + err
            try:
                tmalign_out.parent.mkdir(parents=True, exist_ok=True)
                tmalign_out.write_text(tmalign_text)
            except Exception:
                logger.debug('Failed to write TM-align output into run_logs', exc_info=True)
            # also keep a copy in tmpdir for per-run inspection
            try:
                tmp_tmalign = tmpdir / tmalign_out.name
                tmp_tmalign.write_text(tmalign_text)
            except Exception:
                logger.debug('Failed to write TM-align copy into tmpdir', exc_info=True)
            print('Wrote TM-align -cp output to', str(tmalign_out))

    # Prefer visual-block parsing to extract ordered residue-pair mapping from TM-align
    vis_pairs, vis_pairs_labeled, printed_top_seq, printed_bot_seq = parse_tmalign_visual_pairs(tmalign_text, debug=getattr(args, 'debug_visual_parse', False))
    if vis_pairs:
        # trim to aligned length if available and deduplicated
        parsed_pairs = vis_pairs
        m = re.search(r'Aligned length=\s*(\d+)', tmalign_text)
        if m:
            aln_len = int(m.group(1))
            if len(parsed_pairs) > aln_len:
                logger.info(f'Parsed {len(parsed_pairs)} TM-align pairs from visual blocks; trimming to Aligned length={aln_len}')
                parsed_pairs = parsed_pairs[:aln_len]
        logger.info(f'Parsed {len(parsed_pairs)} TM-align pairs from visual alignment blocks')

        # map printed indices back to original fasta indices, attempting to detect circular shift
        def find_rotation_offset(original: str, printed: str) -> int:
            """Return 0-based rotation offset that maps printed (TM printed subsequence)
            into the original FASTA sequence. If not found return -1.

            This is more robust than a plain substring search: it first strips
            non-letter characters from the printed string (TM output may include
            spaces, '*' markers, or other punctuation), upper-cases both
            sequences, and searches the doubled original. If the full printed
            string is not found, try a small k-mer fallback to locate an
            approximate anchor and compute the offset.
            """
            import re
            if not printed or not original:
                return -1
            # sanitize inputs: keep A-Z letters only and uppercase
            printed_clean = re.sub(r'[^A-Za-z]', '', printed).upper()
            original_clean = re.sub(r'[^A-Za-z]', '', original).upper()
            if not printed_clean or not original_clean:
                return -1

            doubled = original_clean + original_clean
            # 1) try full-match of the cleaned printed string inside doubled original
            pos = doubled.find(printed_clean)
            if pos != -1 and pos < len(original_clean):
                return pos  # 0-based offset

            # 1b) If the printed sequence contains an explicit CP marker '*' then
            # use the residue immediately after '*' as an anchor. The printed
            # index 1 corresponds to the residue after '*', so find candidate
            # anchor positions in the original and validate by short k-mer match.
            if '*' in printed:
                # compute how many letters appear before the '*' (in printed_clean)
                letters_before_star = 0
                for ch in printed:
                    if ch == '*':
                        break
                    if ch.isalpha():
                        letters_before_star += 1
                # anchor in printed_clean is the letter at index letters_before_star
                if letters_before_star < len(printed_clean):
                    anchor_char = printed_clean[letters_before_star]
                    # search for anchor_char occurrences in original_clean
                    candidates = [m.start() for m in re.finditer(re.escape(anchor_char), original_clean)]
                    # choose k for local match
                    k = min(12, max(4, len(printed_clean) // 4))
                    for cand in candidates:
                        # assume printed_clean[letters_before_star] maps to original_clean[cand]
                        # then printed_clean[0] would map to original index (cand - letters_before_star)
                        off = (cand - letters_before_star) % len(original_clean)
                        # validate by comparing up to k residues from that offset
                        ok = True
                        for idx in range(k):
                            p_idx = (idx) % len(printed_clean)
                            o_idx = (off + idx) % len(original_clean)
                            if printed_clean[p_idx] != original_clean[o_idx]:
                                ok = False
                                break
                        if ok:
                            return off

            # 2) fallback: try k-mer anchors from the printed string (no '*' present)
            plen = len(printed_clean)
            k = min(10, max(3, plen // 3))
            for start in range(0, max(1, plen - k + 1)):
                kmer = printed_clean[start:start + k]
                pos_k = doubled.find(kmer)
                if pos_k != -1 and pos_k < len(original_clean):
                    # assume printed starts at pos_k - start
                    off = (pos_k - start) % len(original_clean)
                    return off

            # not found
            return -1

        seqA = read_fasta_sequence(fasta1)
        seqB = read_fasta_sequence(fasta2)

        # Determine which TM-align chain corresponds to the printed top/bot lines.
        # TM-align prints Chain_1 on the top line and Chain_2 on the bottom line.
        # We need to map printed_top_seq -> the FASTA for the chain that was printed on top.
        tm_chain1 = None
        tm_chain2 = None
        try:
            m1 = re.search(r'Name of Chain_1:\s*(\S+)', tmalign_text)
            m2 = re.search(r'Name of Chain_2:\s*(\S+)', tmalign_text)
            if m1:
                tm_chain1 = Path(m1.group(1)).stem
            if m2:
                tm_chain2 = Path(m2.group(1)).stem
        except Exception:
            tm_chain1 = None
            tm_chain2 = None

        # Default: seq_top maps to fasta1 and seq_bot to fasta2 (previous behaviour)
        seq_top = seqA
        seq_bot = seqB
        # If TM-align reports Chain_1/Chain_2, align them to our input pdb stems and swap if needed
        try:
            p1stem = p1_tmp.stem
            p2stem = p2_tmp.stem
            if tm_chain1 and tm_chain2:
                # If Chain_1 matches p2 (second input), the printed top is actually fasta2
                if tm_chain1 == p2stem and tm_chain2 == p1stem:
                    seq_top = seqB
                    seq_bot = seqA
                # If Chain_1 matches p1 as expected, keep default mapping
                elif tm_chain1 == p1stem and tm_chain2 == p2stem:
                    seq_top = seqA
                    seq_bot = seqB
                else:
                    # fallback: try to match any one chain name
                    if tm_chain1 == p2stem:
                        seq_top = seqB
                        seq_bot = seqA
                    elif tm_chain1 == p1stem:
                        seq_top = seqA
                        seq_bot = seqB
        except Exception:
            seq_top = seqA
            seq_bot = seqB

        logger.info(f'Mapping printed TM top/bot to FASTA stems: top->(len={len(seq_top)}) bot->(len={len(seq_bot)})')

        # Build printed->original index maps for both sequences so that parsed
        # visual-block printed indices can be immediately translated to original
        # FASTA coordinates. This follows the same logic as find_rotation_offset
        # but stores a full mapping table and applies '*' anchoring and k-mer
        # fallbacks when necessary.
        def build_printed_to_original_map(original: str, printed: str):
            # returns a list (1-based index -> original_index) of length = len(printed_clean)
            import re
            printed_clean = re.sub(r'[^A-Za-z]', '', printed).upper()
            original_clean = re.sub(r'[^A-Za-z]', '', original).upper()
            n = len(printed_clean)
            if n == 0 or len(original_clean) == 0:
                return None
            L = len(original_clean)
            # try to find a global rotation offset first
            off = find_rotation_offset(original, printed)
            if off >= 0:
                # build mapping
                return [((off + (p - 1)) % L) + 1 for p in range(1, n + 1)]

            # if printed contains a CP marker '*', try star-anchored mapping
            if '*' in printed:
                # count letters before the '*' in the printed (raw) string
                letters_before_star = 0
                for ch in printed:
                    if ch == '*':
                        break
                    if ch.isalpha():
                        letters_before_star += 1
                if letters_before_star < n:
                    anchor_char = printed_clean[letters_before_star]
                    cand_positions = [m.start() for m in re.finditer(re.escape(anchor_char), original_clean)]
                    k = min(12, max(4, n // 4))
                    for cand in cand_positions:
                        off_cand = (cand - letters_before_star) % L
                        # validate k residues
                        ok = True
                        for idx in range(min(k, n)):
                            p_idx = idx
                            o_idx = (off_cand + idx) % L
                            if printed_clean[p_idx] != original_clean[o_idx]:
                                ok = False
                                break
                        if ok:
                            return [((off_cand + (p - 1)) % L) + 1 for p in range(1, n + 1)]

            # fallback: try k-mer anchors across printed_clean
            plen = n
            k = min(10, max(3, plen // 3))
            doubled = original_clean + original_clean
            for start in range(0, max(1, plen - k + 1)):
                kmer = printed_clean[start:start + k]
                pos_k = doubled.find(kmer)
                if pos_k != -1 and pos_k < L:
                    off_cand = (pos_k - start) % L
                    return [((off_cand + (p - 1)) % L) + 1 for p in range(1, n + 1)]

            # last resort: no mapping found
            return None

        # Use seq_top/seq_bot which are selected to correspond to printed_top_seq/printed_bot_seq
        mapA = build_printed_to_original_map(seq_top, printed_top_seq)
        mapB = build_printed_to_original_map(seq_bot, printed_bot_seq)
        if mapA is None and mapB is None:
            logger.info('Could not map printed sequences back to original FASTA positions; keeping printed indices')

        mapped_pairs = []
        for (ia, jb) in parsed_pairs:
            ma = ia
            mb = jb
            # if a map exists and the printed index is within its length, use it
            if mapA and 1 <= ia <= len(mapA):
                ma = mapA[ia - 1]
            if mapB and 1 <= jb <= len(mapB):
                mb = mapB[jb - 1]
            mapped_pairs.append((ma, mb))
            # Validate mapped indices against actual sequence lengths. If a large fraction
            # of mapped indices fall outside the expected ranges, the printed->original
            # mapping likely failed (common with odd TM visual formatting). In that case
            # fall back to the numeric TM-align parser which extracts explicit integer
            # pairs from the TM-align output.
            seqA = read_fasta_sequence(fasta1)
            seqB = read_fasta_sequence(fasta2)
            lenA = len(seqA) if seqA else None
            lenB = len(seqB) if seqB else None
            invalid = 0
            for (ma, mb) in mapped_pairs:
                if ma is None or mb is None:
                    invalid += 1
                    continue
                if (isinstance(lenA, int) and (ma < 1 or ma > lenA)) or (isinstance(lenB, int) and (mb < 1 or mb > lenB)):
                    invalid += 1
            if mapped_pairs and (invalid / len(mapped_pairs) > 0.2):
                logger.warning('Mapped TM indices appear invalid (%.1f%%); falling back to numeric TM-align parser' % (100.0 * invalid / len(mapped_pairs)))
                tm_pairs = parse_tmalign_cp_pairs(tmalign_text)
            else:
                tm_pairs = mapped_pairs
    else:
        # fallback to heuristic cp parser
        tm_pairs = parse_tmalign_cp_pairs(tmalign_text)
        logger.info('Used fallback heuristic TM-align parser')
    print('Parsed %d TM-align pairs (from visual/block parsing if available)' % len(tm_pairs))
    # prepare CP-aware rotation and length info
    seqA = read_fasta_sequence(fasta1) if 'fasta1' in locals() else ''
    seqB = read_fasta_sequence(fasta2) if 'fasta2' in locals() else ''
    lenT = args.len_target if getattr(args, 'len_target', None) else (len(seqB) if seqB else None)
    lenQ = args.len_query if getattr(args, 'len_query', None) else (len(seqA) if seqA else None)
    if lenT is None:
        logger.warning('len_target could not be inferred; defaulting to max j in tm_pairs or 1')
        lenT = max((j for _,j in tm_pairs), default=1)
    if lenQ is None:
        lenQ = max((i for i,_ in tm_pairs), default=1)

    # choose cutpoint s based on tm_pairs and requested mode
    cutmode = getattr(args, 'cutpoint', 'minj')
    try:
        s = choose_cutpoint_from_tm(tm_pairs, cutmode, lenT)
    except Exception:
        s = 1
    logger.info(f'Chosen cutpoint s={s} (mode={cutmode}) len_target={lenT}')

    # rotate TM j indices to cutpoint s
    tm_rot = [(i, rotate_j(j, s, lenT)) for (i,j) in tm_pairs]
    # also keep a sorted version for contiguity checks
    tm_rot_sorted = sorted(tm_rot)
    # write TM pairs to a simple TSV inside tmpdir for inspection
    try:
        tm_pairs_file = tmpdir / 'tm_pairs.tsv'
        with tm_pairs_file.open('w') as tf:
            # include original protein/file names in header
            tf.write(f"#pdb1\t{p1.name}\n#pdb2\t{p2.name}\n#i\tj\n")
            for i,j in tm_pairs:
                tf.write(f"{i}\t{j}\n")
        logger.info(f'Wrote TM pairs to {tm_pairs_file}')
        # also write a TM-align snippet for quick visual inspection
        # Do not write a TM-align snippet file by default (user requested to avoid creating
        # TM_align_cp_result_snippet.txt). Keep the tm_align_text available in-memory.
    except Exception:
        logger.debug('Failed to write TM pairs file into tmpdir', exc_info=True)

    # Run SSW for each generated 3di variant (use full 3di sequences, no trimming)
    # SSW outputs are written only into the run tmp directory under alignment_results
    tmp_align_dir = tmpdir / 'alignment_results'
    tmp_align_dir.mkdir(parents=True, exist_ok=True)
    align_dir = tmp_align_dir

    # If user previously stored alignment_results somewhere, optionally copy those into tmp so
    # we can skip re-running SSW. We will not write new SSW outputs into outdir.
    existing_align_dir = Path(args.alignment_results_dir) if args.alignment_results_dir else None
    if existing_align_dir and existing_align_dir.exists():
        for f in sorted(existing_align_dir.glob('*_ssw_output.txt')):
            if p1.stem in f.name and p2.stem in f.name:
                try:
                    shutil.copy2(f, tmp_align_dir / f.name)
                    logger.info(f'Copied existing SSW output from {f} into {tmp_align_dir}')
                except Exception:
                    logger.debug(f'Failed copying existing SSW {f} into tmpdir', exc_info=True)

    # build a lookup for model tmp dirs
    model_tmp_by_suffix = {s: p for s, p in model_variant_dirs}

    for suffix, (f1, f2) in new3di_by_var.items():
        print(f'Running ssw for variant {suffix}...')
        # locate model tmp for this suffix and global ssw root
        model_tmp = model_tmp_by_suffix.get(suffix)
        ssw_root = Path(getattr(args, 'ssw_root', Path(ROOT) / 'ssw'))
        global_ssw_exec = ssw_root / 'ssw_test'
        # prefer global ssw_test if available, else fall back to model tmp
        if global_ssw_exec.exists():
            ssw_exec = str(global_ssw_exec.resolve())
        else:
            if model_tmp is None:
                logger.warning(f'No model tmp found for variant {suffix}; skipping ssw')
                continue
            ssw_exec = str(Path(model_tmp) / 'ssw_test')

        # Special-case: for 10f we MUST use s_10f.mat from ssw_root and use pre-created 3di FASTA files.
        matrix_abs = None
        if suffix == '10f':
            s10 = ssw_root / 's_10f.mat'
            if not Path(ssw_exec).exists():
                logger.error(f'ssw_test not found at expected location: {ssw_exec}. Provide a working ssw_test under --ssw-root or in the model tmp dir.')
                sys.exit(1)
            if not s10.exists():
                logger.error(f'Required matrix s_10f.mat not found under {ssw_root}. Please place s_10f.mat there and retry.')
                sys.exit(1)
            matrix_abs = str(s10.resolve())
            logger.info(f'Using ssw executable: {ssw_exec} for variant {suffix}')
            logger.info(f'Using substitution matrix (forced for 10f): {matrix_abs} for variant {suffix}')
            # use the pre-created 3di fasta/x2 files from createdb
            try:
                fasta_a = three_10f_map['A']
                fasta_b = three_10f_map['B_x2']
            except Exception:
                logger.error('Pre-created 10f 3Di FASTA files not found; ensure createdb produced *_3di.fasta and *_3di_x2.fasta')
                sys.exit(1)
            cmd_unified = [str(ssw_exec), '-p', '-a', matrix_abs]
        else:
            # For non-10f variants: prefer to use the per-protein FASTA files produced earlier
            # (these live in the per-protein subfolder). Avoid creating duplicate FASTAs in tmpdir.
            try:
                # f1/f2 are paths to per-protein 3di FASTA files registered in new3di_by_var
                fasta_a = Path(f1)
                fasta_b_base = Path(f2)
                # prefer an existing _x2 fasta near fasta_b; if missing, create one alongside fasta_b
                fasta_b_x2 = fasta_b_base.parent / (fasta_b_base.stem + '_x2.fasta')
                if fasta_b_x2.exists():
                    fasta_b = fasta_b_x2
                else:
                    # create x2 file next to fasta_b_base
                    try:
                        seqb = read_fasta_sequence(fasta_b_base)
                        write_fasta_from_string(fasta_b_x2, fasta_b_x2.stem, seqb + seqb)
                        fasta_b = fasta_b_x2
                        logger.info(f'Created doubled query FASTA at {fasta_b_x2}')
                    except Exception:
                        # fallback: create in tmpdir
                        seqb = extract_sequence_from_file_content(f2.read_text())
                        fasta_b = tmpdir / (p2.stem + f'_{suffix}_3di_x2.fasta')
                        write_fasta_from_string(fasta_b, fasta_b.stem, seqb + seqb)

            except Exception:
                # if anything fails, fallback to writing FASTAs into tmpdir (legacy behaviour)
                fasta_a = tmpdir / (p1.stem + f'_{suffix}_3di.fasta')
                fasta_b = tmpdir / (p2.stem + f'_{suffix}_3di_x2.fasta')
                write_fasta_from_string(fasta_a, fasta_a.stem, extract_sequence_from_file_content(f1.read_text()))
                seqb = extract_sequence_from_file_content(f2.read_text())
                write_fasta_from_string(fasta_b, fasta_b.stem, seqb + seqb)

            # choose matrix file: prefer global s_<suffix>.mat in ssw_root, then global s.mat/sub_score.mat,
            # then model_tmp/sub_score.mat or model_tmp/s.mat
            try:
                candidates = []
                candidates.append(ssw_root / f's_{suffix}.mat')
                candidates.append(ssw_root / 's.mat')
                candidates.append(ssw_root / 'sub_score.mat')
                if model_tmp is not None:
                    candidates.append(Path(model_tmp) / 'sub_score.mat')
                    candidates.append(Path(model_tmp) / 's.mat')
                for cand in candidates:
                    if cand.exists():
                        matrix_abs = str(cand.resolve())
                        try:
                            shutil.copy2(cand, _run_logs_dir / cand.name)
                        except Exception:
                            pass
                        break
            except Exception:
                matrix_abs = None

            logger.info(f'Using ssw executable: {ssw_exec} for variant {suffix}')
            logger.info(f'Using substitution matrix: {matrix_abs if matrix_abs else "(none found)"} for variant {suffix}')

            cmd_unified = [ssw_exec, '-p']
            if matrix_abs:
                cmd_unified += ['-a', matrix_abs]
        # use reasonable gap params
        gap_open = getattr(args, 'gap_open', 8)
        gap_ext = getattr(args, 'gap_extend', 2)
        cmd_unified += ['-o', str(gap_open), '-e', str(gap_ext), '-c', str(fasta_a), str(fasta_b)]

        out_ssw_name = f"{fasta_a.stem}_{fasta_b.stem}_{suffix}_ssw_output.txt"
        tmp_out_ssw = tmp_align_dir / out_ssw_name

        # If a persistent alignment_results dir was provided, check there for an existing SSW output and copy it into tmp
        persistent_out_ssw = None
        if existing_align_dir:
            persistent_out_ssw = existing_align_dir / out_ssw_name
        if persistent_out_ssw and persistent_out_ssw.exists():
            try:
                shutil.copy2(persistent_out_ssw, tmp_out_ssw)
                logger.info(f'Found existing SSW output for {suffix} at {persistent_out_ssw}; copied into tmp and skipping SSW run')
            except Exception:
                logger.debug(f'Failed copying existing {persistent_out_ssw} into tmp', exc_info=True)
            # also copy BLAST-like if present in provided alignment dir
            if existing_align_dir:
                out_blast_p = existing_align_dir / f"{fasta_a.stem}_{fasta_b.stem}_{suffix}_blast_output.txt"
                if out_blast_p.exists():
                    try:
                        shutil.copy2(out_blast_p, tmp_align_dir / out_blast_p.name)
                    except Exception:
                        logger.debug(f'Failed copying existing blast {out_blast_p} into tmp', exc_info=True)
            continue

        rc, out, err = run_cmd(cmd_unified, capture=True)
        if rc == 0:
            # write only into tmp alignment_results (do not persist into user-specified outdir)
            tmp_out_ssw.write_text(out)
            logger.info(f'Wrote ssw (-c) output for {suffix} to {tmp_out_ssw}')
            # try converting to BLAST-like using scripts/sam_to_blastlike.py if present
            sam2blast = Path(__file__).parent / 'scripts' / 'sam_to_blastlike.py'
            if sam2blast.exists():
                try:
                    rc_c, out_c, err_c = run_cmd([sys.executable, str(sam2blast), str(fasta_a), str(tmp_out_ssw)], capture=True)
                    if rc_c == 0:
                        # write BLAST-like only to tmp alignment_results
                        out_blast = tmp_align_dir / f"{fasta_a.stem}_{fasta_b.stem}_{suffix}_blast_output.txt"
                        out_blast.write_text(out_c)
                        logger.info(f'Wrote BLAST-like output for {suffix} to {out_blast}')
                    else:
                        logger.warning(f'sam_to_blastlike failed rc={rc_c} stderr={err_c}')
                except Exception:
                    logger.exception('Failed to invoke sam_to_blastlike')
        else:
            logger.warning(f'ssw (-c) failed for {suffix} rc={rc}; stderr: {err}')

    # Ensure we search the alignment_results we just produced unless the user provided a different dir
    if not args.alignment_results_dir:
        args.alignment_results_dir = str(align_dir)

    # 2) Find 3di ssw output files
    ssw_files = []
    if args.alignment_results_dir:
        ar = Path(args.alignment_results_dir)
        if not ar.exists():
            print('alignment-results-dir not found:', ar)
            sys.exit(1)
        for f in sorted(ar.glob('*_ssw_output.txt')):
            # crude filter: require both pdb stems to be in filename
            name = f.name
            if p1.stem in name and p2.stem in name:
                ssw_files.append(f)
    else:
        # search pairwise_3di_* directories
        for d in sorted(Path('.').glob('pairwise_3di_*')):
            ar = d / 'alignment_results'
            if not ar.exists():
                continue
            for f in sorted(ar.glob('*_ssw_output.txt')):
                if p1.stem in f.name and p2.stem in f.name:
                    ssw_files.append(f)
    if not ssw_files:
        print('No *_ssw_output.txt files found for the given PDB stems. Provide --alignment-results-dir or run pairwise pipeline first.')
        sys.exit(1)

    logger.info('Found %d ssw output files to compare' % len(ssw_files))

    # 3) For each ssw file, parse pairs and compare
    rows = []
    for f in ssw_files:
        logger.info(f'Parsing SSW file: {f}')
        text = f.read_text()
        # attempt to infer original (non-doubled) length if not provided: look for a fasta in same dir named with _trim.fasta
        # --query-duped-len CLI removed; always infer from SAM sequence when possible
        q_duped_len = None
        # heuristic: look for a fasta with same stem but without _x2
        # we assume fasta name contains the second pdb stem
        possible_fastas = list(f.parent.glob(f'*{p2.stem}*_{p2.stem}_*_trim.fasta')) if False else []
        # fallback: try to detect seq length from SAM sequence field
        if q_duped_len is None:
            # find first SAM line and get seq length
            for ln in text.splitlines():
                if ln.strip() and not ln.startswith('@'):
                    fields = ln.split('\t')
                    if len(fields) >= 10:
                        seq = fields[9]
                        if seq:
                            q_duped_len = len(seq) // 2  # assume doubled
                    break
        pairs_3di = parse_ssw_sam_pairs(text, query_duped_len=q_duped_len, debug=getattr(args, 'debug_ssw_parse', False))
        # Determine whether the SSW output used query/ref in the same order as our
        # TM-align convention (tm_pairs are written as i->j where i corresponds to pdb1 and j to pdb2).
        # parse_ssw_sam_pairs returns (query_index, ref_index). If the SAM/QNAME indicates
        # that query corresponds to pdb2 and ref corresponds to pdb1, we need to swap each pair
        # to produce (pdb1_index, pdb2_index) ordering for downstream comparison.
        try:
            first_non_header = None
            for ln in text.splitlines():
                if ln.strip() and not ln.startswith('@'):
                    first_non_header = ln.strip(); break
            if first_non_header:
                fields = first_non_header.split('\t')
                if len(fields) >= 3:
                    qname = fields[0]
                    rname = fields[2]
                    # crude stem-match: check if qname contains p2.stem and rname contains p1.stem
                    p1stem = p1.stem
                    p2stem = p2.stem
                    q_is_p2 = p2stem in qname or qname.startswith(p2stem)
                    r_is_p1 = p1stem in rname or rname.startswith(p1stem)
                    q_is_p1 = p1stem in qname or qname.startswith(p1stem)
                    r_is_p2 = p2stem in rname or rname.startswith(p2stem)
                    if q_is_p2 and r_is_p1 and not (q_is_p1 and r_is_p2):
                        # swap each parsed pair (query,ref) -> (ref,query)
                        pairs_3di = [(ref_i, query_i) for (query_i, ref_i) in pairs_3di]
                        logger.info('Detected SSW QNAME= p2 / RNAME= p1 ordering; swapped parsed SSW pairs to match tm (p1,p2) ordering')
        except Exception:
            logger.debug('Failed to inspect/align SSW QNAME/RNAME ordering; proceeding without swap', exc_info=True)
        logger.info(f'  {f.name}: parsed {len(pairs_3di)} 3di aligned pairs (after mod mapping)')

        # write ssw pairs to per-file TSV in tmpdir
        try:
            ssw_pairs_file = tmpdir / (f.name + '.pairs.tsv')
            with ssw_pairs_file.open('w') as sf:
                # pairs_3di are produced as (query_index, ref_index) to match TM (i,j) ordering,
                # so write header accordingly for clarity
                sf.write('#query\tref\n')
                for a,b in pairs_3di:
                    sf.write(f"{a}\t{b}\n")
            logger.info(f'Wrote SSW pairs to {ssw_pairs_file}')
        except Exception:
            logger.debug('Failed to write ssw pairs file', exc_info=True)

        # try to extract an SSW score from SAM or SSW textual output
        def extract_ssw_score(sam_text: str):
            # SAM optional tag AS:i:<score> preferred
            m = re.search(r'AS:i:(-?\d+)', sam_text)
            if m:
                try:
                    return float(m.group(1))
                except Exception:
                    pass
            # look for patterns like 'Score=123' or 'Score: 123' or 'score=123'
            m2 = re.search(r'[Ss]core\s*[=:]\s*(-?\d+)', sam_text)
            if m2:
                try:
                    return float(m2.group(1))
                except Exception:
                    pass
            # fallback: search for 'bit score' or 'bits' (rare)
            m3 = re.search(r'([-+]?\d+)\s+bits', sam_text)
            if m3:
                try:
                    return float(m3.group(1))
                except Exception:
                    pass
            return None

        ssw_score_val = extract_ssw_score(text)

        # prepare visual snippets for CSV: TM-align snippet and BLAST-like snippet if available
        tm_snippet = ''
        try:
            tm_snippet = ' | '.join([f"{a}-{b}" for a,b in tm_pairs[:20]])
        except Exception:
            tm_snippet = ''

        blast_snippet = ''
        # try to find converted blast output in same dir
        blast_path_candidates = []
        try:
            if args.alignment_results_dir:
                blast_path_candidates = list(Path(args.alignment_results_dir).glob(f"*{p1.stem}*{p2.stem}*blast_output.txt"))
        except Exception:
            blast_path_candidates = []
        if blast_path_candidates:
            try:
                btxt = blast_path_candidates[0].read_text()
                blast_snippet = '\n'.join(btxt.splitlines()[:6])
            except Exception:
                blast_snippet = ''

        # Correct SSW target indices to original B coordinates if needed, then rotate by cutpoint s
        pairs_3di_B = []
        try:
            max_jp = max((j for _,j in pairs_3di), default=0)
            for (ia, jp) in pairs_3di:
                # map to original B coordinate jb
                if getattr(args, 'using_bprime', False):
                    if getattr(args, 'window_start', None):
                        # if jp looks like absolute B' index (>lenT) treat as B' absolute
                        if max_jp > lenT:
                            jb = unwrap_bprime_to_b(jp, lenT)
                        else:
                            # assume jp is local window index 1..lenT
                            j_abs = args.window_start + jp - 1
                            jb = unwrap_bprime_to_b(j_abs, lenT)
                    else:
                        jb = unwrap_bprime_to_b(jp, lenT)
                elif getattr(args, 'target_slice', None):
                    jb = slice_to_b(jp, args.target_slice, lenT)
                else:
                    jb = jp
                pairs_3di_B.append((ia, jb))

            # Heuristic: if many SSW pairs appear to be systematically shifted relative to TM j indices,
            # try to detect a global offset and correct it. This handles cases where the SSW target
            # was a circular slice or otherwise offset relative to full-target indexing.
            try:
                # build mapping i -> tm_j (use first occurrence)
                tm_map = {i: j for i, j in tm_pairs}
                deltas = []
                for (ia, jb) in pairs_3di_B:
                    if ia in tm_map:
                        deltas.append(tm_map[ia] - jb)
                if deltas:
                    # Original behaviour: pick the mode (most common delta) and apply it
                    # if it has sufficient support. This is simpler and consistent with
                    # earlier pipeline logic that used the mode of deltas.
                    try:
                        from collections import Counter
                        cnt = Counter(deltas)
                        mode_delta, mode_count = cnt.most_common(1)[0]
                    except Exception:
                        mode_delta, mode_count = None, 0
                    min_support = max(2, int(0.15 * len(deltas)))
                    if mode_delta is not None and mode_count >= min_support:
                        logger.info(f'Detected global SSW->TM j-offset {mode_delta} (mode matches={mode_count}); applying correction')
                        def shift_j(jb, delta):
                            return ((jb + delta - 1) % lenT) + 1 if lenT else jb + delta
                        pairs_3di_B = [(ia, shift_j(jb, mode_delta)) for (ia, jb) in pairs_3di_B]
            except Exception:
                logger.debug('SSW->TM offset detection failed', exc_info=True)
        except Exception:
            pairs_3di_B = pairs_3di[:]

        # rotate SSW target indices
        pairs_3di_rot = [(ia, rotate_j(jb, s, lenT)) for (ia, jb) in pairs_3di_B]
        pairs_3di_rot_sorted = sorted(pairs_3di_rot)

        # TM rotated list (tm_rot) prepared earlier; ensure sorted
        tm_rot_sorted = sorted(tm_rot)

        # pairwise comparison with tolerance window w_pair
        metrics = compare_pairs(tm_rot_sorted, pairs_3di_rot_sorted, w=getattr(args, 'w_pair', 1))

        # range-based coverage metrics (pre-rotation ranges)
        T3_range, Q3_range = infer_ranges_from_pairs(pairs_3di_B)
        Ttm_range, Qtm_range = infer_ranges_from_pairs(tm_pairs)
        # split Ttm if wrapped
        tm_segs = split_wrap_interval(Ttm_range, lenT) if Ttm_range!=(0,0) else []
        overlap_len = circular_overlap_len(T3_range, tm_segs, lenT) if T3_range!=(0,0) else 0
        len_T3 = (T3_range[1]-T3_range[0]+1) if T3_range!=(0,0) else 0
        len_Ttm = sum([(b-a+1) for a,b in tm_segs]) if tm_segs else 0
        precision_T = overlap_len/len_T3 if len_T3>0 else 0.0
        recall_T = overlap_len/len_Ttm if len_Ttm>0 else 0.0
        f1_T = range_f1(precision_T, recall_T)

        # query (linear) side
        q_overlap = 0
        len_Q3 = (Q3_range[1]-Q3_range[0]+1) if Q3_range!=(0,0) else 0
        len_Qtm = (Qtm_range[1]-Qtm_range[0]+1) if Qtm_range!=(0,0) else 0
        if len_Q3>0 and len_Qtm>0:
            lo = max(Q3_range[0], Qtm_range[0])
            hi = min(Q3_range[1], Qtm_range[1])
            q_overlap = max(0, hi-lo+1)
        precision_Q = q_overlap/len_Q3 if len_Q3>0 else 0.0
        recall_Q = q_overlap/len_Qtm if len_Qtm>0 else 0.0
        f1_Q = range_f1(precision_Q, recall_Q)

        # CP sanity checks on rotated SSW pairs
        cp_wrap_ok = True
        try:
            decs = 0
            if pairs_3di_rot_sorted:
                prev_j = pairs_3di_rot_sorted[0][1]
                for _, jj in pairs_3di_rot_sorted[1:]:
                    if jj < prev_j:
                        decs += 1
                    prev_j = jj
            cp_wrap_ok = (decs <= 1)
        except Exception:
            cp_wrap_ok = False
        cp_span_ok = False
        if pairs_3di_rot_sorted:
            span = pairs_3di_rot_sorted[-1][1] - pairs_3di_rot_sorted[0][1] + 1
            cp_span_ok = (span <= lenT)
        cp_ok = cp_wrap_ok and cp_span_ok

        # contiguous run metrics
        lctr = lctr_longest_consistent_run(sorted(pairs_3di_rot_sorted))
        blockrec = block_recall_tau(sorted(tm_rot_sorted), sorted(pairs_3di_rot_sorted), getattr(args, 'tau', 8))

        # log detailed comparison info
        try:
            logger.info(f"Comparison for {f.name}: TM pairs={len(tm_rot_sorted)}, 3Di pairs={len(pairs_3di_rot_sorted)}, TP={metrics['TP']}, FP={metrics['FP']}, FN={metrics['FN']}")
            matched = metrics.get('matched_pairs', [])
            if matched:
                sample = matched[:10]
                sample_str = ','.join([f"tm_idx={t}->ssw_idx={s}" for t, s in sample])
                logger.info(f"  sample matched indices: {sample_str}")
            logger.info(f"CP flags: wrap_ok={cp_wrap_ok} span_ok={cp_span_ok} cp_ok={cp_ok} cutpoint={s} lenT={lenT}")
        except Exception:
            logger.debug('Failed to log detailed comparison info', exc_info=True)

        # include paths to pair files and snippets in CSV
        ssw_pairs_path = str(ssw_pairs_file) if 'ssw_pairs_file' in locals() else ''
        tm_pairs_path = str(tmpdir / 'tm_pairs.tsv') if (tmpdir / 'tm_pairs.tsv').exists() else ''
        # load TM-align snippet (first 200 lines) for CSV
        tm_align_snippet = ''
        try:
            snip = (tmpdir / 'TM_align_cp_result_snippet.txt')
            if snip.exists():
                tm_align_snippet = '\n'.join(snip.read_text().splitlines()[:20])
        except Exception:
            tm_align_snippet = ''

        row = {
            'ssw_file': str(f),
            'n_tm_pairs': len(tm_rot_sorted),
            'n_3di_pairs': len(pairs_3di_rot_sorted),
            'TP': metrics['TP'], 'FP': metrics['FP'], 'FN': metrics['FN'],
            'precision': metrics['precision'], 'recall': metrics['recall'], 'f1': metrics['f1'], 'jaccard': metrics['jaccard'],
            'range_target_precision': precision_T, 'range_target_recall': recall_T, 'range_target_f1': f1_T,
            'range_query_precision': precision_Q, 'range_query_recall': recall_Q, 'range_query_f1': f1_Q,
            # detailed range diagnostics (store explicit values for later inspection)
            'T3_range': T3_range, 'Ttm_range': Ttm_range, 'tm_segs': tm_segs,
            'overlap_T': overlap_len, 'len_T3': len_T3, 'len_Ttm': len_Ttm,
            'Q3_range': Q3_range, 'Qtm_range': Qtm_range, 'overlap_Q': q_overlap, 'len_Q3': len_Q3, 'len_Qtm': len_Qtm,
            'range_target_precision': precision_T, 'range_target_recall': recall_T, 'range_target_f1': f1_T,
            'range_query_precision': precision_Q, 'range_query_recall': recall_Q, 'range_query_f1': f1_Q,
            'block_f1_calc': 0.5 * (f1_T + f1_Q), 'stored_block_f1': None,
            'lctr': lctr, 'block_recall_tau': blockrec,
            # block_f1: compute from matched pairs when available (w-dependent)
            'block_f1': None,
            'matched_pairs': metrics.get('matched_pairs', []),
            'tm_matched_mask': metrics.get('tm_matched_mask', []),
            'three_matched_mask': metrics.get('three_matched_mask', []),
            'cp_wrap_ok': cp_wrap_ok, 'cp_span_ok': cp_span_ok, 'cp_ok': cp_ok,
            'cutpoint_s': s, 'len_target': lenT,
            'tm_snippet': tm_snippet, 'blast_snippet': blast_snippet,
            'tm_pairs_file': tm_pairs_path, 'ssw_pairs_file': ssw_pairs_path, 'tm_align_snippet': tm_align_snippet
        }
        # compute 3di sequence identity for this SSW alignment when possible
        try:
            # infer suffix (8f/9f/10f) from filename
            suf_m = re.search(r'_(8f|9f|10f)_', f.name)
            suffix_guess = suf_m.group(1) if suf_m else None
            if not suffix_guess:
                # fallback: find any variant key that appears in filename
                for k in new3di_by_var.keys():
                    if k in f.name:
                        suffix_guess = k
                        break
            three_q = ''
            three_t = ''
            if suffix_guess and suffix_guess in new3di_by_var:
                try:
                    qpath, tpath = new3di_by_var[suffix_guess]
                    three_q = read_fasta_sequence(Path(qpath))
                    three_t = read_fasta_sequence(Path(tpath))
                except Exception:
                    three_q = ''
                    three_t = ''
            # compute identity on pairs_3di_B (original coordinates) if sequences available
            identity_count = 0
            align_len = 0
            if three_q and three_t and pairs_3di_B:
                for (ia, jb) in pairs_3di_B:
                    if ia and jb and ia <= len(three_q) and jb <= len(three_t) and ia>0 and jb>0:
                        align_len += 1
                        if three_q[ia-1] == three_t[jb-1]:
                            identity_count += 1
            identity_val = (identity_count / align_len) if align_len>0 else None
            row['ssw_score'] = ssw_score_val
            row['3di_identity'] = identity_val
            row['3di_identity_match_count'] = identity_count
            row['3di_identity_alignment_len'] = align_len
        except Exception:
            row['ssw_score'] = ssw_score_val
            row['3di_identity'] = None
            row['3di_identity_match_count'] = 0
            row['3di_identity_alignment_len'] = 0
        rows.append(row)

    # NOTE: Per-user request, we do NOT write the overall CSV summary (tmalign_3di_match_summary.csv)
    # or the per-variant CSV. Instead we produce a single by-variant text file `tmalign_3di_match_summary_by_variant.txt`
    # (see below). Keep per-file rows in-memory for debugging if needed.

    # --- Variant-level summary: include TM-align visual, tm pair list, and per-variant 3di pairs + metrics ---
    variant_rows = []
    # load tm_align visual snippet (multi-line) if available
    tm_align_visual = ''
    try:
        snipf = tmpdir / 'TM_align_cp_result_snippet.txt'
        if snipf.exists():
            tm_align_visual = snipf.read_text()
    except Exception:
        tm_align_visual = ''

    # Enforce desired variant ordering: prefer 10f, 9f, 8f (append any other variants after)
    preferred_order = ['10f', '9f', '8f']
    ordered_suffixes = [s for s in preferred_order if s in new3di_by_var]
    ordered_suffixes += [s for s in sorted(new3di_by_var.keys()) if s not in set(preferred_order)]
    for suffix in ordered_suffixes:
        # try to find ssw file for this suffix in tmp_align_dir or align_dir
        ssw_path = None
        candidate_patterns = [f"*_{suffix}_*_ssw_output.txt", f"*_{suffix}_ssw_output.txt", f"*{suffix}*_ssw_output.txt"]
        for pat in candidate_patterns:
            found = list(tmp_align_dir.glob(pat))
            if not found:
                found = list(align_dir.glob(pat))
            if found:
                # prefer first matching file
                ssw_path = found[0]
                break
        ssw_pairs = []
        ssw_pairs_path = ''
        # Prefer to reuse the per-file parsed metrics we computed earlier (rows) to avoid
        # re-parsing SSW files and missing the global offset correction step.
        matching_rows = [r for r in rows if r.get('ssw_file') and ssw_path and os.path.basename(r.get('ssw_file')) == os.path.basename(ssw_path)]
        if matching_rows:
            # reuse parsed pairs and metrics from the first matching row
            r0 = matching_rows[0]
            try:
                ssw_pairs_path = r0.get('ssw_pairs_file', '')
            except Exception:
                ssw_pairs_path = ''
            # try to load pairs from the pairs.tsv we already wrote earlier
            if ssw_pairs_path and Path(ssw_pairs_path).exists():
                try:
                    for ln in Path(ssw_pairs_path).read_text().splitlines():
                        if not ln or ln.startswith('#'):
                            continue
                        a,b = ln.split('\t')
                        ssw_pairs.append((int(a), int(b)))
                except Exception:
                    ssw_pairs = []
            else:
                # fall back to parsing SSW file directly below (preserve old behaviour)
                pass
            # reuse numeric metrics from the parsed per-file row when available to ensure
            # variant-level summary matches per-file comparisons (avoids re-parsing discrepancies)
            try:
                metrics = {
                    'TP': int(r0.get('TP', 0)), 'FP': int(r0.get('FP', 0)), 'FN': int(r0.get('FN', 0)),
                    'precision': float(r0.get('precision', 0.0)), 'recall': float(r0.get('recall', 0.0)),
                    'f1': float(r0.get('f1', 0.0)), 'jaccard': float(r0.get('jaccard', 0.0)),
                }
            except Exception:
                metrics = None
        if ssw_path and not ssw_pairs:
            # prefer pairs TSV we wrote earlier
            candidate_pairs_tsv = tmp_align_dir / (ssw_path.name + '.pairs.tsv')
            if candidate_pairs_tsv.exists():
                try:
                    ssw_pairs_path = str(candidate_pairs_tsv)
                    for ln in candidate_pairs_tsv.read_text().splitlines():
                        if not ln or ln.startswith('#'):
                            continue
                        a,b = ln.split('\t')
                        ssw_pairs.append((int(a), int(b)))
                except Exception:
                    ssw_pairs = []
            else:
                # parse directly from SSW file
                try:
                    txt = ssw_path.read_text()
                    # attempt to infer query duped len as before (no CLI override available)
                    qlen = None
                    if qlen is None:
                        for ln in txt.splitlines():
                            if ln.strip() and not ln.startswith('@'):
                                flds = ln.split('\t')
                                if len(flds) >= 10:
                                    seq = flds[9]
                                    if seq:
                                        qlen = len(seq) // 2
                                break
                    ssw_pairs = parse_ssw_sam_pairs(txt, query_duped_len=qlen, debug=getattr(args, 'debug_ssw_parse', False))
                except Exception:
                    ssw_pairs = []

        # tm_pairs is global for this comparison
        try:
            tm_list_str = ';'.join([f"{i}-{j}" for i,j in tm_pairs])
        except Exception:
            tm_list_str = ''

        try:
            ssw_list_str = ';'.join([f"{i}-{j}" for i,j in ssw_pairs])
        except Exception:
            ssw_list_str = ''

        # Apply same CP-aware mapping and rotation for variant-level summary
        pairs_3di_B = []
        try:
            max_jp = max((j for _,j in ssw_pairs), default=0)
            for (ia, jp) in ssw_pairs:
                if getattr(args, 'using_bprime', False):
                    if getattr(args, 'window_start', None):
                        if max_jp > lenT:
                            jb = unwrap_bprime_to_b(jp, lenT)
                        else:
                            j_abs = args.window_start + jp - 1
                            jb = unwrap_bprime_to_b(j_abs, lenT)
                    else:
                        jb = unwrap_bprime_to_b(jp, lenT)
                elif getattr(args, 'target_slice', None):
                    jb = slice_to_b(jp, args.target_slice, lenT)
                else:
                    jb = jp
                pairs_3di_B.append((ia, jb))
        except Exception:
            pairs_3di_B = ssw_pairs[:]

        # Heuristic offset correction (same logic as used earlier per-file):
        # detect a common delta between TM j and SSW jb and apply mode-based correction
        try:
            tm_map = {i: j for i, j in tm_pairs}
            deltas = []
            for (ia, jb) in pairs_3di_B:
                if ia in tm_map:
                    deltas.append(tm_map[ia] - jb)
            if deltas:
                from collections import Counter
                cnt = Counter(deltas)
                mode_delta, mode_count = cnt.most_common(1)[0]
                min_support = max(2, int(0.15 * len(deltas)))
                if mode_delta is not None and mode_count >= min_support:
                    logger.info(f'Variant-level: detected global SSW->TM j-offset {mode_delta} (mode matches={mode_count}); applying correction')
                    def shift_j_var(jb, delta):
                        return ((jb + delta - 1) % lenT) + 1 if lenT else jb + delta
                    pairs_3di_B = [(ia, shift_j_var(jb, mode_delta)) for (ia, jb) in pairs_3di_B]
        except Exception:
            logger.debug('Variant-level SSW->TM offset detection failed', exc_info=True)

        pairs_3di_rot = [(ia, rotate_j(jb, s, lenT)) for (ia, jb) in pairs_3di_B]
        pairs_3di_rot_sorted = sorted(pairs_3di_rot)
        tm_rot_sorted = sorted(tm_rot)
        # If metrics were populated from earlier per-file parsing, do not recompute here
        if not ('metrics' in locals() and metrics is not None):
            metrics = compare_pairs(tm_rot_sorted, pairs_3di_rot_sorted, w=getattr(args, 'w_pair', 1))

        # range-based for variant (original, non-matched ranges)
        T3_range, Q3_range = infer_ranges_from_pairs(pairs_3di_B)
        Ttm_range, Qtm_range = infer_ranges_from_pairs(tm_pairs)
        tm_segs = split_wrap_interval(Ttm_range, lenT) if Ttm_range!=(0,0) else []
        overlap_len = circular_overlap_len(T3_range, tm_segs, lenT) if T3_range!=(0,0) else 0
        len_T3 = (T3_range[1]-T3_range[0]+1) if T3_range!=(0,0) else 0
        len_Ttm = sum([(b-a+1) for a,b in tm_segs]) if tm_segs else 0
        precision_T = overlap_len/len_T3 if len_T3>0 else 0.0
        recall_T = overlap_len/len_Ttm if len_Ttm>0 else 0.0
        f1_T = range_f1(precision_T, recall_T)

        # query
        q_overlap = 0
        len_Q3 = (Q3_range[1]-Q3_range[0]+1) if Q3_range!=(0,0) else 0
        len_Qtm = (Qtm_range[1]-Qtm_range[0]+1) if Qtm_range!=(0,0) else 0
        if len_Q3>0 and len_Qtm>0:
            lo = max(Q3_range[0], Qtm_range[0]); hi = min(Q3_range[1], Qtm_range[1])
            q_overlap = max(0, hi-lo+1)
        precision_Q = q_overlap/len_Q3 if len_Q3>0 else 0.0
        recall_Q = q_overlap/len_Qtm if len_Qtm>0 else 0.0
        f1_Q = range_f1(precision_Q, recall_Q)

        decs = 0
        if pairs_3di_rot_sorted:
            prev_j = pairs_3di_rot_sorted[0][1]
            for _, jj in pairs_3di_rot_sorted[1:]:
                if jj < prev_j:
                    decs += 1
                prev_j = jj
        cp_wrap_ok = (decs <= 1)
        cp_span_ok = False
        if pairs_3di_rot_sorted:
            span = pairs_3di_rot_sorted[-1][1] - pairs_3di_rot_sorted[0][1] + 1
            cp_span_ok = (span <= lenT)
        cp_ok = cp_wrap_ok and cp_span_ok
        lctr = lctr_longest_consistent_run(sorted(pairs_3di_rot_sorted))
        blockrec = block_recall_tau(sorted(tm_rot_sorted), sorted(pairs_3di_rot_sorted), getattr(args, 'tau', 8))

        # Produce explicit compact 3-row coverage files for human inspection
        try:
            # create masks of length lenT
            if not isinstance(lenT, int) or lenT <= 0:
                raise ValueError('lenT must be a positive integer to render coverage files')
            tm_cov = [False] * lenT
            ssw_cov = [False] * lenT
            # mark TM coverage (use rotated TM indices)
            for _, j in tm_rot:
                if 1 <= j <= lenT:
                    tm_cov[j-1] = True
            # mark SSW coverage (use rotated SSW indices)
            for _, j in pairs_3di_rot:
                if 1 <= j <= lenT:
                    ssw_cov[j-1] = True

            # load sequences when available
            try:
                aa_seq = read_fasta_sequence(fasta2) if 'fasta2' in locals() else ''
            except Exception:
                aa_seq = ''
            three_seq = ''
            try:
                three_path = new3di_by_var[suffix][1]
                if three_path.exists():
                    three_seq = read_fasta_sequence(three_path)
            except Exception:
                three_seq = ''

            # build compact 3-row per-target view
            indices = list(range(1, lenT + 1))
            # TM chars: show AA letter where TM covers, else space
            tm_chars = [aa_seq[j-1] if (j-1) < len(aa_seq) and tm_cov[j-1] else ' ' for j in indices]
            # SSW chars: show 3Di token where SSW covers, else space
            ssw_chars = [three_seq[j-1] if (j-1) < len(three_seq) and ssw_cov[j-1] else ' ' for j in indices]

            cov_t_path = tmp_align_dir / f"{Path(p2).stem}_{suffix}_result_coverage_target.tsv"
            cov_t_path.write_text(render_three_row_coverage(indices, tm_chars, ssw_chars, label_tm='TM', label_ssw='SSW'))
            logger.info(f'Wrote target coverage view to {cov_t_path}')

            # build per-query compact 3-row view
            # infer lenQ and sequences for query
            seqA = read_fasta_sequence(fasta1) if 'fasta1' in locals() else (read_fasta_sequence(fasta1) if fasta1 else '')
            threeA = ''
            try:
                threeA_path = new3di_by_var[suffix][0]
                if Path(threeA_path).exists():
                    threeA = read_fasta_sequence(Path(threeA_path))
            except Exception:
                threeA = ''
            lenQ_eff = lenQ if isinstance(lenQ, int) and lenQ>0 else len(seqA)
            q_indices = list(range(1, lenQ_eff + 1))
            # show actual query sequence residues where TM covers (not placeholders)
            tm_q_chars = [seqA[q-1] if (q-1) < len(seqA) and any(ti==q for ti,_ in tm_rot) else ' ' for q in q_indices]
            ssw_q_chars = [threeA[q-1] if q-1 < len(threeA) and any(q==ia for ia,_ in pairs_3di_rot_sorted) else ' ' for q in q_indices]
            cov_q_path = tmp_align_dir / f"{Path(p2).stem}_{suffix}_result_coverage_query.tsv"
            cov_q_path.write_text(render_three_row_coverage(q_indices, tm_q_chars, ssw_q_chars, label_tm='TM', label_ssw='SSW'))
            logger.info(f'Wrote query coverage view to {cov_q_path}')

            # --- Residue-level confusion for TARGET ---
            TP_T = FP_T = FN_T = TN_T = 0
            for j in range(lenT):
                tm = tm_cov[j]
                ssw = ssw_cov[j]
                if tm and ssw:
                    TP_T += 1
                elif tm and not ssw:
                    FN_T += 1
                elif (not tm) and ssw:
                    FP_T += 1
                else:
                    TN_T += 1

            sens_T_block = TP_T / (TP_T + FN_T) if (TP_T + FN_T) > 0 else 1.0
            spec_T_block = TN_T / (TN_T + FP_T) if (TN_T + FP_T) > 0 else 1.0

            # --- Residue-level confusion for QUERY ---
            lenQ_eff = lenQ if isinstance(lenQ, int) and lenQ > 0 else lenQ
            if not isinstance(lenQ_eff, int) or lenQ_eff <= 0:
                lenQ_eff = max((i for i,_ in tm_rot), default=0)
                lenQ_eff = max(lenQ_eff, max((i for i,_ in pairs_3di_rot_sorted), default=0))

            tm_cov_Q = [False] * lenQ_eff
            ssw_cov_Q = [False] * lenQ_eff

            for i, _ in tm_rot:
                if 1 <= i <= lenQ_eff:
                    tm_cov_Q[i-1] = True

            for i, _ in pairs_3di_rot_sorted:
                if 1 <= i <= lenQ_eff:
                    ssw_cov_Q[i-1] = True

            TP_Q = FP_Q = FN_Q = TN_Q = 0
            for i in range(lenQ_eff):
                tm = tm_cov_Q[i]
                ssw = ssw_cov_Q[i]
                if tm and ssw:
                    TP_Q += 1
                elif tm and not ssw:
                    FN_Q += 1
                elif (not tm) and ssw:
                    FP_Q += 1
                else:
                    TN_Q += 1

            sens_Q_block = TP_Q / (TP_Q + FN_Q) if (TP_Q + FN_Q) > 0 else 1.0
            spec_Q_block = TN_Q / (TN_Q + FP_Q) if (TN_Q + FP_Q) > 0 else 1.0

            # record paths
            ssw_pairs_path = ssw_pairs_path if 'ssw_pairs_path' in locals() else ''
            cov_t_str = str(cov_t_path)
            cov_q_str = str(cov_q_path)
        except Exception:
            logger.exception(f'Failed to produce coverage-by-target/query TSVs for variant {suffix}')
            # set defaults if computation failed
            TP_T = FP_T = FN_T = TN_T = 0
            sens_T_block = 1.0
            spec_T_block = 1.0
            TP_Q = FP_Q = FN_Q = TN_Q = 0
            sens_Q_block = 1.0
            spec_Q_block = 1.0

        # No UNION / residue-wise generation: we now produce explicit coverage TSVs instead.
        union_text = ''
        union_file_str = ''

        # Compute a w-dependent block-level F1 from matched pairs (matched indices are from
        # compare_pairs which was run on rotated lists tm_rot_sorted and pairs_3di_rot_sorted)
        matched = metrics.get('matched_pairs', []) if metrics else []
        block_f1_val = None
        try:
            if matched:
                matched_tm = [tm_rot_sorted[t_idx] for (t_idx, _) in matched]
                matched_3di = [pairs_3di_rot_sorted[three_idx] for (_, three_idx) in matched]
                # infer ranges from matched sets (rotated coordinates)
                Tm_mat_range, Qm_mat_range = infer_ranges_from_pairs(matched_3di)
                Ttm_mat_range, Qtm_mat_range = infer_ranges_from_pairs(matched_tm)
                # compute target overlap
                tm_mat_segs = split_wrap_interval(Ttm_mat_range, lenT) if Ttm_mat_range!=(0,0) else []
                overlap_len_m = circular_overlap_len(Tm_mat_range, tm_mat_segs, lenT) if Tm_mat_range!=(0,0) else 0
                len_Tm3 = (Tm_mat_range[1]-Tm_mat_range[0]+1) if Tm_mat_range!=(0,0) else 0
                len_Tmtm = sum([(b-a+1) for a,b in tm_mat_segs]) if tm_mat_segs else 0
                precision_Tm = overlap_len_m/len_Tm3 if len_Tm3>0 else 0.0
                recall_Tm = overlap_len_m/len_Tmtm if len_Tmtm>0 else 0.0
                f1_Tm = range_f1(precision_Tm, recall_Tm)
                # query side
                q_overlap_m = 0
                len_Qm3 = (Qm_mat_range[1]-Qm_mat_range[0]+1) if Qm_mat_range!=(0,0) else 0
                len_Qmtm = (Qtm_mat_range[1]-Qtm_mat_range[0]+1) if Qtm_mat_range!=(0,0) else 0
                if len_Qm3>0 and len_Qmtm>0:
                    lo = max(Qm_mat_range[0], Qtm_mat_range[0]); hi = min(Qm_mat_range[1], Qtm_mat_range[1])
                    q_overlap_m = max(0, hi-lo+1)
                precision_Qm = q_overlap_m/len_Qm3 if len_Qm3>0 else 0.0
                recall_Qm = q_overlap_m/len_Qmtm if len_Qmtm>0 else 0.0
                f1_Qm = range_f1(precision_Qm, recall_Qm)
                block_f1_val = 0.5 * (f1_Tm + f1_Qm)
            else:
                block_f1_val = 0.0
        except Exception:
            block_f1_val = None

        # find any per-ssw row we computed earlier for this variant to extract ssw_score and identity
        ssw_score_for_variant = None
        identity_for_variant = None
        identity_match_count_for_variant = 0
        identity_align_len_for_variant = 0
        try:
            matching_rows = [r for r in rows if r.get('ssw_file') and (suffix in r.get('ssw_file'))]
            if matching_rows:
                r0 = matching_rows[0]
                ssw_score_for_variant = r0.get('ssw_score') if r0.get('ssw_score') is not None else None
                identity_for_variant = r0.get('3di_identity') if r0.get('3di_identity') is not None else None
                identity_match_count_for_variant = r0.get('3di_identity_match_count', 0)
                identity_align_len_for_variant = r0.get('3di_identity_alignment_len', 0)
        except Exception:
            ssw_score_for_variant = None

        # Finalize block_f1: prefer matched-based block_f1_val; if missing or zero, fall back
        # to the simpler range-based average of range_target_f1 and range_query_f1.
        # Compute block-level F1 using range-based method (always):
        try:
            block_f1_final = 0.5 * (float(f1_T) if f1_T is not None else 0.0) + 0.5 * (float(f1_Q) if f1_Q is not None else 0.0)
        except Exception:
            block_f1_final = 0.0

        variant_rows.append({
            'variant': suffix,
            'pair': f"{p1.stem}_vs_{p2.stem}",
            'query_name': p1.stem,
            'target_name': p2.stem,
            'w': getattr(args, 'w_pair', 1),
            'tm_align_visual': tm_align_visual,
            'tm_pairs_count': len(tm_pairs),
            'tm_pairs_list': tm_list_str,
            'ssw_pairs_count': len(ssw_pairs),
            'TP': metrics['TP'], 'FP': metrics['FP'], 'FN': metrics['FN'],
            'precision': metrics['precision'], 'recall': metrics['recall'], 'f1': metrics['f1'], 'jaccard': metrics['jaccard'],
            # residue-level confusion for target
            'TP_T': TP_T if 'TP_T' in locals() else 0,
            'FP_T': FP_T if 'FP_T' in locals() else 0,
            'FN_T': FN_T if 'FN_T' in locals() else 0,
            'TN_T': TN_T if 'TN_T' in locals() else 0,
            'sens_T_block': sens_T_block if 'sens_T_block' in locals() else 1.0,
            'spec_T_block': spec_T_block if 'spec_T_block' in locals() else 1.0,
            # residue-level confusion for query
            'TP_Q': TP_Q if 'TP_Q' in locals() else 0,
            'FP_Q': FP_Q if 'FP_Q' in locals() else 0,
            'FN_Q': FN_Q if 'FN_Q' in locals() else 0,
            'TN_Q': TN_Q if 'TN_Q' in locals() else 0,
            'sens_Q_block': sens_Q_block if 'sens_Q_block' in locals() else 1.0,
            'spec_Q_block': spec_Q_block if 'spec_Q_block' in locals() else 1.0,
            'range_target_precision': precision_T, 'range_target_recall': recall_T, 'range_target_f1': f1_T,
            'range_query_precision': precision_Q, 'range_query_recall': recall_Q, 'range_query_f1': f1_Q,
            # detailed range diagnostics
            'T3_range': T3_range, 'Ttm_range': Ttm_range, 'tm_segs': tm_segs,
            'overlap_T': overlap_len, 'len_T3': len_T3, 'len_Ttm': len_Ttm,
            'Q3_range': Q3_range, 'Qtm_range': Qtm_range, 'overlap_Q': q_overlap, 'len_Q3': len_Q3, 'len_Qtm': len_Qtm,
            'precision_T': precision_T, 'recall_T': recall_T, 'range_target_f1': f1_T,
            'precision_Q': precision_Q, 'recall_Q': recall_Q, 'range_query_f1': f1_Q,
            'block_f1_calc': 0.5 * (f1_T + f1_Q),
            'lctr': lctr, 'block_recall_tau': blockrec,
            'cp_wrap_ok': cp_wrap_ok, 'cp_span_ok': cp_span_ok, 'cp_ok': cp_ok,
            'cutpoint_s': s, 'len_target': lenT,
            'ssw_file': str(ssw_path) if ssw_path else '', 'ssw_pairs_file': ssw_pairs_path,
            'union_file': union_file_str, 'union_text': union_text,
            'coverage_by_target': cov_t_str if 'cov_t_str' in locals() else '',
            'coverage_by_query': cov_q_str if 'cov_q_str' in locals() else '',
            # block_f1 is intentionally omitted here; use range_target_f1/range_query_f1 for downstream plots
            'ssw_score': ssw_score_for_variant,
            '3di_identity': identity_for_variant,
            '3di_identity_match_count': identity_match_count_for_variant,
            '3di_identity_alignment_len': identity_align_len_for_variant
        })

    # Also produce a horizontal (wide) summary file similar to older outputs for quick glance
    try:
        horiz_path = tmpdir / f"tmalign_3di_match_summary_by_variant_w{getattr(args, 'w_pair', 1)}.txt"
        cols = ['variant', 'tm_pairs_count', 'ssw_pairs_count', 'TP', 'FP', 'FN', 'precision', 'recall', 'f1', 'jaccard', 'range_target_f1', 'range_query_f1', 'lctr', 'block_recall_tau', 'cp_ok', 'cutpoint_s', 'len_target']
        # header
        lines = []
        lines.append('Variant summary:')
        # table header line
        hdr = ' '.join([c.rjust(10) for c in cols])
        lines.append(hdr)
        lines.append('-' * max(80, len(hdr)))
        # each variant as a row
        for vr in variant_rows:
            row_vals = []
            for c in cols:
                v = vr.get(c, '')
                if v is None:
                    s = ''
                elif isinstance(v, float):
                    s = f"{v:.6g}"
                else:
                    s = str(v)
                row_vals.append(s.rjust(10))
            lines.append(' '.join(row_vals))

        # append coverage blocks per variant (target then query)
        lines.append('\n')
        for vr in variant_rows:
            suf = vr.get('variant', '')
            lines.append(f"============================================================\nCoverage by target ({suf}):")
            cov_t = vr.get('coverage_by_target', '')
            if cov_t and Path(cov_t).exists():
                try:
                    lines.append(Path(cov_t).read_text())
                except Exception:
                    lines.append(f"(coverage file exists but could not be read: {cov_t})")
            else:
                lines.append('(no coverage target file)')
            lines.append(f'\nCoverage by query ({suf}):')
            cov_q = vr.get('coverage_by_query', '')
            if cov_q and Path(cov_q).exists():
                try:
                    lines.append(Path(cov_q).read_text())
                except Exception:
                    lines.append(f"(coverage file exists but could not be read: {cov_q})")
            else:
                lines.append('(no coverage query file)')
            lines.append('\n')

        horiz_path.write_text('\n'.join(lines))
        logger.info(f'Wrote horizontal variant summary to {horiz_path}')
    except Exception:
        logger.debug('Failed to write horizontal summary file', exc_info=True)



    # Produce a CSV variant summary suitable for plotting (one row per variant)
    try:
        wval = getattr(args, 'w_pair', 1)
        csv_path = tmpdir / f"{p1.stem}_vs_{p2.stem}_variant_summary_w{wval}.csv"
        # choose a stable column order for CSV (include ssw_score & 3di identity fields)
        fields = [
            'pair', 'query_name', 'target_name', 'variant', 'w',
            'tm_pairs_count', 'ssw_pairs_count', 'TP', 'FP', 'FN', 'precision', 'recall', 'f1', 'jaccard',
            # residue-level confusion for target
            'TP_T', 'FP_T', 'FN_T', 'TN_T', 'sens_T_block', 'spec_T_block',
            # residue-level confusion for query
            'TP_Q', 'FP_Q', 'FN_Q', 'TN_Q', 'sens_Q_block', 'spec_Q_block',
            # detailed range diagnostics
            'T3_range', 'Ttm_range', 'tm_segs', 'overlap_T', 'len_T3', 'len_Ttm', 'precision_T', 'recall_T', 'range_target_f1',
            'Q3_range', 'Qtm_range', 'overlap_Q', 'len_Q3', 'len_Qtm', 'precision_Q', 'recall_Q', 'range_query_f1',
            'lctr', 'block_recall_tau', 'block_f1_calc',
            'ssw_score', '3di_identity', '3di_identity_match_count', '3di_identity_alignment_len', 'cp_ok', 'cutpoint_s', 'len_target',
            'tm_pairs_list', 'ssw_pairs_file', 'coverage_by_target', 'coverage_by_query'
        ]
        with csv_path.open('w', newline='') as cf:
            writer = csv.DictWriter(cf, fieldnames=fields)
            writer.writeheader()
            for vr in variant_rows:
                # ensure only stringable scalar values are written
                out = {k: vr.get(k, '') for k in fields}
                writer.writerow(out)
        print('Wrote variant summary CSV to', str(csv_path))
    except Exception:
        logger.exception('Failed to write variant summary CSV')

    # Also write a verbose human-readable text summary for quick inspection
    try:
        pretty_path = tmpdir / f"{p1.stem}_vs_{p2.stem}_variant_summary_w{wval}.txt"
        lines = []
        lines.append(f"Pair: {p1.stem} vs {p2.stem}")
        lines.append(f"Run dir: {tmpdir}")
        lines.append(f"w: {wval}")
        lines.append('')
        lines.append('TM pairs (first 200 shown):')
        try:
            lines.append('\n'.join([f"{i}\t{j}" for i,j in (tm_pairs[:200] if isinstance(tm_pairs, list) else [])]))
        except Exception:
            pass
        lines.append('')
        lines.append('Per-variant summaries:')
        for vr in variant_rows:
            lines.append('------------------------------------------------------------')
            for k in sorted(vr.keys()):
                v = vr.get(k)
                # truncate long lists/strings for readability but include full ones separately
                if isinstance(v, (list, tuple)):
                    lines.append(f"{k}: {v}")
                else:
                    lines.append(f"{k}: {v}")
            lines.append('')
        pretty_path.write_text('\n'.join(lines))
        print('Wrote verbose variant summary to', str(pretty_path))
    except Exception:
        logger.exception('Failed to write verbose pretty summary text')

    # (pretty printing removed per user request)


if __name__ == '__main__':
    main()
