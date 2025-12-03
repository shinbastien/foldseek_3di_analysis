#!/usr/bin/env python3
from pathlib import Path
import csv
import sys
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
TMP = ROOT / 'tmp'
CSV_GLOB = list(TMP.glob('*/*_variant_summary*.csv'))

if not CSV_GLOB:
    print('No per-run variant_summary CSVs under tmp/. Run pipeline first.')
    sys.exit(2)


def read_tm_pairs(tm_pairs_path):
    pairs = []
    if not tm_pairs_path.exists():
        return pairs
    with tm_pairs_path.open() as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            a,b = ln.split()[:2]
            pairs.append((int(a), int(b)))
    return pairs


def read_ssw_pairs(pairs_path):
    pairs = []
    p = Path(pairs_path)
    if not p.exists():
        return pairs
    with p.open() as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            toks = ln.split('\t') if '\t' in ln else ln.split()
            try:
                a = int(toks[0]); b = int(toks[1])
                pairs.append((a,b))
            except Exception:
                continue
    return pairs


def compare_pairs(tm_pairs, three_pairs, w=0):
    # Use ref-sorted sliding-window matching to avoid ordering/shift issues
    tm_matched = [False]*len(tm_pairs)
    three_matched = [False]*len(three_pairs)
    matched_pairs = []
    tm_indexed = [(t_idx, it, jt) for t_idx, (it, jt) in enumerate(tm_pairs)]
    three_indexed = [(s_idx, i3, j3) for s_idx, (i3, j3) in enumerate(three_pairs)]
    tm_indexed.sort(key=lambda x: x[2])
    three_indexed.sort(key=lambda x: x[2])
    left = 0
    n_tm = len(tm_indexed)
    for s_idx, i3, j3 in three_indexed:
        while left < n_tm and tm_indexed[left][2] < j3 - w:
            left += 1
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
    FN = max(0, len(tm_pairs)-TP)
    FP = max(0, len(three_pairs)-matched_three)
    precision = TP/(TP+FP) if (TP+FP)>0 else 0.0
    recall = TP/(TP+FN) if (TP+FN)>0 else 0.0
    f1 = (2*precision*recall/(precision+recall)) if (precision+recall)>0 else 0.0
    return {'TP':TP,'FP':FP,'FN':FN,'precision':precision,'recall':recall,'f1':f1,'matched_pairs':matched_pairs}


def split_wrap_interval(seg, lenB):
    l,r = seg
    if l<=r:
        return [(l,r)]
    return [(l,lenB),(1,r)]


def overlap_len(seg, segs):
    # seg is (l,r), segs is list of (l,r), compute total overlap length
    l,r = seg
    total=0
    for bl,br in segs:
        lo = max(l,bl)
        hi = min(r,br)
        if hi>=lo:
            total += (hi-lo+1)
    return total


def infer_range(pairs):
    if not pairs:
        return (0,0)
    is_, js_ = zip(*pairs)
    return (min(js_), max(js_)), (min(is_), max(is_))


reports = []
for csvp in sorted(CSV_GLOB):
    try:
        import pandas as pd
        df = pd.read_csv(csvp)
    except Exception:
        # fallback to csv module
        with open(csvp) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        df = None
    # iterate rows
    rows = df.to_dict('records') if df is not None else rows
    for r in rows:
        try:
            block_f1 = float(r.get('block_f1', 0.0))
        except Exception:
            block_f1 = 0.0
        if block_f1 > 1e-8:
            continue
        pair = r.get('pair')
        variant = r.get('variant')
        w = int(r.get('w', 1))
        tm_pairs_file = None
        # derive run dir from csvp path
        run_dir = csvp.parent
        # try standard tm_pairs.tsv
        cand_tm = run_dir / 'tm_pairs.tsv'
        if not cand_tm.exists():
            # also try parent dir
            cand_tm = run_dir.parent / 'tm_pairs.tsv'
        tm_pairs = read_tm_pairs(cand_tm)
        ssw_pairs_file = r.get('ssw_pairs_file','')
        ssw_pairs = []
        if ssw_pairs_file and Path(ssw_pairs_file).exists():
            ssw_pairs = read_ssw_pairs(Path(ssw_pairs_file))
            ssw_output_path = Path(ssw_pairs_file).with_suffix('')
            # try to locate original ssw output by removing .pairs.tsv
            ssw_output_candidate = Path(ssw_pairs_file).with_suffix('')
        else:
            # try alignment_results for file matching variant
            ar = run_dir / 'alignment_results'
            ssw_files = list(ar.glob(f'*_{variant}_ssw_output.txt')) if ar.exists() else []
            if ssw_files:
                ssw_output_candidate = ssw_files[0]
                # try to parse pairs from it (not necessary but try)
                # fallback: try to find pairs.tsv created earlier
                pss = Path(str(ssw_output_candidate) + '.pairs.tsv')
                if pss.exists():
                    ssw_pairs = read_ssw_pairs(pss)
                    ssw_pairs_file = str(pss)
        # diagnostics
        print('\n==== DIAG:', pair, 'variant=', variant, 'w=', w, '====')
        print('tm_pairs_file:', cand_tm if cand_tm.exists() else '(missing)')
        print('n_tm_pairs:', len(tm_pairs))
        print('ssw_pairs_file:', ssw_pairs_file)
        print('n_ssw_pairs:', len(ssw_pairs))
        # matched pairs sample
        comp = compare_pairs(tm_pairs, ssw_pairs, w=w)
        matched = comp.get('matched_pairs', [])
        print('Computed match metrics (recomputed): TP=%d FP=%d FN=%d precision=%.3g recall=%.3g f1=%.3g' % (comp['TP'], comp['FP'], comp['FN'], comp['precision'], comp['recall'], comp['f1']))
        print('Matched pairs sample (up to 20):')
        for t_idx,s_idx in matched[:20]:
            try:
                print(' tm_pair:', tm_pairs[t_idx], 'ssw_pair:', ssw_pairs[s_idx])
            except Exception:
                pass
        # compute range-based F1s per user's requested formula
        lenT = int(float(r.get('len_target', 0))) if r.get('len_target') not in (None,'') else (max((j for _,j in tm_pairs), default=0))
        # infer ranges
        T3_range, Q3_range = infer_range(ssw_pairs)
        Ttm_range, Qtm_range = infer_range(tm_pairs)
        tm_segs = split_wrap_interval(Ttm_range, lenT) if Ttm_range!=(0,0) else []
        overlap_T = overlap_len(T3_range, tm_segs)
        len_T3 = (T3_range[1]-T3_range[0]+1) if T3_range!=(0,0) else 0
        len_Ttm = sum((b-a+1) for a,b in tm_segs) if tm_segs else 0
        precision_T = overlap_T/len_T3 if len_T3>0 else 0.0
        recall_T = overlap_T/len_Ttm if len_Ttm>0 else 0.0
        def f1(p, r):
            return (2*p*r/(p+r)) if (p+r)>0 else 0.0
        range_target_f1 = f1(precision_T, recall_T)
        # query side
        Q3 = Q3_range; Qtm = Qtm_range
        len_Q3 = (Q3[1]-Q3[0]+1) if Q3!=(0,0) else 0
        len_Qtm = (Qtm[1]-Qtm[0]+1) if Qtm!=(0,0) else 0
        overlap_Q = 0
        if len_Q3>0 and len_Qtm>0:
            lo = max(Q3[0], Qtm[0]); hi = min(Q3[1], Qtm[1]); overlap_Q = max(0, hi-lo+1)
        precision_Q = overlap_Q/len_Q3 if len_Q3>0 else 0.0
        recall_Q = overlap_Q/len_Qtm if len_Qtm>0 else 0.0
        range_query_f1 = f1(precision_Q, recall_Q)
        block_f1_calc = 0.5*(range_target_f1 + range_query_f1)
        print('\nRange-based diagnostics (user formula):')
        print(' T3_range=', T3_range, ' Ttm_range=', Ttm_range, ' tm_segs=', tm_segs)
        print(' overlap_T=', overlap_T, ' len_T3=', len_T3, ' len_Ttm=', len_Ttm)
        print(' precision_T=%.6g recall_T=%.6g range_target_f1=%.6g' % (precision_T, recall_T, range_target_f1))
        print(' Q3_range=', Q3_range, ' Qtm_range=', Qtm_range)
        print(' overlap_Q=', overlap_Q, ' len_Q3=', len_Q3, ' len_Qtm=', len_Qtm)
        print(' precision_Q=%.6g recall_Q=%.6g range_query_f1=%.6g' % (precision_Q, recall_Q, range_query_f1))
        print(' block_f1 (calc)=%.6g stored_block_f1=%s' % (block_f1_calc, str(r.get('block_f1'))))

        # compute SSW->TM delta mode and counts
        deltas = []
        tm_map = {i:j for i,j in tm_pairs}
        for (ia,jb) in ssw_pairs:
            if ia in tm_map:
                deltas.append(tm_map[ia] - jb)
        if deltas:
            cnt = Counter(deltas)
            mode_delta, mode_count = cnt.most_common(1)[0]
            min_support = max(2, int(0.15*len(deltas)))
            print('\nSSW->TM delta stats: n=%d modes:{%s:%d} min_support=%d' % (len(deltas), mode_delta, mode_count, min_support))
            if mode_count>=min_support:
                print(' mode_delta selected ->', mode_delta)
            else:
                print(' no dominant mode (mode_count < min_support)')
        else:
            print('\nNo deltas computed (no matching i indices between SSW and TM)')

        # inspect ssw output first non-header line for QNAME/RNAME
        try:
            # try to find original ssw output file
            cand = None
            if ssw_pairs_file:
                cand = Path(ssw_pairs_file).with_suffix('')
                if not cand.exists():
                    # try replacing .pairs.tsv with _ssw_output.txt
                    cand = Path(ssw_pairs_file.replace('.pairs.tsv',''))
            else:
                ar = run_dir / 'alignment_results'
                if ar.exists():
                    L = list(ar.glob('*_ssw_output.txt'))
                    cand = L[0] if L else None
            if cand and cand.exists():
                first=None
                with cand.open() as fh:
                    for ln in fh:
                        if ln.strip() and not ln.startswith('@'):
                            first=ln.strip(); break
                if first:
                    flds = first.split('\t')
                    qname = flds[0] if len(flds)>0 else ''
                    rname = flds[2] if len(flds)>2 else ''
                    print('\nSSW first alignment line: qname=%s rname=%s (file=%s)' % (qname, rname, cand))
                    p1stem = r.get('query_name'); p2stem = r.get('target_name')
                    q_is_p2 = p2stem in qname or qname.startswith(p2stem)
                    r_is_p1 = p1stem in rname or rname.startswith(p1stem)
                    q_is_p1 = p1stem in qname or qname.startswith(p1stem)
                    r_is_p2 = p2stem in rname or rname.startswith(p2stem)
                    if q_is_p2 and r_is_p1 and not (q_is_p1 and r_is_p2):
                        print('Detected SSW QNAME=p2 / RNAME=p1 ordering (SSW used p2 as query and p1 as ref)')
                    else:
                        print('No QNAME/RNAME swap detected (SSW ordering appears standard or ambiguous)')
                else:
                    print('\nCould not find non-header alignment line in ssw output file', cand)
            else:
                print('\nOriginal SSW output file not found for deeper inspection')
        except Exception as e:
            print('SSW file inspection failed:', e)

print('\nDiagnostics complete')
