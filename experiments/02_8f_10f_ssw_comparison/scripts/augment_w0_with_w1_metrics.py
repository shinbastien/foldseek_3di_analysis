#!/usr/bin/env python3
"""Augment per-run variant_summary_w0.csv files with metrics copied from the matching _w1 CSV.

For each tmp/*/*_variant_summary_w0.csv, find tmp/<same-run>/*_variant_summary_w1.csv and copy TP/FP/FN/precision/recall/f1
for each variant into new columns TP_w1, FP_w1, FN_w1, precision_w1, recall_w1, f1_w1.
"""
import csv
from pathlib import Path
import sys

root = Path('.').resolve()
tmp = root / 'tmp'
if not tmp.exists():
    print('No tmp/ directory found; nothing to do')
    sys.exit(0)

for run_dir in sorted(tmp.iterdir()):
    if not run_dir.is_dir():
        continue
    # find w0 csv in this dir
    csv_w0_list = list(run_dir.glob('*_variant_summary_w0.csv'))
    if not csv_w0_list:
        continue
    csv_w0 = csv_w0_list[0]
    # corresponding w1 dir: replace trailing '_w0' with '_w1' in dir name if present
    run_name = run_dir.name
    if run_name.endswith('_w0'):
        sibling_name = run_name[:-2] + '1'
    else:
        sibling_name = run_name + '_w1'
    sibling_dir = tmp / sibling_name
    csv_w1 = None
    if sibling_dir.exists() and sibling_dir.is_dir():
        w1_list = list(sibling_dir.glob('*_variant_summary_w1.csv'))
        if w1_list:
            csv_w1 = w1_list[0]
    # fallback: try to find any w1 CSV anywhere under tmp with same pair prefix
    if csv_w1 is None:
        pair_prefix = csv_w0.name.replace('_variant_summary_w0.csv', '')
        # search for any file matching prefix + '_variant_summary_w1.csv'
        candidates = list(tmp.rglob(pair_prefix + '_variant_summary_w1.csv'))
        if candidates:
            csv_w1 = candidates[0]
    if csv_w1 is None:
        print('No matching w1 CSV found for', csv_w0)
        continue
    print('Processing', csv_w0.name, 'and', csv_w1.name)
    # load w1 metrics into dict by variant
    w1_map = {}
    with csv_w1.open('r', newline='') as f1:
        r1 = csv.DictReader(f1)
        for row in r1:
            var = row.get('variant')
            if not var:
                continue
            w1_map[var] = {
                'TP_w1': row.get('TP', ''), 'FP_w1': row.get('FP', ''), 'FN_w1': row.get('FN', ''),
                'precision_w1': row.get('precision', ''), 'recall_w1': row.get('recall', ''), 'f1_w1': row.get('f1', '')
            }

    # read w0 and write augmented temp file
    out_path = run_dir / (csv_w0.stem + '_aug.csv')
    with csv_w0.open('r', newline='') as f0, out_path.open('w', newline='') as fo:
        r0 = csv.DictReader(f0)
        fieldnames = list(r0.fieldnames or [])
        # add new columns if not present
        extras = ['TP_w1', 'FP_w1', 'FN_w1', 'precision_w1', 'recall_w1', 'f1_w1']
        for ex in extras:
            if ex not in fieldnames:
                fieldnames.append(ex)
        w = csv.DictWriter(fo, fieldnames=fieldnames)
        w.writeheader()
        for row in r0:
            var = row.get('variant')
            if var and var in w1_map:
                row.update(w1_map[var])
            else:
                for ex in extras:
                    row.setdefault(ex, '')
            w.writerow(row)
    # replace original with augmented file (keep backup)
    bak = run_dir / (csv_w0.name + '.bak')
    csv_w0.rename(bak)
    out_path.rename(csv_w0)
    print('Wrote augmented CSV', csv_w0)

print('Done')
