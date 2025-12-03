"""
plot_all_runs.py

Collect variant_summary CSVs under tmp/, merge them (adding a 'w' column if missing),
create a combined CSV, and draw the three plots using plot_3di_results.

Outputs go to scripts/plots/ by default.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
TMP = ROOT / 'tmp'
OUT_DIR = Path(__file__).resolve().parent / 'plots'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# find variant summary CSVs (match both old and new filenames, e.g. *_variant_summary.csv and *_variant_summary_w{w}.csv)
csvs = sorted(TMP.glob('*/*_variant_summary*.csv'))
if not csvs:
    print('No variant_summary CSVs found under tmp/. Run the pipeline first.')
    sys.exit(2)

print('Found CSVs:')
for c in csvs:
    print(' -', c)

try:
    import pandas as pd
except Exception as e:
    print('pandas is required to merge CSVs and plot. Install it in your environment.', file=sys.stderr)
    raise

# load and merge, add w if missing (assume w=0 for these runs)
dfs = []
for p in csvs:
    df = pd.read_csv(p)
    if 'w' not in df.columns:
        df['w'] = 0
    # add source file column
    df['source_csv'] = str(p)
    dfs.append(df)
combined = pd.concat(dfs, ignore_index=True, sort=False)
combined_path = Path(__file__).resolve().parent / 'combined_variant_summary.csv'
# compute/repair block_f1 if range columns exist
if 'range_target_f1' in combined.columns and 'range_query_f1' in combined.columns:
    combined['block_f1'] = 0.5 * (combined['range_target_f1'].fillna(0.0) + combined['range_query_f1'].fillna(0.0))
else:
    # fallback: preserve existing block_f1 or set zero
    if 'block_f1' not in combined.columns:
        combined['block_f1'] = 0.0

combined.to_csv(combined_path, index=False)
print('Wrote combined CSV to', combined_path)

# --- enforce global variant ordering: prefer the canonical helper if available
try:
    # variant_order.py lives next to this script (same directory)
    from variant_order import apply_variant_order
    # apply_variant_order returns a new dataframe with the variant column set to an
    # ordered Categorical using the canonical order (e.g. ['tm','8f','9f','10f']).
    combined = apply_variant_order(combined, variant_col='variant')
except Exception:
    # fallback: build numeric variants ascending (8f,9f,10f,...) and place 'tm' at front
    try:
        variants_all = list(dict.fromkeys(combined['variant'].astype(str).tolist()))
        import re
        numeric_variants = [v for v in variants_all if re.search(r"(\d+)", v)]
        # sort by numeric part ascending
        def _num(v):
            m = re.search(r"(\d+)", str(v))
            return int(m.group(1)) if m else float('inf')
        ordered_numeric = sorted(set(numeric_variants), key=_num)
        ordered = [f for f in ordered_numeric]
        # if 'tm' appears anywhere, put it at the leftmost position
        if any(str(v).lower() == 'tm' for v in variants_all):
            ordered.insert(0, 'tm')
        # also include any other non-numeric variants at the end preserving appearance order
        others = [v for v in variants_all if v not in ordered]
        ordered.extend(others)
        combined['variant'] = pd.Categorical(combined['variant'].astype(str), categories=ordered, ordered=True)
    except Exception:
        pass

# now draw the two new consolidated plots using plot_3di_results
try:
    from plot_3di_results import plot_pair_f1_all, plot_block_f1_w1_all
except Exception as e:
    print('Failed to import plotting utilities from scripts/plot_3di_results.py:', e, file=sys.stderr)
    raise

# 1) 모든 pair + w를 한 플롯에: pair-level F1
out_png = OUT_DIR / 'pair_f1_allpairs.png'
print('Drawing pair-level F1 for all pairs ->', out_png)
try:
    # aggregated (all w on same plot)
    plot_pair_f1_all(combined, save_path=str(out_png), show=False, per_w=False)
    # per-w separate files
    plot_pair_f1_all(combined, save_path=str(OUT_DIR / 'pair_f1_allpairs.png'), show=False, per_w=True)
except Exception as e:
    print('Failed to draw pair_f1_all:', e, file=sys.stderr)

# 2) w=1에서 block-level F1: variant vs pair
out_png2 = OUT_DIR / 'block_f1_w1_allpairs.png'
print('Drawing block-level F1 (w=1) for all pairs ->', out_png2)
try:
    # default w=1 aggregated (do NOT draw per-w files)
    plot_block_f1_w1_all(combined, w_value=1, save_path=str(out_png2), show=False, per_w=False)
except Exception as e:
    print('Failed to draw block_f1_w1_all:', e, file=sys.stderr)

print('Done. Plots saved to', OUT_DIR)

# --- additional metric plots requested: ssw_identity/identity, alignment length, ssw_score, lctr
try:
    from plot_3di_results import plot_pair_metrics_over_pairs
    metrics = ['ssw_score', '3di_identity', '3di_identity_alignment_len', 'lctr', 'ssw_pairs_count']
    print('Drawing additional per-pair metric barplots ->', OUT_DIR)
    plot_pair_metrics_over_pairs(combined, metrics=metrics, agg='mean', save_dir=str(OUT_DIR), show=False)
except Exception as e:
    print('Failed to draw additional metric plots:', e, file=sys.stderr)

# also produce grouped per-variant barplots (w=1) and line plots per pair
try:
    from plot_3di_results import plot_metric_by_variant_pairs, plot_metric_variant_lines
    per_variant_metrics = ['ssw_score', '3di_identity', '3di_identity_alignment_len', 'lctr']
    for metric in per_variant_metrics:
        out_bar = OUT_DIR / f"{metric}_by_variant_bar.png"
        try:
            # keep w filtering but do not write 'w1' in title or filename
            plot_metric_by_variant_pairs(combined, metric=metric, w_value=1, save_path=str(out_bar), show=False)
        except Exception as e:
            print(f'Could not draw barplot for {metric}:', e, file=sys.stderr)
        out_line = OUT_DIR / f"{metric}_by_variant_lines.png"
        try:
            plot_metric_variant_lines(combined, metric=metric, save_path=str(out_line), show=False)
        except Exception as e:
            print(f'Could not draw line plot for {metric}:', e, file=sys.stderr)
except Exception as e:
    print('Failed to draw per-variant metric plots:', e, file=sys.stderr)

# --- draw sequence length per pair (target vs query) if possible
try:
    from plot_3di_results import plot_sequence_lengths_by_pair
    out_len = OUT_DIR / 'sequence_lengths_by_pair.png'
    print('Drawing sequence lengths per pair ->', out_len)
    try:
        plot_sequence_lengths_by_pair(combined, save_path=str(out_len), show=False)
    except Exception as e:
        print('Could not draw sequence lengths plot:', e, file=sys.stderr)
except Exception as e:
    print('Sequence-length plotting utility not available:', e, file=sys.stderr)
