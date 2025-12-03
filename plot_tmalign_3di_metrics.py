#!/usr/bin/env python3
"""
plot_tmalign_3di_metrics.py

Plotting script for tmalign_3di_match_pipeline.py output CSVs.

This script reads variant summary CSVs and produces:
1. Pair-level precision–recall scatter plots (one per w value)
2. Block-level sensitivity–specificity scatter plots (target and query sides)

Uses only pandas, numpy, and matplotlib (no seaborn).

Usage:
    python plot_tmalign_3di_metrics.py --csv-dir /path/to/tmp/run_dir

or provide explicit CSV paths:
    python plot_tmalign_3di_metrics.py --csvs file1.csv file2.csv ...

The script will:
- Read all *_variant_summary_w*.csv files from the provided directory
- Generate pair-level P/R scatter plots (one per w value)
- Generate block-level sensitivity/specificity plots (per pair, target/query sides)
- Save all plots to the current directory or --output-dir
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import argparse
import itertools
import sys


def load_all_csvs(csv_paths):
    """Load and concatenate multiple variant summary CSVs."""
    dfs = []
    for p in csv_paths:
        try:
            df = pd.read_csv(p)
            df['source_csv'] = str(p)
            dfs.append(df)
            print(f'Loaded {len(df)} rows from {p}')
        except Exception as e:
            print(f'Warning: Failed to load {p}: {e}', file=sys.stderr)
    if not dfs:
        print('Error: No CSVs loaded successfully', file=sys.stderr)
        sys.exit(1)
    combined = pd.concat(dfs, ignore_index=True, sort=False)
    print(f'Combined total: {len(combined)} rows from {len(dfs)} CSV(s)')
    return combined


def plot_pair_level_precision_recall(df, output_dir='.'):
    """
    Plot pair-level precision–recall scatter for each w value.

    x-axis: recall (pair-level)
    y-axis: precision (pair-level)
    Points colored by pair, marker shape by variant.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract unique pairs and variants
    pairs = sorted(df['pair'].dropna().unique())
    variants = sorted(df['variant'].dropna().unique())

    # Color palette for pairs (cycle through matplotlib default colors)
    pair_to_color = {p: c for p, c in zip(pairs, itertools.cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']))}
    
    # Marker shapes for variants
    variant_to_marker = {'8f': 'o', '9f': '^', '10f': 's'}
    # fallback for any other variants
    for v in variants:
        if v not in variant_to_marker:
            variant_to_marker[v] = 'D'

    # Group by w and produce one plot per w value
    w_values = sorted(df['w'].dropna().unique())
    for w_val in w_values:
        sub = df[df['w'] == w_val].copy()
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(7, 6))

        # Plot each (pair, variant) group
        for (pair, variant), g in sub.groupby(['pair', 'variant']):
            color = pair_to_color.get(pair, 'C0')
            marker = variant_to_marker.get(variant, 'o')
            ax.scatter(g['recall'], g['precision'],
                       facecolors='none', edgecolors=color, marker=marker, s=80, linewidths=1.5, alpha=0.9)

        ax.set_xlabel("Recall (pair-level)", fontsize=12)
        ax.set_ylabel("Precision (pair-level)", fontsize=12)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, linestyle="--", alpha=0.4)

        # Build two legends: one for pair (color), one for variant (marker)
        pair_handles = [Line2D([0], [0], marker='o', color='none',
                               markeredgecolor=pair_to_color[p], markersize=8, label=p, linestyle='')
                        for p in pairs]
        var_handles = [Line2D([0], [0], marker=variant_to_marker[v], color='k', markersize=8, label=v, linestyle='')
                       for v in variants]

        # Place legends side by side at bottom
        first_legend = ax.legend(handles=pair_handles, title="Pair (color)", loc="lower left", fontsize=8, framealpha=0.9)
        ax.add_artist(first_legend)
        ax.legend(handles=var_handles, title="Variant (marker)", loc="lower right", fontsize=8, framealpha=0.9)

        ax.set_title(f"Pair-level Precision–Recall (w={w_val})", fontsize=14, weight='bold')
        fig.tight_layout()
        
        out_path = output_dir / f"pair_level_precision_recall_w{w_val}.png"
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved pair-level P/R plot: {out_path}')


def plot_block_level_sens_spec(df, output_dir='.'):
    """
    Plot block-level sensitivity vs specificity scatter plots.

    For each pair, produce two plots:
    - Target side: spec_T_block vs sens_T_block
    - Query side: spec_Q_block vs sens_Q_block

    Points are colored by pair and shaped by variant.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Use a single w value (block metrics should not depend on w conceptually,
    # but we may have multiple w runs; pick the first for simplicity)
    if 'w' in df.columns and not df['w'].isna().all():
        w_ref = sorted(df['w'].dropna().unique())[0]
        df_filtered = df[df['w'] == w_ref].copy()
        print(f'Using w={w_ref} for block-level plots')
    else:
        df_filtered = df.copy()
        w_ref = 'NA'

    pairs = sorted(df_filtered['pair'].dropna().unique())
    variants = sorted(df_filtered['variant'].dropna().unique())
    variant_to_marker = {'8f': 'o', '9f': '^', '10f': 's'}
    for v in variants:
        if v not in variant_to_marker:
            variant_to_marker[v] = 'D'
    pair_to_color = {p: c for p, c in zip(pairs, itertools.cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']))}

    for pair in pairs:
        sub = df_filtered[df_filtered['pair'] == pair].copy()
        if sub.empty:
            continue

        # --- Target side plot ---
        fig, ax = plt.subplots(figsize=(6, 6))
        for variant, g in sub.groupby('variant'):
            color = pair_to_color[pair]
            marker = variant_to_marker.get(variant, 'o')
            # spec_T_block on x-axis, sens_T_block on y-axis
            ax.scatter(g['spec_T_block'], g['sens_T_block'],
                       marker=marker, edgecolors=color, facecolors='none', s=100, linewidths=2, label=variant)
        ax.set_xlabel("Specificity (target-side, residue-level)", fontsize=12)
        ax.set_ylabel("Sensitivity (target-side, residue-level)", fontsize=12)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(title="Variant", fontsize=9, loc='lower right', framealpha=0.9)
        ax.set_title(f"Target Sensitivity vs Specificity — {pair}\n(w={w_ref})", fontsize=12, weight='bold')
        fig.tight_layout()
        
        out_path = output_dir / f"block_target_sens_spec_{pair}_w{w_ref}.png"
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved target block plot: {out_path}')

        # --- Query side plot ---
        fig, ax = plt.subplots(figsize=(6, 6))
        for variant, g in sub.groupby('variant'):
            color = pair_to_color[pair]
            marker = variant_to_marker.get(variant, 'o')
            # spec_Q_block on x-axis, sens_Q_block on y-axis
            ax.scatter(g['spec_Q_block'], g['sens_Q_block'],
                       marker=marker, edgecolors=color, facecolors='none', s=100, linewidths=2, label=variant)
        ax.set_xlabel("Specificity (query-side, residue-level)", fontsize=12)
        ax.set_ylabel("Sensitivity (query-side, residue-level)", fontsize=12)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(title="Variant", fontsize=9, loc='lower right', framealpha=0.9)
        ax.set_title(f"Query Sensitivity vs Specificity — {pair}\n(w={w_ref})", fontsize=12, weight='bold')
        fig.tight_layout()
        
        out_path = output_dir / f"block_query_sens_spec_{pair}_w{w_ref}.png"
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved query block plot: {out_path}')

    # --- Combined target side plot (all pairs in one) ---
    fig, ax = plt.subplots(figsize=(8, 8))
    for pair in pairs:
        sub = df_filtered[df_filtered['pair'] == pair].copy()
        if sub.empty:
            continue
        for variant, g in sub.groupby('variant'):
            color = pair_to_color[pair]
            marker = variant_to_marker.get(variant, 'o')
            ax.scatter(g['spec_T_block'], g['sens_T_block'],
                       marker=marker, edgecolors=color, facecolors='none', s=100, linewidths=2)
    
    # Create two legends: one for pairs (colors), one for variants (markers)
    pair_handles = [plt.Line2D([0], [0], marker='o', color='w', markeredgecolor=pair_to_color[p],
                               markerfacecolor='none', markersize=8, linewidth=2, label=p) for p in pairs]
    variant_handles = [plt.Line2D([0], [0], marker=variant_to_marker[v], color='w', markeredgecolor='k',
                                  markerfacecolor='none', markersize=8, linewidth=2, label=v) for v in sorted(variant_to_marker.keys())]
    
    legend1 = ax.legend(handles=pair_handles, title="Pair", fontsize=7, loc='lower left', framealpha=0.9)
    ax.add_artist(legend1)
    ax.legend(handles=variant_handles, title="Variant", fontsize=8, loc='lower right', framealpha=0.9)
    
    ax.set_xlabel("Specificity (target-side, residue-level)", fontsize=12)
    ax.set_ylabel("Sensitivity (target-side, residue-level)", fontsize=12)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.set_title(f"Target Sensitivity vs Specificity — All Pairs\n(w={w_ref})", fontsize=12, weight='bold')
    fig.tight_layout()
    
    out_path = output_dir / f"block_target_sens_spec_all_pairs_w{w_ref}.png"
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved combined target block plot: {out_path}')

    # --- Combined query side plot (all pairs in one) ---
    fig, ax = plt.subplots(figsize=(8, 8))
    for pair in pairs:
        sub = df_filtered[df_filtered['pair'] == pair].copy()
        if sub.empty:
            continue
        for variant, g in sub.groupby('variant'):
            color = pair_to_color[pair]
            marker = variant_to_marker.get(variant, 'o')
            ax.scatter(g['spec_Q_block'], g['sens_Q_block'],
                       marker=marker, edgecolors=color, facecolors='none', s=100, linewidths=2)
    
    # Create two legends
    pair_handles = [plt.Line2D([0], [0], marker='o', color='w', markeredgecolor=pair_to_color[p],
                               markerfacecolor='none', markersize=8, linewidth=2, label=p) for p in pairs]
    variant_handles = [plt.Line2D([0], [0], marker=variant_to_marker[v], color='w', markeredgecolor='k',
                                  markerfacecolor='none', markersize=8, linewidth=2, label=v) for v in sorted(variant_to_marker.keys())]
    
    legend1 = ax.legend(handles=pair_handles, title="Pair", fontsize=7, loc='lower left', framealpha=0.9)
    ax.add_artist(legend1)
    ax.legend(handles=variant_handles, title="Variant", fontsize=8, loc='lower right', framealpha=0.9)
    
    ax.set_xlabel("Specificity (query-side, residue-level)", fontsize=12)
    ax.set_ylabel("Sensitivity (query-side, residue-level)", fontsize=12)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.set_title(f"Query Sensitivity vs Specificity — All Pairs\n(w={w_ref})", fontsize=12, weight='bold')
    fig.tight_layout()
    
    out_path = output_dir / f"block_query_sens_spec_all_pairs_w{w_ref}.png"
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved combined query block plot: {out_path}')


def main():
    parser = argparse.ArgumentParser(description='Plot tmalign_3di_match_pipeline variant summary metrics')
    parser.add_argument('--csv-dir', help='Directory containing *_variant_summary_w*.csv files')
    parser.add_argument('--csvs', nargs='+', help='Explicit list of CSV files to load')
    parser.add_argument('--output-dir', default='.', help='Directory to save output plots (default: current directory)')
    args = parser.parse_args()

    csv_paths = []
    if args.csvs:
        csv_paths = [Path(c) for c in args.csvs]
    elif args.csv_dir:
        csv_dir = Path(args.csv_dir)
        if not csv_dir.exists():
            print(f'Error: CSV directory not found: {csv_dir}', file=sys.stderr)
            sys.exit(1)
        csv_paths = sorted(csv_dir.glob('*_variant_summary_w*.csv'))
        if not csv_paths:
            print(f'Warning: No *_variant_summary_w*.csv files found in {csv_dir}', file=sys.stderr)
            # fallback: try any CSV in that directory
            csv_paths = sorted(csv_dir.glob('*.csv'))
    else:
        print('Error: Provide either --csv-dir or --csvs', file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if not csv_paths:
        print('Error: No CSV files found to process', file=sys.stderr)
        sys.exit(1)

    print(f'Found {len(csv_paths)} CSV file(s) to load')
    df = load_all_csvs(csv_paths)

    # Check for required columns
    required_cols = ['pair', 'variant', 'w', 'precision', 'recall',
                     'sens_T_block', 'spec_T_block', 'sens_Q_block', 'spec_Q_block']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f'Error: Missing required columns in CSV: {missing}', file=sys.stderr)
        print(f'Available columns: {list(df.columns)}', file=sys.stderr)
        sys.exit(1)

    # Generate plots
    print('\n--- Generating pair-level precision–recall plots ---')
    plot_pair_level_precision_recall(df, output_dir=args.output_dir)

    print('\n--- Generating block-level sensitivity–specificity plots ---')
    plot_block_level_sens_spec(df, output_dir=args.output_dir)

    print(f'\nAll plots saved to: {Path(args.output_dir).resolve()}')


if __name__ == "__main__":
    main()
