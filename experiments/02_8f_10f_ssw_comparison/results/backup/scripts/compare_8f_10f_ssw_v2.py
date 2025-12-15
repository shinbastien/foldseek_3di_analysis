#!/usr/bin/env python3
"""
Compare SSW alignment scores between 8f (learned) and 10f (official) alphabets
on identical query-target pairs from SCOP validation set.

Usage:
    python compare_8f_10f.py [--plot-range {full|95pct}]
    
    --plot-range: Control plot axis ranges
        full   - Show all data points (may be dominated by outliers)
        95pct  - Show 95% percentile range for better visibility (default)
"""

import os
import sys
import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Paths
PATH_8F_ALIGN = "/mnt/scratch/jugipalace/foldseek_new_3di/training_3di_gpu_8f/tmp/alignments"
PATH_10F_ALIGN = "/mnt/scratch/jugipalace/foldseek_new_3di/foldseek_10f/tmp/alignments"
OUTPUT_DIR = "/mnt/scratch/jugipalace/foldseek_new_3di/foldseek_10f"

def load_alignments(align_dir):
    """Load all .m8 files and return dict of {(query, target): score}"""
    pairs = {}
    m8_files = glob.glob(os.path.join(align_dir, "*.m8"))
    print(f"Loading alignments from {align_dir}...")
    print(f"Found {len(m8_files)} .m8 files")
    
    for m8_file in m8_files:
        with open(m8_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    query, target, score = parts[0], parts[1], int(parts[2])
                    pairs[(query, target)] = score
    
    print(f"Loaded {len(pairs)} alignment pairs")
    return pairs

def compare_alignments(pairs_8f, pairs_10f):
    """Find common pairs and compute statistics"""
    keys_8f = set(pairs_8f.keys())
    keys_10f = set(pairs_10f.keys())
    
    common_keys = keys_8f & keys_10f
    print(f"\nCommon pairs: {len(common_keys)}")
    print(f"8f-only: {len(keys_8f - keys_10f)}")
    print(f"10f-only: {len(keys_10f - keys_8f)}")
    
    # Extract scores for common pairs
    scores_8f = [pairs_8f[k] for k in common_keys]
    scores_10f = [pairs_10f[k] for k in common_keys]
    
    df = pd.DataFrame({
        '8f_score': scores_8f,
        '10f_score': scores_10f
    })
    
    df['score_diff'] = df['10f_score'] - df['8f_score']
    df['pct_change'] = (df['score_diff'] / df['8f_score'].replace(0, 1)) * 100
    
    return df, common_keys

def compute_statistics(df):
    """Compute detailed statistics"""
    stats_dict = {
        '8f_mean': df['8f_score'].mean(),
        '8f_std': df['8f_score'].std(),
        '8f_median': df['8f_score'].median(),
        '8f_min': df['8f_score'].min(),
        '8f_max': df['8f_score'].max(),
        '10f_mean': df['10f_score'].mean(),
        '10f_std': df['10f_score'].std(),
        '10f_median': df['10f_score'].median(),
        '10f_min': df['10f_score'].min(),
        '10f_max': df['10f_score'].max(),
        'mean_diff': df['score_diff'].mean(),
        'median_diff': df['score_diff'].median(),
        'pct_change_mean': df['pct_change'].mean(),
        'pct_change_median': df['pct_change'].median(),
    }
    
    # Correlation and paired t-test
    correlation = df['8f_score'].corr(df['10f_score'])
    t_stat, p_val = stats.ttest_rel(df['10f_score'], df['8f_score'])
    
    stats_dict['correlation'] = correlation
    stats_dict['t_stat'] = t_stat
    stats_dict['p_value'] = p_val
    
    return stats_dict

def create_plots(df, output_dir, plot_range='95pct'):
    """Create visualization plots
    
    Args:
        df: DataFrame with comparison data
        output_dir: Output directory
        plot_range: 'full' for all data, '95pct' for 95% percentile range
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Calculate ranges based on plot_range parameter
    if plot_range == '95pct':
        score_low = df[['8f_score', '10f_score']].quantile(0.025).min()
        score_high = df[['8f_score', '10f_score']].quantile(0.975).max()
        diff_low = df['score_diff'].quantile(0.025)
        diff_high = df['score_diff'].quantile(0.975)
        pct_low = df['pct_change'][np.isfinite(df['pct_change'])].quantile(0.025)
        pct_high = df['pct_change'][np.isfinite(df['pct_change'])].quantile(0.975)
        range_label = '95% range'
    else:  # 'full'
        score_low = df[['8f_score', '10f_score']].min().min()
        score_high = df[['8f_score', '10f_score']].max().max()
        diff_low = df['score_diff'].min()
        diff_high = df['score_diff'].max()
        pct_low = df['pct_change'][np.isfinite(df['pct_change'])].min()
        pct_high = df['pct_change'][np.isfinite(df['pct_change'])].max()
        range_label = 'full range'
    
    fig = plt.figure(figsize=(16, 12))
    
    # 1. Scatter plot: 8f vs 10f
    ax1 = plt.subplot(2, 3, 1)
    df_filtered = df[(df['8f_score'] >= score_low) & (df['8f_score'] <= score_high) &
                     (df['10f_score'] >= score_low) & (df['10f_score'] <= score_high)]
    ax1.scatter(df_filtered['8f_score'], df_filtered['10f_score'], alpha=0.3, s=1)
    ax1.plot([score_low, score_high], [score_low, score_high], 'r--', lw=2, label='y=x')
    ax1.set_xlabel('8f SSW Score', fontsize=11)
    ax1.set_ylabel('10f SSW Score', fontsize=11)
    ax1.set_title('8f vs 10f SSW Scores\n({}, n={:,})'.format(range_label, len(df_filtered)), fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Score difference distribution
    ax2 = plt.subplot(2, 3, 2)
    diff_filtered = df['score_diff'][(df['score_diff'] >= diff_low) & (df['score_diff'] <= diff_high)]
    ax2.hist(diff_filtered, bins=100, edgecolor='black', alpha=0.7)
    ax2.axvline(df['score_diff'].mean(), color='r', linestyle='--', lw=2, label=f'Mean: {df["score_diff"].mean():.2f}')
    ax2.axvline(df['score_diff'].median(), color='g', linestyle='--', lw=2, label=f'Median: {df["score_diff"].median():.2f}')
    ax2.set_xlabel('10f - 8f Score Difference', fontsize=11)
    ax2.set_ylabel('Frequency', fontsize=11)
    ax2.set_title('Score Difference Distribution\n({}, n={:,})'.format(range_label, len(diff_filtered)), fontsize=12, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    # 3. Percentage change distribution
    ax3 = plt.subplot(2, 3, 3)
    valid_pct = df['pct_change'][np.isfinite(df['pct_change'])]
    pct_filtered = valid_pct[(valid_pct >= pct_low) & (valid_pct <= pct_high)]
    ax3.hist(pct_filtered, bins=100, edgecolor='black', alpha=0.7, color='orange')
    ax3.axvline(valid_pct.mean(), color='r', linestyle='--', lw=2, label=f'Mean: {valid_pct.mean():.2f}%')
    ax3.set_xlabel('Percentage Change (%)', fontsize=11)
    ax3.set_ylabel('Frequency', fontsize=11)
    ax3.set_title('10f vs 8f Percentage Change\n({}, n={:,})'.format(range_label, len(pct_filtered)), fontsize=12, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3, axis='y')
    
    # 4. Distribution comparison (8f vs 10f)
    ax4 = plt.subplot(2, 3, 4)
    scores_8f_filtered = df['8f_score'][(df['8f_score'] >= score_low) & (df['8f_score'] <= score_high)]
    scores_10f_filtered = df['10f_score'][(df['10f_score'] >= score_low) & (df['10f_score'] <= score_high)]
    ax4.hist(scores_8f_filtered, bins=100, alpha=0.5, label='8f', edgecolor='black')
    ax4.hist(scores_10f_filtered, bins=100, alpha=0.5, label='10f', edgecolor='black')
    ax4.set_xlabel('SSW Score', fontsize=11)
    ax4.set_ylabel('Frequency', fontsize=11)
    ax4.set_title('SSW Score Distribution Comparison\n({})'.format(range_label), fontsize=12, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3, axis='y')
    
    # 5. Box plot comparison
    ax5 = plt.subplot(2, 3, 5)
    box_data = [df['8f_score'], df['10f_score']]
    bp = ax5.boxplot(box_data, labels=['8f', '10f'], patch_artist=True)
    for patch, color in zip(bp['boxes'], ['lightblue', 'lightcoral']):
        patch.set_facecolor(color)
    ax5.set_ylabel('SSW Score', fontsize=11)
    ax5.set_title('SSW Score Box Plot Comparison', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3, axis='y')
    
    # 6. Statistics text box
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')
    
    stats_text = f"""
    STATISTICS (Common Pairs: {len(df)})
    
    8f:
      Mean: {df['8f_score'].mean():.2f}
      Median: {df['8f_score'].median():.2f}
      Std Dev: {df['8f_score'].std():.2f}
      Min-Max: {df['8f_score'].min():.0f}-{df['8f_score'].max():.0f}
    
    10f:
      Mean: {df['10f_score'].mean():.2f}
      Median: {df['10f_score'].median():.2f}
      Std Dev: {df['10f_score'].std():.2f}
      Min-Max: {df['10f_score'].min():.0f}-{df['10f_score'].max():.0f}
    
    Difference (10f - 8f):
      Mean: {df['score_diff'].mean():.2f}
      Median: {df['score_diff'].median():.2f}
      % Change: {valid_pct.mean():.2f}%
    
    Correlation: {df['8f_score'].corr(df['10f_score']):.4f}
    Paired t-test p-value: {stats.ttest_rel(df['10f_score'], df['8f_score'])[1]:.2e}
    """
    
    ax6.text(0.1, 0.95, stats_text, transform=ax6.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    # Save with range suffix
    suffix = '_95pct' if plot_range == '95pct' else '_full'
    output_filename = f'ssw_comparison_8f_vs_10f{suffix}.png'
    plt.savefig(os.path.join(output_dir, output_filename), dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {os.path.join(output_dir, output_filename)}")
    plt.close()

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Compare SSW alignment scores between 8f and 10f alphabets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python compare_8f_10f.py                    # Default: 95% range
  python compare_8f_10f.py --plot-range 95pct # 95% percentile range
  python compare_8f_10f.py --plot-range full  # Full data range
        """
    )
    parser.add_argument('--plot-range', choices=['full', '95pct'], default='95pct',
                        help='Plot axis range: full (all data) or 95pct (95%% percentile, default)')
    
    args = parser.parse_args()
    
    print(f"Plot range mode: {args.plot_range}")
    print("")
    
    # Load alignments
    pairs_8f = load_alignments(PATH_8F_ALIGN)
    pairs_10f = load_alignments(PATH_10F_ALIGN)
    
    # Compare
    df, common_keys = compare_alignments(pairs_8f, pairs_10f)
    
    # Statistics
    stats_dict = compute_statistics(df)
    
    # Print summary
    print("\n" + "="*70)
    print("8f vs 10f SSW SCORE COMPARISON")
    print("="*70)
    print(f"\n8f Statistics:")
    print(f"  Mean:   {stats_dict['8f_mean']:.2f}")
    print(f"  Median: {stats_dict['8f_median']:.2f}")
    print(f"  Std:    {stats_dict['8f_std']:.2f}")
    print(f"  Range:  {stats_dict['8f_min']:.0f} - {stats_dict['8f_max']:.0f}")
    
    print(f"\n10f Statistics:")
    print(f"  Mean:   {stats_dict['10f_mean']:.2f}")
    print(f"  Median: {stats_dict['10f_median']:.2f}")
    print(f"  Std:    {stats_dict['10f_std']:.2f}")
    print(f"  Range:  {stats_dict['10f_min']:.0f} - {stats_dict['10f_max']:.0f}")
    
    print(f"\nDifference (10f - 8f):")
    print(f"  Mean:           {stats_dict['mean_diff']:.2f}")
    print(f"  Median:         {stats_dict['median_diff']:.2f}")
    print(f"  % Change:       {stats_dict['pct_change_mean']:.2f}%")
    print(f"  Correlation:    {stats_dict['correlation']:.4f}")
    
    print(f"\nPaired t-test:")
    print(f"  t-statistic:    {stats_dict['t_stat']:.4f}")
    print(f"  p-value:        {stats_dict['p_value']:.2e}")
    if stats_dict['p_value'] < 0.001:
        print(f"  Result:         *** HIGHLY SIGNIFICANT (p < 0.001) ***")
    elif stats_dict['p_value'] < 0.05:
        print(f"  Result:         ** SIGNIFICANT (p < 0.05) **")
    else:
        print(f"  Result:         Not significant")
    
    print("\n" + "="*70)
    
    # Save CSV
    csv_path = os.path.join(OUTPUT_DIR, 'ssw_comparison_8f_vs_10f.csv')
    df.to_csv(csv_path, index=False)
    print(f"CSV saved to: {csv_path}")
    
    # Create plots with specified range
    create_plots(df, OUTPUT_DIR, plot_range=args.plot_range)
    
    print("\nComparison complete!")

if __name__ == '__main__':
    main()
