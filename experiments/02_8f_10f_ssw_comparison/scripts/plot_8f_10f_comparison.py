import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
import argparse

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (16, 12)

def load_results(csv_file):
    """Load results from CSV file."""
    df = pd.read_csv(csv_file)
    return df

def create_plots(df, output_dir):
    """Create multiple plots for 8f and 10f comparison."""
    
    # Separate data by encoding
    df_8f = df[df['encoding'] == '8f'].copy()
    df_10f = df[df['encoding'] == '10f'].copy()
    
    # Create figure with subplots
    fig = plt.figure(figsize=(18, 14))
    
    # 1. Distribution of X1 scores for 8f and 10f
    ax1 = plt.subplot(3, 3, 1)
    ax1.hist(df_8f['score_x1'], bins=20, alpha=0.6, label='8f', color='blue', edgecolor='black')
    ax1.hist(df_10f['score_x1'], bins=20, alpha=0.6, label='10f', color='red', edgecolor='black')
    ax1.set_xlabel('X1 Score')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Distribution of X1 Scores (8f vs 10f)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Distribution of X2 scores for 8f and 10f
    ax2 = plt.subplot(3, 3, 2)
    ax2.hist(df_8f['score_x2'], bins=20, alpha=0.6, label='8f', color='blue', edgecolor='black')
    ax2.hist(df_10f['score_x2'], bins=20, alpha=0.6, label='10f', color='red', edgecolor='black')
    ax2.set_xlabel('X2 Score')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Distribution of X2 Scores (8f vs 10f)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Distribution of Δssw (score_diff) for 8f and 10f
    ax3 = plt.subplot(3, 3, 3)
    ax3.hist(df_8f['score_diff'], bins=20, alpha=0.6, label='8f', color='blue', edgecolor='black')
    ax3.hist(df_10f['score_diff'], bins=20, alpha=0.6, label='10f', color='red', edgecolor='black')
    ax3.set_xlabel('Δssw (X2 - X1)')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Distribution of Δssw (8f vs 10f)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Box plot of X1 scores
    ax4 = plt.subplot(3, 3, 4)
    data_x1 = [df_8f['score_x1'], df_10f['score_x1']]
    bp4 = ax4.boxplot(data_x1, labels=['8f', '10f'], patch_artist=True)
    for patch, color in zip(bp4['boxes'], ['lightblue', 'lightcoral']):
        patch.set_facecolor(color)
    ax4.set_ylabel('X1 Score')
    ax4.set_title('Box Plot of X1 Scores')
    ax4.grid(True, alpha=0.3, axis='y')
    
    # 5. Box plot of X2 scores
    ax5 = plt.subplot(3, 3, 5)
    data_x2 = [df_8f['score_x2'], df_10f['score_x2']]
    bp5 = ax5.boxplot(data_x2, labels=['8f', '10f'], patch_artist=True)
    for patch, color in zip(bp5['boxes'], ['lightblue', 'lightcoral']):
        patch.set_facecolor(color)
    ax5.set_ylabel('X2 Score')
    ax5.set_title('Box Plot of X2 Scores')
    ax5.grid(True, alpha=0.3, axis='y')
    
    # 6. Box plot of Δssw
    ax6 = plt.subplot(3, 3, 6)
    data_diff = [df_8f['score_diff'], df_10f['score_diff']]
    bp6 = ax6.boxplot(data_diff, labels=['8f', '10f'], patch_artist=True)
    for patch, color in zip(bp6['boxes'], ['lightblue', 'lightcoral']):
        patch.set_facecolor(color)
    ax6.set_ylabel('Δssw (X2 - X1)')
    ax6.set_title('Box Plot of Δssw')
    ax6.grid(True, alpha=0.3, axis='y')
    
    # 7. Scatter plot: X1 vs X2 for 8f
    ax7 = plt.subplot(3, 3, 7)
    ax7.scatter(df_8f['score_x1'], df_8f['score_x2'], alpha=0.6, s=50, color='blue')
    # Add diagonal line (y=x)
    max_val = max(df_8f['score_x1'].max(), df_8f['score_x2'].max())
    ax7.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='y=x')
    ax7.set_xlabel('X1 Score')
    ax7.set_ylabel('X2 Score')
    ax7.set_title('X1 vs X2 Scores (8f)')
    ax7.legend()
    ax7.grid(True, alpha=0.3)
    
    # 8. Scatter plot: X1 vs X2 for 10f
    ax8 = plt.subplot(3, 3, 8)
    ax8.scatter(df_10f['score_x1'], df_10f['score_x2'], alpha=0.6, s=50, color='red')
    # Add diagonal line (y=x)
    max_val = max(df_10f['score_x1'].max(), df_10f['score_x2'].max())
    ax8.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='y=x')
    ax8.set_xlabel('X1 Score')
    ax8.set_ylabel('X2 Score')
    ax8.set_title('X1 vs X2 Scores (10f)')
    ax8.legend()
    ax8.grid(True, alpha=0.3)
    
    # 9. Violin plot of Δssw
    ax9 = plt.subplot(3, 3, 9)
    data_for_violin = pd.concat([
        df_8f[['score_diff']].assign(encoding='8f'),
        df_10f[['score_diff']].assign(encoding='10f')
    ])
    sns.violinplot(data=data_for_violin, x='encoding', y='score_diff', ax=ax9, palette=['lightblue', 'lightcoral'])
    ax9.set_ylabel('Δssw (X2 - X1)')
    ax9.set_title('Violin Plot of Δssw')
    ax9.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/8f_10f_comparison_plots.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/8f_10f_comparison_plots.png")
    
    # Generate summary statistics
    print("\n" + "="*60)
    print("Summary Statistics for 8f encoding:")
    print("="*60)
    print(f"X1 Score - Mean: {df_8f['score_x1'].mean():.2f}, Std: {df_8f['score_x1'].std():.2f}")
    print(f"X2 Score - Mean: {df_8f['score_x2'].mean():.2f}, Std: {df_8f['score_x2'].std():.2f}")
    print(f"Δssw     - Mean: {df_8f['score_diff'].mean():.2f}, Std: {df_8f['score_diff'].std():.2f}")
    print(f"Δssw     - Min: {df_8f['score_diff'].min():.2f}, Max: {df_8f['score_diff'].max():.2f}")
    
    print("\n" + "="*60)
    print("Summary Statistics for 10f encoding:")
    print("="*60)
    print(f"X1 Score - Mean: {df_10f['score_x1'].mean():.2f}, Std: {df_10f['score_x1'].std():.2f}")
    print(f"X2 Score - Mean: {df_10f['score_x2'].mean():.2f}, Std: {df_10f['score_x2'].std():.2f}")
    print(f"Δssw     - Mean: {df_10f['score_diff'].mean():.2f}, Std: {df_10f['score_diff'].std():.2f}")
    print(f"Δssw     - Min: {df_10f['score_diff'].min():.2f}, Max: {df_10f['score_diff'].max():.2f}")
    
    # Count non-zero Δssw values
    print("\n" + "="*60)
    print("Non-zero Δssw values:")
    print("="*60)
    non_zero_8f = (df_8f['score_diff'] != 0).sum()
    non_zero_10f = (df_10f['score_diff'] != 0).sum()
    print(f"8f: {non_zero_8f}/{len(df_8f)} pairs ({100*non_zero_8f/len(df_8f):.1f}%)")
    print(f"10f: {non_zero_10f}/{len(df_10f)} pairs ({100*non_zero_10f/len(df_10f):.1f}%)")

def main():
    parser = argparse.ArgumentParser(description="Plot 8f and 10f SSW score comparisons.")
    parser.add_argument("csv_file", type=str, help="Path to the comparison_results.csv file.")
    parser.add_argument("--output_dir", type=str, default="./", help="Output directory for plots.")
    
    args = parser.parse_args()
    
    df = load_results(args.csv_file)
    create_plots(df, args.output_dir)

if __name__ == "__main__":
    main()
