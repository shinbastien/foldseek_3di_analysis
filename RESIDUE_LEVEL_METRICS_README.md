# Residue-Level Confusion Metrics — Implementation Summary

## Overview

This document describes the changes made to `tmalign_3di_match_pipeline.py` to compute residue-level (block-level) confusion metrics and the new plotting script `plot_tmalign_3di_metrics.py`.

## Changes to `tmalign_3di_match_pipeline.py`

### New Residue-Level Confusion Metrics

The pipeline now computes **residue-level confusion matrices** for both TARGET and QUERY sequences:

#### Target (B) Side Metrics
- **TP_T**: Number of residues covered by both TM-align and SSW/3Di
- **FP_T**: Number of residues covered by SSW/3Di but not TM-align
- **FN_T**: Number of residues covered by TM-align but not SSW/3Di
- **TN_T**: Number of residues not covered by either
- **sens_T_block**: Sensitivity = TP_T / (TP_T + FN_T) if denominator > 0 else 1.0
- **spec_T_block**: Specificity = TN_T / (TN_T + FP_T) if denominator > 0 else 1.0

#### Query (A) Side Metrics
- **TP_Q**: Number of query residues covered by both TM-align and SSW/3Di
- **FP_Q**: Number of query residues covered by SSW/3Di but not TM-align
- **FN_Q**: Number of query residues covered by TM-align but not SSW/3Di
- **TN_Q**: Number of query residues not covered by either
- **sens_Q_block**: Sensitivity = TP_Q / (TP_Q + FN_Q) if denominator > 0 else 1.0
- **spec_Q_block**: Specificity = TN_Q / (TN_Q + FP_Q) if denominator > 0 else 1.0

### Implementation Details

The computation happens in the per-variant loop after coverage arrays (`tm_cov`, `ssw_cov`) are built:

```python
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
```

Similarly for query side using `tm_cov_Q` and `ssw_cov_Q` arrays.

### Updated CSV Output

The per-variant CSV now includes these additional columns (in order):

```
pair, query_name, target_name, variant, w,
tm_pairs_count, ssw_pairs_count, TP, FP, FN, precision, recall, f1, jaccard,
TP_T, FP_T, FN_T, TN_T, sens_T_block, spec_T_block,
TP_Q, FP_Q, FN_Q, TN_Q, sens_Q_block, spec_Q_block,
T3_range, Ttm_range, tm_segs, overlap_T, len_T3, len_Ttm, precision_T, recall_T, range_target_f1,
Q3_range, Qtm_range, overlap_Q, len_Q3, len_Qtm, precision_Q, recall_Q, range_query_f1,
lctr, block_recall_tau, block_f1_calc,
ssw_score, 3di_identity, 3di_identity_match_count, 3di_identity_alignment_len, cp_ok, cutpoint_s, len_target,
tm_pairs_list, ssw_pairs_file, coverage_by_target, coverage_by_query
```

## New Plotting Script: `plot_tmalign_3di_metrics.py`

### Purpose
A standalone plotting script that reads variant summary CSVs and generates:
1. **Pair-level precision–recall scatter plots** (one per w value)
2. **Block-level sensitivity–specificity scatter plots** (per pair, separate for target and query)

### Dependencies
- pandas
- numpy
- matplotlib
- No seaborn required

### Usage

#### Basic Usage
```bash
# Plot all CSVs from a directory (e.g., tmp_2)
python plot_tmalign_3di_metrics.py --csv-dir /path/to/tmp_2

# Plot specific CSV files
python plot_tmalign_3di_metrics.py --csvs file1.csv file2.csv file3.csv

# Specify output directory
python plot_tmalign_3di_metrics.py --csv-dir tmp_2 --output-dir my_plots
```

### Generated Plots

#### 1. Pair-Level Precision–Recall Scatter
- **Filename**: `pair_level_precision_recall_w{w}.png`
- **X-axis**: Recall (pair-level)
- **Y-axis**: Precision (pair-level)
- **Points**: Each point represents one (pair, variant) combination
- **Color**: Represents pair (protein comparison)
- **Marker shape**: Represents variant (8f=circle, 9f=triangle, 10f=square)
- **One plot per w value**

#### 2. Block-Level Sensitivity–Specificity (Target Side)
- **Filename**: `block_target_sens_spec_{pair}_w{w}.png`
- **X-axis**: Specificity (target-side, residue-level)
- **Y-axis**: Sensitivity (target-side, residue-level)
- **Points**: Each point represents one variant
- **Color**: Pair color
- **Marker shape**: Variant
- **One plot per pair**

#### 3. Block-Level Sensitivity–Specificity (Query Side)
- **Filename**: `block_query_sens_spec_{pair}_w{w}.png`
- **X-axis**: Specificity (query-side, residue-level)
- **Y-axis**: Sensitivity (query-side, residue-level)
- **Points**: Each point represents one variant
- **Color**: Pair color
- **Marker shape**: Variant
- **One plot per pair**

### Example Workflow

```bash
# Step 1: Run the pipeline for multiple w values (generates CSVs with new metrics)
python tmalign_3di_match_pipeline.py pdb1.pdb pdb2.pdb --w-pair 0 --run-dir tmp/run_w0
python tmalign_3di_match_pipeline.py pdb1.pdb pdb2.pdb --w-pair 1 --run-dir tmp/run_w1
python tmalign_3di_match_pipeline.py pdb1.pdb pdb2.pdb --w-pair 2 --run-dir tmp/run_w2

# Step 2: Generate plots from all runs
python plot_tmalign_3di_metrics.py --csvs \
    tmp/run_w0/*_variant_summary_w0.csv \
    tmp/run_w1/*_variant_summary_w1.csv \
    tmp/run_w2/*_variant_summary_w2.csv \
    --output-dir plots/
```

## Interpretation of Metrics

### Pair-Level Metrics (existing)
- **Precision**: Of the 3Di/SSW pairs predicted, what fraction match TM-align pairs (within window w)?
- **Recall**: Of the TM-align pairs, what fraction are recovered by 3Di/SSW (within window w)?

### Residue-Level Block Metrics (new)
- **Sensitivity (Target)**: Of all target residues covered by TM-align, what fraction are also covered by 3Di/SSW?
- **Specificity (Target)**: Of all target residues NOT covered by TM-align, what fraction are also NOT covered by 3Di/SSW?
- **Sensitivity (Query)**: Of all query residues covered by TM-align, what fraction are also covered by 3Di/SSW?
- **Specificity (Query)**: Of all query residues NOT covered by TM-align, what fraction are also NOT covered by 3Di/SSW?

### Key Differences
- **Pair-level** metrics measure alignment pair matching (with tolerance window w)
- **Residue-level** metrics measure coverage agreement at individual residue positions (binary coverage, no window tolerance)

## Testing the Changes

### Verify Pipeline Syntax
```bash
python3 -m py_compile tmalign_3di_match_pipeline.py
```

### Verify Plotting Script Syntax
```bash
python3 -m py_compile plot_tmalign_3di_metrics.py
```

### Run a Test Pipeline
```bash
# Run with existing PDB files
python tmalign_3di_match_pipeline.py \
    path/to/pdb1.pdb \
    path/to/pdb2.pdb \
    --w-pair 1 \
    --run-dir test_run

# Check that CSV has new columns
head -1 test_run/*_variant_summary_w1.csv
```

### Generate Test Plots
```bash
# After pipeline run, generate plots
python plot_tmalign_3di_metrics.py --csv-dir test_run --output-dir test_plots
```

## Files Modified/Created

### Modified
- `tmalign_3di_match_pipeline.py`
  - Added residue-level confusion computation for target (lines ~2080-2095)
  - Added residue-level confusion computation for query (lines ~2096-2120)
  - Updated `variant_rows` dict to include new metrics (lines ~2160-2175)
  - Updated CSV field list (lines ~2230-2240)

### Created
- `plot_tmalign_3di_metrics.py` (new standalone plotting script)
- `RESIDUE_LEVEL_METRICS_README.md` (this document)

## Notes

- The residue-level metrics are computed **after** circular permutation rotation, so they measure coverage agreement in the rotated coordinate frame
- Default values (if computation fails): TP=FP=FN=TN=0, sensitivity=specificity=1.0
- Block-level plots use only one w value (typically the first found) since residue-level coverage is conceptually independent of the pair-matching window w
- All plotting uses only pandas/numpy/matplotlib (no seaborn)

## Future Enhancements

Potential extensions:
1. Add ROC-like curves by varying a threshold parameter
2. Compute Matthews correlation coefficient (MCC) from the confusion matrices
3. Add per-variant aggregate statistics to plot titles
4. Generate combined multi-panel figures with all pairs and metrics

---

**Contact**: If you encounter issues or have questions, check the pipeline logs in `run_logs/` and ensure all required dependencies (pandas, numpy, matplotlib) are installed.
