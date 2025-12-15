# Foldseek 3Di Alphabet Benchmark Results

## Overview
Comprehensive performance comparison between **Foldseek deployed 10f alphabet** and **custom-learned 8f alphabet** on SCOP validation set.

## Benchmark Configuration

| Parameter | Value |
|-----------|-------|
| **Query Set** | SCOP Validation (pdbs_val.txt) |
| **Total Domains** | 4,705 |
| **SCOP-labeled Queries** | 1,433 (30.5%) |
| **Alignment Method** | Smith-Waterman (SSW) |
| **Reference Baseline** | TMalign structure alignment |

---

## Raw Performance Metrics

### Foldseek 10f (Official Deployed Alphabet)

| Metric | AUC Value | TMalign % | vs 8f Improvement |
|--------|-----------|-----------|-------------------|
| **Family** | 0.8260 | 89.99% | +8.81% |
| **Superfamily** | 0.4536 | 68.52% | +18.14% |
| **Fold** | 0.1457 | 52.88% | +33.16% |

### Custom 8f (VQ-VAE Learned Alphabet)

| Metric | AUC Value | TMalign % |
|--------|-----------|-----------|
| **Family** | 0.7534 | 81.18% |
| **Superfamily** | 0.3335 | 50.38% |
| **Fold** | 0.0543 | 19.72% |

---

## Normalized Scores (TMalign Reference)

**Normalization Formula:**
```
Normalized Score = (Family/0.928162 + Superfamily/0.662063 + Fold/0.275436) / 3
```

| Alphabet | Normalized Score | Improvement vs 8f |
|----------|------------------|-------------------|
| **10f (Official)** | **0.7013** | **+39.08%** ⬆️ |
| **8f (Custom)** | 0.5042 | Baseline |

---

## Key Findings

### 🎯 Overall Performance
- **Foldseek 10f outperforms custom 8f by 39.1%** (normalized score)
- Official deployed 10f alphabet is substantially superior across all SCOP levels

### 📊 Category-wise Analysis

1. **Fold Classification** (+33.16% improvement)
   - 10f: 52.88% of TMalign performance
   - 8f: 19.72% of TMalign performance
   - **Largest gap**: Fold prediction is particularly weak in 8f; 10f demonstrates significant advantage

2. **Superfamily Classification** (+18.14% improvement)
   - 10f: 68.52% of TMalign performance
   - 8f: 50.38% of TMalign performance
   - Consistent improvement across mid-level hierarchies

3. **Family Classification** (+8.81% improvement)
   - 10f: 89.99% of TMalign performance
   - 8f: 81.18% of TMalign performance
   - Both perform reasonably well; gap narrows at finest granularity

### 💡 Interpretation
- **10f alphabet design** is optimized for structural hierarchy (Fold > Superfamily > Family)
- **8f model** captures family-level details well but struggles with coarser fold distinctions
- The 39% overall improvement suggests 10f better leverages 3Di information across scales

---

## Technical Details

### Data Processing
- **Extraction**: Foldseek `createdb` → `*_ss` files → 3Di tokens
- **Query Target**: Consolidated `tmp/target.fasta` (4,705 sequences)
- **Alignment**: SSW with `s_10f.mat` substitution matrix (gap open=8, extend=2)
- **ROC Calculation**: Per-query Family/Superfamily/Fold detection rates

### Evaluation Criteria
- **Family**: Same SCOP Family classification
- **Superfamily**: Same Superfamily classification  
- **Fold**: Same Fold classification
- AUC = average detection rate across queries

---

## Recommendations

1. **For Production Use**: Deploy **Foldseek 10f** (39% better performance)
2. **For Research**: 
   - Use 10f as strong baseline for new 3Di alphabets
   - 8f remains valuable for interpretability studies (simpler, trained alphabet)
3. **Future Optimization**: Explore 10f + fine-tuning on specific fold families

---

## Files Referenced

- **Benchmark Script**: `foldseek_10f/run-benchmark.sh`
- **3Di Extraction**: `foldseek_10f/extract_3di_from_createdb.py`
- **Result Data**: 
  - `foldseek_10f/tmp/result.rocx` (1,433 query results)
  - `foldseek_10f/tmp/target.fasta` (4,705 sequences)
  - `foldseek_10f/tmp/alignments/*.m8` (SSW outputs)

---

## Baseline References

**TMalign Performance (Structure Alignment)**
- Family: 0.928162
- Superfamily: 0.662063
- Fold: 0.275436

Source: `training_3di_gpu_8f/learnAlphabet.sh` (lines 81-84)

---

**Generated**: December 12, 2025  
**Benchmark Type**: 3Di Sequence Alignment on SCOP v1.75
