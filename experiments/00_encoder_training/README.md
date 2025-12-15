# Experiment 01: 3Di Encoder Training

## Overview
Training and evaluation of three custom 3Di sequence encoders designed for fold structure representation learning.

## Objectives
- Train 3Di encoders with different feature dimensions (8f, 9f, 10f)
- Evaluate encoder quality through CP detection and structural alignment tasks
- Compare custom encoders with Foldseek's official 10f encoder

## Encoder Variants

### 8f (8-feature custom encoder)
- **Dimensions**: 8 features per position
- **Status**: OPTIMAL for CP detection
- **Performance**: Highest SSW scores across all benchmarks (AS: 332 on Sap_C CP)
- **Use Case**: Recommended for circular permutation detection
- **Location**: `encoders_and_tools/training_3di_gpu_8f/`

### 9f (9-feature custom encoder)
- **Dimensions**: 9 features per position
- **Status**: Per-domain optimized
- **Performance**: Good on trained domains (173 on Sap_C), poor generalization (46 on 2hga vs 2vsv)
- **Use Case**: Specialized for known domain families
- **Location**: `encoders_and_tools/training_3di_gpu_9f/`
- **Note**: Exhibits overfitting characteristics

### 10f (10-feature encoders)
- **Custom 10f**:
  - **Dimensions**: 10 features per position
  - **Performance**: Moderate, stable (136-127 AS range)
  - **Location**: `encoders_and_tools/training_3di_gpu_10f/`
  
- **Foldseek 10f**:
  - **Dimensions**: 10 features per position
  - **Status**: Official reference implementation
  - **Source**: ESMFold/Foldseek official
  - **Use Case**: General-purpose structural representation

## Model Structure
All encoders consist of:
- **encoder.pt**: Trained encoder weights
- **decoder.pt**: Decoder weights for downstream tasks
- **states.txt**: Training state checkpoint

## Performance Comparison

### SSW Alignment Scores (Higher is Better)
```
Sap_C vs Sap_C_CP:
  8f:  AS=332 (ratio: 115.5, very clear)
  9f:  AS=173 (ratio: ~1.9, ambiguous)
  10f: AS=136 (ratio: 17.9, ambiguous)

2hga_trim vs 2vsv_PDZ:
  8f:  AS=142
  9f:  AS=46 (significant drop, overfitting)
  10f: AS=80

2vsv_PDZ vs 2z9i_pdz:
  8f:  AS=168
  9f:  AS=114
  10f: AS=127
```

### Key Findings
1. **8f superiority**: 2.4x better than 9f and 10f on Sap_C CP (332 vs 173 vs 136)
2. **9f overfitting**: Trained well on Sap_C but fails on other domains
3. **10f stability**: Consistent but suboptimal performance across domains
4. **Generalization**: 8f > 10f (custom) > 9f > 10f (Foldseek)

## Training Details

### Data
- **Positive CP pairs**: 451 pairs (cp_positive_pairs.tsv)
- **Non-CP homologs**: 1,593 pairs (noncp_homolog_pairs.tsv)
- **Domain structures**: 13,712 domains from SCOPE 2.08 (scope_pdb/)

### Training Configuration
- **Feature dimensions**: 8, 9, 10
- **Architecture**: Multi-layer perceptron with residual connections
- **Optimization**: SGD with momentum
- **Validation**: Stratified 80/20 split

### Loss Functions
- Primary: Classification loss (CP vs non-CP)
- Secondary: Reconstruction loss (decoder)
- Regularization: L2 weight decay

## Validation Strategy

### Benchmarks
1. **Circular Permutation Detection** (pairwise_3di_pipeline.py)
   - Input: PDB structures
   - Output: SSW alignment scores with 3Di sequences
   - Evaluation: Classification accuracy (CP vs non-CP)

2. **Structural Alignment** (TM-align)
   - Compares 3Di-based alignment with structure-based metrics
   - TM-score validation

3. **Trimmed Sequence Comparison**
   - Tests robustness to sequence variations
   - Handles domain boundary artifacts

## Usage

### Loading a Trained Encoder
```python
import torch
from pathlib import Path

# Detect available encoders
encoder_path = Path('encoders_and_tools/training_3di_gpu_8f')
encoder = torch.load(encoder_path / 'encoder.pt')
decoder = torch.load(encoder_path / 'decoder.pt')
```

### Converting PDB to 3Di
```bash
python pdb_to_3di.py --pdb_file 2hga.pdb --encoder 8f --output 2hga.3di
```

### Running SSW Alignment
```bash
# Requires s_8f.mat (8f substitution matrix)
ssw_binary -q query.3di -d subject.3di -a s_8f.mat
```

## Directory Structure
```
encoders_and_tools/
├── training_3di_gpu_8f/
│   ├── encoder.pt       # 8-feature encoder weights
│   ├── decoder.pt       # Decoder for 8f
│   └── states.txt       # Training checkpoint
├── training_3di_gpu_9f/
│   ├── encoder.pt       # 9-feature encoder weights
│   ├── decoder.pt
│   └── states.txt
├── training_3di_gpu_10f/
│   ├── encoder.pt       # 10-feature encoder weights
│   ├── decoder.pt
│   └── states.txt
├── ssw/                 # Smith-Waterman binary & libraries
├── s_8f.mat            # 8f 3Di substitution matrix
├── s_10f.mat           # 10f 3Di substitution matrix
└── sub_score.mat       # General substitution matrix
```

## Next Steps
1. ✅ Train 3 encoder variants (8f, 9f, 10f)
2. ✅ Benchmark on CP detection (Exp 02, 04)
3. ✅ Compare with official Foldseek 10f
4. ⚠️ Fine-tune 9f to improve generalization
5. ⚠️ Test on additional domain families

## References
- Exp 02: 8f/10f SSW comparison results
- Exp 03: X1/X2 CP detection with 8f
- Exp 04: PDZ domain sanity check validation
