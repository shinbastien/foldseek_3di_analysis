# PDB to 3Di Converter

A script to convert PDB files to 3Di sequences.

## Key Features

1. ✅ **PDB File Input**: Accepts a single PDB file as input
2. ✅ **HETATM Filtering**: Automatically filters out unnecessary HETATM records
3. ✅ **3Di Sequence Generation**: Saves output in `{protein_name}_3di.txt` format
4. ✅ **Detailed Logging**: Progress logs for each processing step

## Usage

### Basic Usage

```bash
# Activate fold3di environment (required)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate fold3di

# Run the script
python pdb_to_3di.py <pdb_file> <model_directory>
```

### Example 1: Using training_3di_gpu_new model

```bash
python pdb_to_3di.py permutation_examples/2hga_pdz.pdb training_3di_gpu_new/tmp
```

Output: `permutation_examples/2hga_pdz_3di.txt`

### Example 2: Using training_3di_gpu_original model

```bash
python pdb_to_3di.py permutation_examples/Sap_C.pdb training_3di_gpu_original/tmp
```

Output: `permutation_examples/Sap_C_3di.txt`

### Example 3: Specifying custom output file

```bash
python pdb_to_3di.py permutation_examples/2z9i_pdz.pdb training_3di_gpu_new/tmp -o my_output.txt
```

Output: `my_output.txt`

## Requirements

### Model Directory Structure

The model directory **must** contain the following files:
- `encoder.pt` - Trained encoder model
- `states.txt` - Cluster centroids

### Error Handling

If files are missing, a clear error message will be displayed:

```
❌ Error: Required model files are missing:
  - encoder.pt not found in training_3di_gpu_new/tmp
  - states.txt not found in training_3di_gpu_new/tmp
```

## Output Format

Format of the generated `_3di.txt` file:

```
>protein_name
ABCdeFghIjKLMnoPQrSTuvWYZ...
```

- First line: FASTA format header (`>protein_name`)
- Second line: 3Di sequence (50 alphabet characters, X marks invalid positions)

## Sample Log Output

```
============================================================
PDB to 3Di Converter
============================================================
✓ Found PDB file: permutation_examples/2hga_pdz.pdb
✓ Found encoder model: training_3di_gpu_new/tmp/encoder.pt
✓ Found centroids: training_3di_gpu_new/tmp/states.txt

[Loading Models]
  - Loading encoder...
  - Loading centroids...
  - Centroids shape: (50, 16)

[Step 1] Reading PDB structure...
  - Input file: permutation_examples/2hga_pdz.pdb
  - Total residues: 95
  - Valid residues (non-HETATM): 95

[Step 2] Adjusting virtual C-beta positions...

[Step 3] Finding nearest neighbor residues...

[Step 4] Calculating angular features...
  - Residues with complete features: 93

[Step 5] Encoding features with neural network...
  - Encoded feature shape: (93, 16)

[Step 6] Discretizing to 3Di states...
  - Number of unique states used: 38

[Step 7] Creating final 3Di sequence...

[Step 8] Saving results...
  - Output saved to: permutation_examples/2hga_pdz_3di.txt

============================================================
✓ Conversion completed successfully!
============================================================
  Protein name: 2hga_pdz
  Total residues: 95
  Valid 3Di states: 93
  Invalid positions (X): 2
  3Di sequence length: 95
  Output file: permutation_examples/2hga_pdz_3di.txt
============================================================
```

## Batch Processing Example

To process multiple PDB files at once:

```bash
#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate fold3di

MODEL_DIR="training_3di_gpu_new/tmp"

for pdb_file in permutation_examples/*.pdb; do
    echo "Processing: $pdb_file"
    python pdb_to_3di.py "$pdb_file" "$MODEL_DIR"
    echo "---"
done
```
