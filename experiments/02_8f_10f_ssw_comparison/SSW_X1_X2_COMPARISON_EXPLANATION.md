# SSW X1 vs X2 Comparison: Detailed Explanation

## Overview
This document explains how we calculate and compare SSW (Smith-Waterman) alignment scores between X1 (single sequence) and X2 (doubled sequence) to detect circular permutations.

---

## 1. X1 Calculation (Single Sequence Alignment)

### Pseudo Code
```
INPUT:
  - query_seq: Original 3Di sequence (e.g., "VLLVVCCVPP...")
  - target_seq: Original 3Di sequence
  - matrix: Substitution matrix (8f or 10f)
  - gap_open: 8
  - gap_extend: 2

PROCESS:
  1. Create temporary FASTA files:
     - query.fasta: ">query_name\nVLLVVCCVPP..."
     - target.fasta: ">target_name\nVLLVVCCVPP..."
  
  2. Run SSW (Smith-Waterman) local alignment:
     ssw_test -o 8 -e 2 -a matrix.mat query.fasta target.fasta
  
  3. Parse output to extract optimal alignment score:
     optimal_alignment_score: 480
  
  4. Return score_x1 = 480

OUTPUT:
  - score_x1: Integer (SSW optimal alignment score)
  - alignment visualization with match/mismatch positions
```

### Real Python Code (from compare_8f_10f_ssw.py)
```python
def compute_ssw_score_with_binary(ssw_binary: str, query_seq: str, target_seq: str,
                                   query_name: str, target_name: str, encoding: str,
                                   is_x2: bool, output_dir: str, gap_open: int, gap_extend: int) -> int:
    """
    Compute SSW alignment score.
    
    Parameters:
    -----------
    query_seq : str
        Query 3Di sequence (single or doubled)
    target_seq : str
        Target 3Di sequence (always single)
    is_x2 : bool
        If True, this is X2 (doubled query); if False, this is X1 (single)
    gap_open : int
        Gap opening penalty (default: 8)
    gap_extend : int
        Gap extension penalty (default: 2)
    
    Returns:
    --------
    int : Optimal alignment score from SSW
    """
    
    # Step 1: Write sequences to temporary FASTA files
    query_file = os.path.join(output_dir, "query_tmp.fasta")
    target_file = os.path.join(output_dir, "target_tmp.fasta")
    
    with open(query_file, "w") as qf:
        qf.write(f">{query_name}_{encoding}_3di\n{query_seq}")
    
    with open(target_file, "w") as tf:
        tf.write(f">{target_name}_{encoding}_3di\n{target_seq}")
    
    # Step 2: Run SSW binary with parameters
    cmd = [ssw_binary, '-o', str(gap_open), '-e', str(gap_extend), 
           '-p', '-c', query_file, target_file]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"SSW binary failed: {result.stderr}")
    
    # Step 3: Parse output to extract score
    # SSW output contains: "optimal_alignment_score: <score>"
    match = re.search(r"optimal_alignment_score:\s*(\d+)", result.stdout)
    
    if match:
        score = int(match.group(1))
        logging.info(f"Optimal alignment score: {score}")
        return score
    
    raise ValueError("Could not find alignment score in SSW output")
```

---

## 2. X2 Calculation (Doubled Query Sequence)

### Pseudo Code
```
INPUT:
  - query_seq: Original 3Di sequence (e.g., "VLLVVCCVPP...")
  - target_seq: Original 3Di sequence
  - matrix: Substitution matrix (8f or 10f)

PROCESS:
  1. Double the query sequence:
     query_seq_x2 = query_seq + query_seq
     Example: "VLLVVCCVPP..." → "VLLVVCCVPP...VLLVVCCVPP..."
  
  2. Create temporary FASTA files:
     - query_x2.fasta: ">query_name_x2\nVLLVVCCVPP...VLLVVCCVPP..."
     - target.fasta: ">target_name\nVLLVVCCVPP..."
  
  3. Run SSW with doubled query:
     ssw_test -o 8 -e 2 -a matrix.mat query_x2.fasta target.fasta
     
     Note: Target remains SINGLE sequence, only query is doubled
  
  4. Parse output to extract score:
     optimal_alignment_score: 501
  
  5. Return score_x2 = 501

OUTPUT:
  - score_x2: Integer (SSW optimal alignment score with doubled query)
```

... (content unchanged; moved from repo root)
