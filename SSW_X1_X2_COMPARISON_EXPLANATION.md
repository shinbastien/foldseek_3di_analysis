Moved to `experiments/02_8f_10f_ssw_comparison/SSW_X1_X2_COMPARISON_EXPLANATION.md`.

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

### Real Python Code
```python
# From process_pair() in compare_8f_10f_ssw.py

# First, get sequences for given encoding
query_seq, query_seq_x2, query_output, query_output_x2 = get_3di_sequence(
    query, pdb_dir, encoding, MODEL_DIR, OUTPUT_DIR
)
target_seq, target_seq_x2, target_output, target_output_x2 = get_3di_sequence(
    target, pdb_dir, encoding, MODEL_DIR, OUTPUT_DIR
)

# Key: Notice that X2 sequences are created by doubling
# query_seq_x2 = query_seq + query_seq  (doubled)
# target_seq_x2 = target_seq + target_seq  (also doubled, but we use single target)

# Compute X1 score (single sequences)
logging.info(f"Computing SSW X1 score...")
score_x1 = compute_ssw_score_with_binary(
    ssw_binary, 
    query_seq,        # ← SINGLE sequence
    target_seq,       # ← SINGLE sequence
    query, target, encoding, 
    is_x2=False,      # Flag indicating this is X1
    output_dir=OUTPUT_DIR, 
    gap_open=8, gap_extend=2
)

# Compute X2 score (doubled query, single target)
logging.info(f"Computing SSW X2 score...")
score_x2 = compute_ssw_score_with_binary(
    ssw_binary, 
    query_seq_x2,     # ← DOUBLED sequence (query_seq + query_seq)
    target_seq,       # ← SINGLE sequence (not doubled!)
    query, target, encoding, 
    is_x2=True,       # Flag indicating this is X2
    output_dir=OUTPUT_DIR, 
    gap_open=8, gap_extend=2
)

# Calculate difference
score_diff = score_x2 - score_x1  # ← Delta SSW = X2 - X1
```

---

## 3. Delta SSW Calculation

### Pseudo Code
```
INPUT:
  - score_x1: Integer (single sequence alignment score)
  - score_x2: Integer (doubled query alignment score)

PROCESS:
  score_diff = score_x2 - score_x1

OUTPUT:
  - score_diff: Integer (positive = X2 better, negative = X1 better, 0 = no difference)
```

### Real Python Code
```python
# From process_pair() in compare_8f_10f_ssw.py

# After computing both scores
score_diff = score_x2 - score_x1

logging.info(f"\nSummary for {query} vs {target} ({encoding}):")
logging.info(f"  Score X1: {score_x1}")
logging.info(f"  Score X2: {score_x2}")
logging.info(f"  Difference (X2 - X1): {score_diff}")

# Store in results dataframe
results.append({
    "query": query,
    "target": target,
    "encoding": encoding,
    "score_x1": score_x1,
    "score_x2": score_x2,
    "score_diff": score_diff
})
```

---

## 4. Complete Example with Real Data

### Input Data
```
Query protein:  d1byka_  (SCOP domain)
Target protein: d8abpa_  (SCOP domain)
Encoding:       10f (Foldseek official)

Query 3Di sequence (example):
  DQLLVLLLQLLVVLVVCVVVPDDLVRQLVCLVCSLVVDDPVRSVSNVVCCVPPVSVSNVVVVVVDDSVCSCVVVPSD

Target 3Di sequence:
  VLLVVCCVPPVLVLLVVVVVVDDLVCSCVVVPLDCPRSSVSSNQSSVVSVVCVVVPDDSVRRNVCSVCSSV
```

### X1 Calculation
```
Step 1: Create FASTA files
  query.fasta:
    >d1byka_10f_3di
    DQLLVLLLQLLVVLVVCVVVPDDLVRQLVCLVCSLVVDDPVRSVSNVVCCVPPVSVSNVVVVVVDDSVCSCVVVPSD
  
  target.fasta:
    >d8abpa_10f_3di
    VLLVVCCVPPVLVLLVVVVVVDDLVCSCVVVPLDCPRSSVSSNQSSVVSVVCVVVPDDSVRRNVCSVCSSV

Step 2: Run SSW
  ssw_test -o 8 -e 2 -a s_10f.mat query.fasta target.fasta

Step 3: SSW Output (simplified)
  target_name: d8abpa_10f_3di
  query_name: d1byka_10f_3di
  optimal_alignment_score: 480
  suboptimal_alignment_score: 15
  
  Alignment visualization:
  DQLLVLLLQLLVVLVVCVVVPDDLVRQLVCLVCSLVVDDPVRSVSNVVCCVPPVSVSNVVVVVVDDSVCSCVVVPSD
  |  || |||| | || ||||||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||
  VLLVVCCVPPVLVLLVVVVVVDDLVCSCVVVPLDCPRSSVSSNQSSVVSVVCVVVPDDSVRRNVCSVCSSV

Step 4: Extract Score
  score_x1 = 480
```

### X2 Calculation
```
Step 1: Double the query sequence
  query_x2 = DQLLVLLLQLLVVLVVCVVVPDDLVRQLVCLVCSLVVDDPVRSVSNVVCCVPPVSVSNVVVVVVDDSVCSCVVVPSD
           + DQLLVLLLQLLVVLVVCVVVPDDLVRQLVCLVCSLVVDDPVRSVSNVVCCVPPVSVSNVVVVVVDDSVCSCVVVPSD
           
           = DQLLVLLLQLLVVLVVCVVVPDDLVRQLVCLVCSLVVDDPVRSVSNVVCCVPPVSVSNVVVVVVDDSVCSCVVVPSDDQLLVLLLQLLVVLVVCVVVPDDLVRQLVCLVCSLVVDDPVRSVSNVVCCVPPVSVSNVVVVVVDDSVCSCVVVPSD
  
  (Length doubled: 79 → 158 characters)

Step 2: Create FASTA files
  query_x2.fasta:
    >d1byka_10f_3di_x2
    DQLLVLLLQLLVVLVVCVVVPDDLVRQLVCLVCSLVVDDPVRSVSNVVCCVPPVSVSNVVVVVVDDSVCSCVVVPSDDQLLVLLLQLLVVLVVCVVVPDDLVRQLVCLVCSLVVDDPVRSVSNVVCCVPPVSVSNVVVVVVDDSVCSCVVVPSD
  
  target.fasta:
    >d8abpa_10f_3di
    VLLVVCCVPPVLVLLVVVVVVDDLVCSCVVVPLDCPRSSVSSNQSSVVSVVCVVVPDDSVRRNVCSVCSSV
    
  (Target remains SINGLE)

Step 3: Run SSW
  ssw_test -o 8 -e 2 -a s_10f.mat query_x2.fasta target.fasta

Step 4: SSW Output
  optimal_alignment_score: 501
  
  Note: Higher score (501 > 480) because SSW can now find alignment
  at different starting positions in the doubled query

Step 5: Extract Score
  score_x2 = 501
```

### Delta SSW
```
score_diff = score_x2 - score_x1
           = 501 - 480
           = 21 points

Interpretation:
  - Positive difference (21) indicates X2 is better
  - This suggests the query might have a circular permutation
  - The doubled query allows SSW to find a better alignment at a different starting point
```

---

## 5. Interpretation of Results

### Case 1: No Circular Permutation
```
Query: d3d02a_
Target: d8abpa_
Encoding: 10f

score_x1 = 691
score_x2 = 691
score_diff = 0

→ X1 and X2 scores are identical
→ No circular permutation detected
→ Both alignments equally good
```

### Case 2: Circular Permutation Detected
```
Query: d1byka_
Target: d8abpa_
Encoding: 10f

score_x1 = 480
score_x2 = 501
score_diff = 21

→ X2 score is 21 points higher than X1
→ Circular permutation likely present
→ SSW found better alignment when query was doubled
```

### Case 3: Different between 8f and 10f
```
Query: d1jyea_
Target: d8abpa_

8f encoding:
  score_x1 = 204
  score_x2 = 204
  score_diff = 0

10f encoding:
  score_x1 = 638
  score_x2 = 640
  score_diff = 2

→ 10f alphabet gives higher absolute scores
→ But circular permutation effect (diff=2) is subtle in 10f
```

---

## 6. How This Works (Biological Insight)

Circular permutations are proteins with the same fold but different starting points in the sequence:

```
Original:     A---B---C---D---E---F---A (circular protein)

Same protein, different starting point:
              D---E---F---A---B---C---D

As 1D sequence:
  Original:   ABCDEFABCDEFABCDEF...
  Permuted:   DEFABCDEFABCDEFABC...

When we DOUBLE the original:
  Original×2: ABCDEFABCDEFABCDEFABCDEFABCDEFABCDEF...
                         ↑ permutation alignment starts here

This allows SSW to find the permuted region in the doubled sequence,
resulting in higher alignment score.
```

---

## 7. Code Flow Summary

```
START
  ↓
FOR EACH (query, target, encoding) PAIR:
  ├─ Load 3Di sequences
  │  ├─ query_seq = get_3di_for_encoding(query, encoding)
  │  ├─ query_seq_x2 = query_seq + query_seq  ← Double it
  │  └─ target_seq = get_3di_for_encoding(target, encoding)
  │
  ├─ COMPUTE X1 (single sequence)
  │  ├─ Run SSW(query_seq, target_seq)
  │  └─ score_x1 = extract_score(ssw_output)
  │
  ├─ COMPUTE X2 (doubled query)
  │  ├─ Run SSW(query_seq_x2, target_seq)  ← Doubled query!
  │  └─ score_x2 = extract_score(ssw_output)
  │
  ├─ COMPUTE DELTA
  │  └─ score_diff = score_x2 - score_x1
  │
  └─ SAVE RESULTS
     └─ {query, target, encoding, score_x1, score_x2, score_diff}

END
```

---

## 8. Key Parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| Gap Open | 8 | Penalty for opening a gap |
| Gap Extend | 2 | Penalty for extending a gap |
| Matrix | s_8f.mat or s_10f.mat | Amino acid/3Di substitution matrix |
| Query (X1) | Original | Single 3Di sequence |
| Query (X2) | Doubled | Concatenated sequence (seq + seq) |
| Target | Original | Always single (never doubled) |

---

## References

- **SSW Library**: Smith-Waterman with SIMD optimization
- **3Di**: 3-letter structure alphabet from Foldseek
- **Circular Permutation**: Protein with same fold but different sequence order

