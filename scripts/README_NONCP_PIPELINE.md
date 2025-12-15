# Non-CP Homolog Pair Set Construction Pipeline

This pipeline constructs a **high-confidence Non-CP (non-circular permutation) homolog pair set** from SCOPe40 for benchmarking circular permutation detection with 3Di + SSW alignment.

## Overview

### Goal
Build a benchmark set of **~1,000–2,000 Non-CP homolog pairs** that:
- Share the **same structural fold** (SCOPe superfamily/SCCS)
- Have **high structural similarity** (TM-score ≥ 0.5)
- Are **unlikely to be CP** (|ΔCP| ≤ 0.02, where ΔCP = TM_cp − TM_nocp)
- Match the **length distribution of the CP-positive set**
- Exclude all domains that appear in the CP-positive set

### Method
1. **Stratified sampling**: Select query domains with length distribution matching CP-positive set
2. **Foldseek search**: Find homologs in the full SCOPe40 database
3. **Filtering**: Keep same-SF, different-species hits
4. **TM-align validation**: Run structural alignment with and without CP option
5. **Final filtering**: Apply TM-score and ΔCP thresholds
6. **Downsampling**: Balance to target size (~2,000 pairs) with max hits per query

## Input Files Required

### 1. `scope40_domains.tsv` (auto-generated)
TSV file with all SCOPe40 domains.

**Columns:**
- `domain_id`: Domain identifier (e.g., `d1jx6a_`)
- `length`: Number of amino acids
- `sccs`: SCOPe S/F code (e.g., `c.93.1.1`)
- `superfamily`: Superfamily code (e.g., `c.93.1`)
- `species_id`: Species identifier (proxy: PDB code)

**Generation:**
```bash
python3 scripts/create_scope40_domains_tsv.py
```

### 2. `input_data/datasets/cp_positive_pairs.tsv`
CP-positive pairs to exclude from query set.

**Columns:**
- `domA`, `domB`: Domain IDs

**Example:**
```
domA    domB    source
1qdm    2gtg    benchmark
1puj    1zbd    benchmark
```

### 3. `scope_pdb/` directory
PDB files for all SCOPe40 domains.

**Format:** `scope_pdb/<domain_id>.pdb`

Example: `scope_pdb/d1jx6a_.pdb`

### 4. External Tools
- **Foldseek**: For homology search
- **TM-align**: For structural validation

Ensure these are installed and in your `$PATH`.

## Pipeline Scripts

### Script 1: `sample_queries.py`
**Purpose**: Select query domains using stratified sampling over length bins.

**Key Features:**
- Excludes all CP-positive domains
- Stratifies by length bins: `<50`, `50-100`, `100-150`, `150-200`, `200-300`, `300-500`, `>500`
- Matches CP-positive distribution: 36.6% in 50–100aa, 35.0% in 100–150aa, etc.

**Usage:**
```bash
python3 scripts/sample_queries.py \
  --scope40_domains scope40_domains.tsv \
  --cp_positive_pairs input_data/datasets/cp_positive_pairs.tsv \
  --output selected_queries.tsv \
  --n_query 2000
```

**Output:** `selected_queries.tsv`

---

### Script 2: `build_foldseek_dbs.py`
**Purpose**: Build Foldseek databases and run homology search.

**Key Features:**
- Converts domains to FASTA (uses CA atoms from PDB)
- Builds query DB (from selected queries) and target DB (all ~19k domains)
- Runs `foldseek search` to find homologs
- Converts output to TSV with columns: `query, target, qlen, tlen, evalue, bits`

**Usage:**
```bash
python3 scripts/build_foldseek_dbs.py \
    --queries selected_queries.tsv \
    --scope40_domains scope40_domains.tsv \
    --pdb_dir scope_pdb \
    --work_dir foldseek_work
```

**Output:** `foldseek_work/foldseek_results.tsv`

**Timing**: 10–60 minutes depending on system.

---

### Script 3: `parse_foldseek_and_select_hits.py`
**Purpose**: Parse Foldseek results and select candidate Non-CP pairs.

**Filtering Logic:**
1. Keep only selected queries
2. Remove self-hits (query ≠ target)
3. Add metadata (SCCS, species, length) from `scope40_domains.tsv`
4. **Keep only same-SF** hits: `sccs(query) == sccs(target)`
5. Optionally **require different species**: `species_id(query) != species_id(target)`
6. Sort by bitscore (descending)
7. Select top **K hits per query** (default K=2)

**Usage:**
```bash
python3 scripts/parse_foldseek_and_select_hits.py \
    --foldseek_results foldseek_work/foldseek_results.tsv \
    --scope40_domains scope40_domains.tsv \
    --selected_queries selected_queries.tsv \
    --output candidate_noncp_pairs.tsv \
    --k_hits 2 \
    --require_different_species True
```

**Output:** `candidate_noncp_pairs.tsv`

**Columns:**
- `query_id`, `target_id`
- `query_sccs`, `target_sccs`
- `query_species`, `target_species`
- `query_length`, `target_length`
- `foldseek_bitscore`, `foldseek_evalue`

---

### Script 4: `run_tmalign_and_filter_noncp.py`
**Purpose**: Run TM-align structural validation and final filtering.

**Key Features:**
- For each pair:
  1. Run `TM-align A.pdb B.pdb` → get TM_nocp
  2. Run `TM-align A.pdb B.pdb -cp` → get TM_cp
  3. Compute ΔCP = TM_cp − TM_nocp
- Apply filters:
  - **TM_nocp ≥ 0.5** (fold-level similarity)
  - **|ΔCP| ≤ 0.02** (not benefiting from CP)
- Downsample to target size (default 2,000 pairs)
  - Random selection at pair level
  - Max 2 pairs per query (avoid query bias)
  - Approximate length distribution matching

**Usage:**
```bash
python3 scripts/run_tmalign_and_filter_noncp.py \
    --candidate_pairs candidate_noncp_pairs.tsv \
    --pdb_dir scope_pdb \
    --output_with_tm noncp_with_tm.tsv \
    --output_final noncp_homolog_pairs.tsv \
    --tm_threshold 0.5 \
    --delta_cp_threshold 0.02 \
    --n_final 2000 \
    --max_pairs_per_query 2
```

**Outputs:**
- `noncp_with_tm.tsv`: All pairs with TM-align results (for inspection)
- `noncp_homolog_pairs.tsv`: **Final filtered set** (ready for benchmarking)

**Timing**: 1–6 hours (depends on number of candidate pairs).

**Final Output Columns:**
- `query_id`, `target_id`
- `TM_nocp`, `TM_cp`, `ΔCP`
- `foldseek_bitscore`, `foldseek_evalue`
- `query_sccs`, `target_sccs`
- `query_species`, `target_species`
- `query_length`, `target_length`

---

## Running the Full Pipeline

### Option 1: Automated (Recommended)

```bash
bash scripts/run_noncp_pipeline.sh \
    --n_query 2000 \
    --n_final 2000 \
    --tm_threshold 0.5 \
    --delta_cp_threshold 0.02 \
    --k_hits 2 \
    --work_dir noncp_work
```

This runs all 4 scripts in sequence. Output files go to `noncp_work/`.

### Option 2: Manual Step-by-Step

```bash
# Step 0: Create scope40_domains.tsv (if not already done)
python3 scripts/create_scope40_domains_tsv.py

# Step 1: Sample queries
mkdir -p noncp_work
python3 scripts/sample_queries.py \
  --scope40_domains scope40_domains.tsv \
  --cp_positive_pairs input_data/datasets/cp_positive_pairs.tsv \
  --output noncp_work/selected_queries.tsv \
  --n_query 2000

# Step 2: Build Foldseek DBs and search
python3 scripts/build_foldseek_dbs.py \
    --queries noncp_work/selected_queries.tsv \
    --scope40_domains scope40_domains.tsv \
    --pdb_dir scope_pdb \
    --work_dir noncp_work

# Step 3: Parse Foldseek and select candidates
python3 scripts/parse_foldseek_and_select_hits.py \
    --foldseek_results noncp_work/foldseek_results.tsv \
    --scope40_domains scope40_domains.tsv \
    --selected_queries noncp_work/selected_queries.tsv \
    --output noncp_work/candidate_noncp_pairs.tsv \
    --k_hits 2

# Step 4: Run TM-align and filter
python3 scripts/run_tmalign_and_filter_noncp.py \
    --candidate_pairs noncp_work/candidate_noncp_pairs.tsv \
    --pdb_dir scope_pdb \
    --output_with_tm noncp_work/noncp_with_tm.tsv \
    --output_final noncp_work/noncp_homolog_pairs.tsv \
    --tm_threshold 0.5 \
    --delta_cp_threshold 0.02 \
    --n_final 2000
```

## Output Files

### Final Benchmark Set
**File:** `noncp_work/noncp_homolog_pairs.tsv`

This is the **main deliverable** for benchmarking.

**Statistics (example):**
- Total pairs: ~2,000
- Unique queries: ~1,000
- Mean TM_nocp: 0.68 (±0.12)
- Mean ΔCP: −0.001 (±0.015)
- Length distribution: Matches CP-positive set

### Intermediate Files (for debugging/inspection)
- `noncp_work/selected_queries.tsv`: 2,000 sampled query domains
- `noncp_work/candidate_noncp_pairs.tsv`: Foldseek hits after filtering
- `noncp_work/noncp_with_tm.tsv`: All pairs with TM-align results (before final thresholds)
- `noncp_work/foldseek_results.tsv`: Raw Foldseek output

## Key Parameters and Thresholds

### Stratified Sampling
Based on CP-positive set (738 domains):
- `<50aa`: 0.0%
- `50–100aa`: 36.6% ← most frequent
- `100–150aa`: 35.0%
- `150–200aa`: 9.9%
- `200–300aa`: 13.6%
- `300–500aa`: 5.0%
- `>500aa`: 0.0%

### Non-CP Filtering Thresholds
- **TM_threshold** (default: 0.5)
  - Standard cutoff for fold-level structural similarity
  - Adjust down to 0.4–0.45 for more lenient filtering
  
- **DELTA_CP_THRESHOLD** (default: 0.02)
  - A pair is "non-CP" if |ΔCP| ≤ 0.02
  - This means CP alignment gains <2% TM-score
  - Adjust up to 0.03–0.05 for stricter/looser definitions

### Query Selection
- **n_query** (default: 2,000)
  - Adjust based on compute budget
  - Typical range: 1,000–3,000
  
- **k_hits** (default: 2)
  - Maximum Foldseek hits per query retained
  - Adjust up to 3–5 for more candidates

### Final Downsampling
- **n_final** (default: 2,000)
  - Target size of final set
  - Adjust based on benchmarking needs
  
- **max_pairs_per_query** (default: 2)
  - Prevents single query from dominating final set
  - Adjust to 1–3 depending on diversity preference

## Expected Runtime

| Step | Typical Duration | Notes |
|------|-------------------|-------|
| Script 1 (Sampling) | <1 min | Fast |
| Script 2 (Foldseek) | 10–60 min | Depends on Foldseek DB size |
| Script 3 (Parsing) | <5 min | Fast |
| Script 4 (TM-align) | 1–6 hours | TM-align is slow; ~500 candidate pairs → 4–5 hours |
| **Total** | **2–7 hours** | Sequential execution |

## Troubleshooting

### Issue: `foldseek: command not found`
**Solution**: Install Foldseek or add to `$PATH`.
```bash
export PATH="/path/to/foldseek:$PATH"
```

### Issue: `TM-align: command not found`
**Solution**: Install TM-align or add to `$PATH`.
```bash
export PATH="/path/to/tmalign:$PATH"
```

### Issue: "Missing PDB files" in Script 1
**Solution**: Ensure `scope_pdb/` directory exists and contains all domain PDBs.

### Issue: Very few candidate pairs after Foldseek
**Possible causes:**
- Try increasing `k_hits` (e.g., 3–5 instead of 2)
- Relax `tm_threshold` or `delta_cp_threshold` in Script 4
- Check that Foldseek search ran successfully

## Citation & References

This pipeline is designed for benchmarking CP detection using:
- **3Di encoding**: Foldseek protein language model
- **SSW alignment**: Smith-Waterman for 3Di comparison
- **TM-align validation**: Structural similarity assessment

Key papers:
- SCOPe40: [Citation]
- CIRPIN CP detection: [Citation]
- Foldseek: [Citation]

## Contact

For questions or issues, please refer to the repository or contact the maintainers.
