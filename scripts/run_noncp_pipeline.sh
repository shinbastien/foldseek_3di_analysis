#!/bin/bash
# 
# Non-CP Homolog Pair Set Construction Pipeline
# 
# This script runs the full pipeline to construct a high-confidence Non-CP homolog
# set for benchmarking circular permutation (CP) detection.
#
# Prerequisites:
#   - scope40_domains.tsv (all SCOPe40 domains with metadata)
#   - input_data/datasets/cp_positive_pairs.tsv (CP-positive pairs to exclude)
#   - scope_pdb/ directory (PDB files)
#   - Foldseek installed and in PATH
#   - TM-align installed and in PATH
#
# Usage:
#   bash run_noncp_pipeline.sh [options]
#
# Options:
#   --n_query N          Number of query domains (default: 2000)
#   --n_final N          Target number of final pairs (default: 2000)
#   --tm_threshold T     Minimum TM_nocp (default: 0.5)
#   --delta_cp_threshold T  Maximum |ΔCP| (default: 0.02)
#   --k_hits K           Max Foldseek hits per query (default: 2)
#   --work_dir DIR       Working directory (default: noncp_work)
#   --help               Show this message
#

set -e  # Exit on error

# Default parameters
N_QUERY=2000
N_FINAL=2000
TM_THRESHOLD=0.5
DELTA_CP_THRESHOLD=0.02
K_HITS=2
WORK_DIR="noncp_work"
SCOPE40_DOMAINS="scope40_domains.tsv"
CP_POSITIVE_PAIRS="input_data/datasets/cp_positive_pairs.tsv"
PDB_DIR="scope_pdb"
SCRIPTS_DIR="scripts"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --n_query)
            N_QUERY=$2
            shift 2
            ;;
        --n_final)
            N_FINAL=$2
            shift 2
            ;;
        --tm_threshold)
            TM_THRESHOLD=$2
            shift 2
            ;;
        --delta_cp_threshold)
            DELTA_CP_THRESHOLD=$2
            shift 2
            ;;
        --k_hits)
            K_HITS=$2
            shift 2
            ;;
        --work_dir)
            WORK_DIR=$2
            shift 2
            ;;
        --help)
            grep "^#" "$0" | grep -v "^#!/bin/bash" | sed 's/^# //'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create working directory
mkdir -p "$WORK_DIR"

echo "================================================================================"
echo "Non-CP Homolog Pair Set Construction Pipeline"
echo "================================================================================"
echo ""
echo "Parameters:"
echo "  N_query: $N_QUERY"
echo "  N_final: $N_FINAL"
echo "  TM_threshold: $TM_THRESHOLD"
echo "  Delta_CP_threshold: $DELTA_CP_THRESHOLD"
echo "  K_hits (per query): $K_HITS"
echo "  Work directory: $WORK_DIR"
echo ""

# Step 0: Create scope40_domains.tsv if not exists
if [ ! -f "$SCOPE40_DOMAINS" ]; then
    echo "Step 0: Creating scope40_domains.tsv..."
    python3 "$SCRIPTS_DIR/create_scope40_domains_tsv.py"
    echo ""
fi

# Step 1: Sample query domains
echo "Step 1: Sampling query domains with stratified sampling..."
SELECTED_QUERIES="$WORK_DIR/selected_queries.tsv"
python3 "$SCRIPTS_DIR/sample_queries.py" \
    --scope40_domains "$SCOPE40_DOMAINS" \
    --cp_positive_pairs "$CP_POSITIVE_PAIRS" \
    --output "$SELECTED_QUERIES" \
    --n_query "$N_QUERY"
echo ""

# Step 2: Build Foldseek databases and run search
echo "Step 2: Building Foldseek databases and running search..."
echo "(This may take 10-60 minutes depending on database size)"
FOLDSEEK_RESULTS="$WORK_DIR/foldseek_results.tsv"
python3 "$SCRIPTS_DIR/build_foldseek_dbs.py" \
    --queries "$SELECTED_QUERIES" \
    --scope40_domains "$SCOPE40_DOMAINS" \
    --pdb_dir "$PDB_DIR" \
    --work_dir "$WORK_DIR"
echo ""

# Step 3: Parse Foldseek results and select candidate pairs
echo "Step 3: Parsing Foldseek results and selecting candidate pairs..."
CANDIDATE_PAIRS="$WORK_DIR/candidate_noncp_pairs.tsv"
python3 "$SCRIPTS_DIR/parse_foldseek_and_select_hits.py" \
    --foldseek_results "$FOLDSEEK_RESULTS" \
    --scope40_domains "$SCOPE40_DOMAINS" \
    --selected_queries "$SELECTED_QUERIES" \
    --output "$CANDIDATE_PAIRS" \
    --k_hits "$K_HITS"
echo ""

# Step 4: Run TM-align and filter
echo "Step 4: Running TM-align and filtering Non-CP homolog pairs..."
echo "(This may take 1-6 hours depending on number of pairs)"
NONCP_WITH_TM="$WORK_DIR/noncp_with_tm.tsv"
NONCP_FINAL="$WORK_DIR/noncp_homolog_pairs.tsv"
python3 "$SCRIPTS_DIR/run_tmalign_and_filter_noncp.py" \
    --candidate_pairs "$CANDIDATE_PAIRS" \
    --pdb_dir "$PDB_DIR" \
    --output_with_tm "$NONCP_WITH_TM" \
    --output_final "$NONCP_FINAL" \
    --tm_threshold "$TM_THRESHOLD" \
    --delta_cp_threshold "$DELTA_CP_THRESHOLD" \
    --n_final "$N_FINAL"
echo ""

echo "================================================================================"
echo "Pipeline Completed Successfully!"
echo "================================================================================"
echo ""
echo "Output files:"
echo "  Intermediate results:"
echo "    $SELECTED_QUERIES          - Selected query domains"
echo "    $CANDIDATE_PAIRS           - Candidate pairs from Foldseek"
echo "    $NONCP_WITH_TM             - All pairs with TM-align scores"
echo ""
echo "  Final results:"
echo "    $NONCP_FINAL               - FINAL Non-CP homolog pairs (ready for benchmarking)"
echo ""
echo "Next steps:"
echo "  1. Review the final pairs in $NONCP_FINAL"
echo "  2. Compare CP-positive and Non-CP sets using 3Di-SSW analysis"
echo "================================================================================"
