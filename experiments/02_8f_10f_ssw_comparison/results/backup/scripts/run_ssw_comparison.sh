#!/bin/bash -x

# Usage:
#   ./run-benchmark.sh /path/to/s_10f.mat data/pdbs_val.txt data/scop_lookup.tsv
# Requires:
#   - foldseek in PATH (e.g., conda activate torchclean)
#   - ssw_test built (we will build if missing like 8f flow)

SUBMAT_SRC=$1
PDBS=$2           # list of SIDs (validation set)
SCOPLOOKUP=$3

if [ "$#" -ne 3 ]; then
  echo "Illegal number of parameters"
  echo "Usage: $0 /path/to/s_10f.mat data/pdbs_val.txt data/scop_lookup.tsv"
  exit 2
fi

# Sanity
if [ ! -f "$SUBMAT_SRC" ]; then
  echo "substitution matrix not found: $SUBMAT_SRC" >&2
  exit 3
fi
if [ ! -f "$PDBS" ]; then
  echo "pdbs list not found: $PDBS" >&2
  exit 3
fi
if [ ! -f "$SCOPLOOKUP" ]; then
  echo "scop_lookup.tsv not found: $SCOPLOOKUP" >&2
  exit 3
fi

PDB_DIR=../scope_pdb
mkdir -p tmp tmp/splits tmp/alignments tmp/fa tmp/logs

# Ensure foldseek is available (expect conda activate torchclean outside)
FOLDSEEK_BIN=$(command -v foldseek || true)
if [ -z "$FOLDSEEK_BIN" ]; then
  # try known path
  if [ -x /home/jugipalace/miniconda3/envs/torchclean/bin/foldseek ]; then
    FOLDSEEK_BIN=/home/jugipalace/miniconda3/envs/torchclean/bin/foldseek
  else
    echo "ERROR: foldseek not found in PATH. Please 'conda activate torchclean' first." >&2
    exit 4
  fi
fi
export FOLDSEEK_BIN

# Build ssw_test if missing (same as 8f flow)
if [ ! -f tmp/ssw_test ]; then
    git clone --depth 1 https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library tmp/ssw || true
    (cd tmp/ssw/src && make)
    cp tmp/ssw/src/ssw_test tmp/ssw_test
fi

# Filter scop_lookup.tsv to validation set
awk 'FNR==NR {pdbs[$1]=1; next} ($1 in pdbs) {print $0}' \
    "$PDBS" "$SCOPLOOKUP" > tmp/scop_lookup_filtered.tsv

# Copy submat to expected name
cp "$SUBMAT_SRC" tmp/sub_score.mat

# Build 3Di target.fasta using foldseek createdb per PDB
# Skip if already extracted (resume-friendly)
while read -r sid; do
  [ -z "$sid" ] && continue
  pdb_path="$PDB_DIR/${sid}.pdb"
  if [ ! -f "$pdb_path" ]; then
    echo "WARN: missing PDB for $sid at $pdb_path, skipping" >&2
    continue
  fi
  out_fa="tmp/fa/${sid}.fasta"
  
  # Skip if already extracted
  if [ -f "$out_fa" ]; then
    echo "SKIP: $sid (already extracted)"
    continue
  fi
  
  python3 ./extract_3di_from_createdb.py "$pdb_path" "$out_fa" \
    > tmp/logs/${sid}.log 2>&1 || { echo "WARN: failed 3Di for $sid" >&2; continue; }
  echo "DONE: $sid"

done < "$PDBS"

# Now rebuild target.fasta from all extracted fa/*.fasta
echo "Rebuilding target.fasta from all fa/*.fasta..."
cat tmp/fa/*.fasta > tmp/target.fasta
echo "target.fasta regenerated with $(grep -c '^>' tmp/target.fasta) sequences"

# Split FASTA into 30 chunks like 8f flow
split -n 30 -d tmp/target.fasta tmp/splits/split_ --additional-suffix=.fasta

# Reuse existing smith-waterman runner (expects tmp/sub_score.mat and tmp/target.fasta)
../training_3di_gpu_8f/run-smithwaterman.sh 8 2

# ROC and AUC (reuse existing awk)
../training_3di_gpu_8f/roc1.awk tmp/scop_lookup_filtered.tsv \
    <(cat tmp/alignments/*.m8) > tmp/result.rocx

awk '{famsum+=$3; supfamsum+=$4; foldsum+=$5} END{print famsum/NR,supfamsum/NR,foldsum/NR}' \
    tmp/result.rocx
