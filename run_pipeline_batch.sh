#!/usr/bin/env bash
set -euo pipefail

# Sequential foreground runner for tmalign_3di_match_pipeline.py
# Creates tmp/<p1>_vs_<p2>_w<w> directories and writes runner.out/runner.err there.

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP_DIR="$ROOT_DIR/tmp"
PYTHON="python3"
PIPELINE="$ROOT_DIR/tmalign_3di_match_pipeline.py"

mkdir -p "$TMP_DIR"

# Start from a clean tmp state per user request: remove any existing tmp contents so runs start fresh
if [ -d "$TMP_DIR" ]; then
	echo "Cleaning existing tmp directory: $TMP_DIR"
	rm -rf "$TMP_DIR"/* || true
fi
mkdir -p "$TMP_DIR"

PAIRS=(
	# Use original PDB files (user prefers providing pdb inputs rather than .3di text files)
	"$ROOT_DIR/permutation_examples/2hga_trim.pdb:$ROOT_DIR/permutation_examples/2vsv_PDZ.pdb"
	"$ROOT_DIR/permutation_examples/2vsv_PDZ.pdb:$ROOT_DIR/permutation_examples/2z9i_pdz.pdb"
	"$ROOT_DIR/permutation_examples/2z9i_pdz.pdb:$ROOT_DIR/permutation_examples/2hga_trim.pdb"
	"$ROOT_DIR/permutation_examples/Sap_C.pdb:$ROOT_DIR/permutation_examples/Sap_C_circular_permutation.pdb"
)

WS=(0 1 2)

echo "Running pipeline sequentially for ${#PAIRS[@]} pairs × ${#WS[@]} w values..."

for pair in "${PAIRS[@]}"; do
	IFS=':' read -r p1 p2 <<< "$pair"
		# compute stems and strip common extensions so batch and pipeline agree on run name
		p1stem_raw="$(basename "$p1")"
		p2stem_raw="$(basename "$p2")"
		# strip final extension (e.g. .3di.txt, .txt, .pdb) by removing suffix after last dot
		p1stem="${p1stem_raw%.*}"
		p2stem="${p2stem_raw%.*}"

	for w in "${WS[@]}"; do
		rundir="$TMP_DIR/${p1stem}_vs_${p2stem}_w${w}"
		# Make a placeholder logs dir so we can store stdout/stderr there; do not pre-create
		# the pipeline's tmpdir (pipeline will create "$TMP_DIR/${p1stem}_vs_${p2stem}" itself).
		mkdir -p "$rundir"/logs
		echo "==> Running: $p1stem vs $p2stem  (w=$w)  -> $rundir"

		out="$rundir/logs/runner.out"
		err="$rundir/logs/runner.err"

		# Run pipeline in foreground, capture stdout/stderr to files while also streaming to terminal
		( 
			echo "START: $(date)" > "$out"
			echo "START: $(date)" > "$err"
			# Pass an explicit --run-dir so the pipeline will write into our expected rundir
			"$PYTHON" "$PIPELINE" "$p1" "$p2" --w-pair "$w" --run-dir "$rundir" >> "$out" 2>> "$err"
			echo "EXIT: $?  $(date)" >> "$out"
		) &

		pid=$!
		echo "Started PID $pid for ${p1stem}_vs_${p2stem}_w${w}; waiting for it to finish..."
		wait $pid || echo "Warning: run exited with non-zero status (pair ${p1stem}_vs_${p2stem}, w=$w)" >> "$err"

		# Find the tmp directory the pipeline actually created (may be without the _w suffix)
		created_tmp="$(ls -td "$TMP_DIR/${p1stem}_vs_${p2stem}"* 2>/dev/null | head -n1 || true)"
		if [ -n "$created_tmp" ]; then
			# if pipeline created directory different from our rundir, move it under the rundir
			if [ "$created_tmp" != "$rundir" ]; then
				# move contents from created_tmp into the expected rundir (preserve rundir placeholder)
				mv "$created_tmp"/* "$rundir"/ 2>/dev/null || true
				# cleanup the now-empty created_tmp directory if possible
				rmdir "$created_tmp" 2>/dev/null || true
			fi
		fi
		echo "Finished run ${p1stem}_vs_${p2stem}_w${w}"
		echo
	done
done

echo "All runs attempted. Check $TMP_DIR for per-run runner.out/runner.err and outputs."
