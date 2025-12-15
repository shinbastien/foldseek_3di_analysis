#!/bin/bash

# SSW test script for basic (non-x2) 3di.fasta files
# This script runs ssw_test on all 4 protein pairs with 3 variants each

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=========================================="
echo "SSW Test for Basic 3di FASTA Files"
echo "=========================================="
echo ""

# Define pairs
pairs=(
    "2hga_trim:2vsv_PDZ"
    "2vsv_PDZ:2z9i_pdz"
    "2z9i_pdz:2hga_trim"
    "Sap_C:Sap_C_circular_permutation"
)

variants=("8f" "9f" "10f")
w=0  # Using w=0 for consistency

output_dir="ssw_basic_3di_results"
mkdir -p "$output_dir"

total_runs=$((${#pairs[@]} * ${#variants[@]}))
current_run=0

for pair_str in "${pairs[@]}"; do
    IFS=':' read -r query_name target_name <<< "$pair_str"
    
    echo -e "${YELLOW}=====================================${NC}"
    echo -e "${YELLOW}Pair: ${query_name} vs ${target_name}${NC}"
    echo -e "${YELLOW}=====================================${NC}"
    
    for variant in "${variants[@]}"; do
        current_run=$((current_run + 1))
        
        # Construct file paths
        query_fasta="tmp/${query_name}_vs_${target_name}_w${w}/${query_name}/${query_name}_${variant}_3di.fasta"
        target_fasta="tmp/${query_name}_vs_${target_name}_w${w}/${target_name}/${target_name}_${variant}_3di.fasta"
        matrix="ssw/s_${variant}.mat"
        
        # Check if files exist
        if [[ ! -f "$query_fasta" ]]; then
            echo -e "${RED}[${current_run}/${total_runs}] SKIP: Query file not found: $query_fasta${NC}"
            continue
        fi
        if [[ ! -f "$target_fasta" ]]; then
            echo -e "${RED}[${current_run}/${total_runs}] SKIP: Target file not found: $target_fasta${NC}"
            continue
        fi
        if [[ ! -f "$matrix" ]]; then
            echo -e "${RED}[${current_run}/${total_runs}] SKIP: Matrix not found: $matrix${NC}"
            continue
        fi
        
        # Output file
        output_file="${output_dir}/${query_name}_vs_${target_name}_${variant}_w${w}.txt"
        
        echo -e "${GREEN}[${current_run}/${total_runs}] Running: ${query_name} vs ${target_name} (${variant})${NC}"
        echo "  Query:  $query_fasta"
        echo "  Target: $target_fasta"
        echo "  Matrix: $matrix"
        echo "  Output: $output_file"
        
        # Run ssw_test
        ./ssw/ssw_test -c "$query_fasta" "$target_fasta" "$matrix" > "$output_file" 2>&1
        
        if [[ $? -eq 0 ]]; then
            # Extract score from output
            score=$(grep -oP 'optimal_alignment_score:\s+\K\d+' "$output_file" | head -1)
            if [[ -n "$score" ]]; then
                echo -e "  ${GREEN}✓ Score: $score${NC}"
            else
                echo -e "  ${YELLOW}✓ Done (score not found in output)${NC}"
            fi
        else
            echo -e "  ${RED}✗ Failed${NC}"
        fi
        echo ""
    done
done


echo ""
echo "=========================================="
echo "Summary: All results saved to $output_dir/"
echo "=========================================="
echo ""
echo "To view all scores:"
echo "  grep 'optimal_alignment_score:' $output_dir/*.txt"
echo ""
echo "To create a summary table:"
echo "  for f in $output_dir/*.txt; do echo -n \"$(basename $f .txt): \"; grep 'optimal_alignment_score:' $f | head -1; done"
