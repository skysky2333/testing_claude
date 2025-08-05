#!/bin/bash

# Rosetta Protein Folding and Mutation Analysis Pipeline
# Author: Claude Code Assistant
# Description: Complete pipeline for ab initio folding, relaxation, and ΔΔG calculation

set -euo pipefail

# Configuration variables
ROSETTA_PATH="/usr/local/rosetta.source.release-371"
ROSETTA_BIN="${ROSETTA_PATH}/main/source/bin"
ROSETTA_DB="${ROSETTA_PATH}/main/database"

# Input files
PROTEIN_FASTA="protein.faa"
MUTATIONS_FILE="mutations.txt"

# Output directories
OUTPUT_DIR="rosetta_results"
FRAGMENTS_DIR="${OUTPUT_DIR}/fragments"
FOLDING_DIR="${OUTPUT_DIR}/folding"
RELAX_DIR="${OUTPUT_DIR}/relax"
DDG_DIR="${OUTPUT_DIR}/ddg"
IMAGES_DIR="${OUTPUT_DIR}/images"

# Parameters
NSTRUCT_FOLDING=1000
NSTRUCT_RELAX=10
DDG_ITERATIONS=50

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
    exit 1
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

# Check dependencies
check_dependencies() {
    log "Checking dependencies..."
    
    if [[ ! -d "$ROSETTA_PATH" ]]; then
        error "Rosetta installation not found at $ROSETTA_PATH"
    fi
    
    if [[ ! -f "$PROTEIN_FASTA" ]]; then
        error "Protein FASTA file not found: $PROTEIN_FASTA"
    fi
    
    if [[ ! -f "$MUTATIONS_FILE" ]]; then
        error "Mutations file not found: $MUTATIONS_FILE"
    fi
    
    # Check for required Rosetta executables
    local required_exes=(
        "fragment_picker.linuxgccrelease"
        "AbinitioRelax.linuxgccrelease" 
        "relax.linuxgccrelease"
        "ddg_monomer.linuxgccrelease"
        "extract_pdbs.linuxgccrelease"
    )
    
    for exe in "${required_exes[@]}"; do
        if [[ ! -f "${ROSETTA_BIN}/${exe}" ]]; then
            error "Required Rosetta executable not found: ${exe}"
        fi
    done
    
    # Check for PyMOL
    if ! command -v pymol &> /dev/null; then
        warning "PyMOL not found. Structure visualization will be skipped."
    fi
    
    success "All dependencies checked successfully"
}

# Create directory structure
setup_directories() {
    log "Setting up directory structure..."
    
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$FRAGMENTS_DIR"
    mkdir -p "$FOLDING_DIR" 
    mkdir -p "$RELAX_DIR"
    mkdir -p "$DDG_DIR"
    mkdir -p "$IMAGES_DIR"
    
    success "Directory structure created"
}

# Extract protein name from FASTA
get_protein_name() {
    local protein_name=$(grep "^>" "$PROTEIN_FASTA" | head -1 | sed 's/^>//' | tr ' ' '_' | cut -d'|' -f1)
    echo "${protein_name:-protein}"
}

# Generate fragment files
generate_fragments() {
    log "Generating fragment files..."
    
    local protein_name=$(get_protein_name)
    
    # Create fragment picker configuration
    cat > "${FRAGMENTS_DIR}/fragment_picker_config" << EOF
# Fragment picker configuration
database ${ROSETTA_DB}
in:file:fasta ${PROTEIN_FASTA}
frags:describe_fragments ${FRAGMENTS_DIR}/${protein_name}_frags.fsc
out:file:frag_prefix ${FRAGMENTS_DIR}/${protein_name}
frags:n_candidates 200
frags:n_frags 25
frags:frag_sizes 3 9
EOF

    cd "$FRAGMENTS_DIR"
    
    "${ROSETTA_BIN}/fragment_picker.linuxgccrelease" \
        @fragment_picker_config \
        > fragment_picker.log 2>&1
    
    cd - > /dev/null
    
    if [[ -f "${FRAGMENTS_DIR}/${protein_name}.200.3mers" && -f "${FRAGMENTS_DIR}/${protein_name}.200.9mers" ]]; then
        success "Fragment files generated successfully"
    else
        error "Fragment generation failed. Check ${FRAGMENTS_DIR}/fragment_picker.log"
    fi
}

# Perform ab initio folding
run_ab_initio_folding() {
    log "Running ab initio folding..."
    
    local protein_name=$(get_protein_name)
    
    cd "$FOLDING_DIR"
    
    "${ROSETTA_BIN}/AbinitioRelax.linuxgccrelease" \
        -database "$ROSETTA_DB" \
        -in:file:fasta "../${PROTEIN_FASTA}" \
        -in:file:frag3 "../${FRAGMENTS_DIR}/${protein_name}.200.3mers" \
        -in:file:frag9 "../${FRAGMENTS_DIR}/${protein_name}.200.9mers" \
        -abinitio:relax \
        -relax:fast \
        -nstruct "$NSTRUCT_FOLDING" \
        -out:file:silent "folding_results.out" \
        -out:nooutput \
        > folding.log 2>&1
    
    cd - > /dev/null
    
    if [[ -f "${FOLDING_DIR}/folding_results.out" ]]; then
        success "Ab initio folding completed"
    else
        error "Ab initio folding failed. Check ${FOLDING_DIR}/folding.log"
    fi
}

# Select best folding model
select_best_model() {
    log "Selecting best folding model..."
    
    cd "$FOLDING_DIR"
    
    # Extract scores and find best model
    grep "^SCORE:" folding_results.out | sort -k2,2n | head -1 > best_model_score.txt
    local best_tag=$(awk '{print $NF}' best_model_score.txt)
    
    # Extract best model as PDB
    "${ROSETTA_BIN}/extract_pdbs.linuxgccrelease" \
        -in:file:silent folding_results.out \
        -in:file:tags <(echo "$best_tag") \
        -out:prefix "best_"
    
    mv "best_${best_tag}.pdb" best_model.pdb
    
    cd - > /dev/null
    
    success "Best model selected: $best_tag"
}

# Relax the best model
relax_structure() {
    log "Relaxing structure..."
    
    cd "$RELAX_DIR"
    
    "${ROSETTA_BIN}/relax.linuxgccrelease" \
        -database "$ROSETTA_DB" \
        -in:file:s "../${FOLDING_DIR}/best_model.pdb" \
        -relax:fast \
        -relax:constrain_relax_to_start_coords \
        -relax:ramp_constraints false \
        -ex1 \
        -ex2 \
        -use_input_sc \
        -flip_HNQ \
        -no_optH false \
        -nstruct "$NSTRUCT_RELAX" \
        -out:prefix "relaxed_" \
        -out:suffix "_wt" \
        > relax.log 2>&1
    
    cd - > /dev/null
    
    # Select best relaxed structure
    cd "$RELAX_DIR"
    local best_relaxed=$(ls relaxed_*_wt_*.pdb | head -1)
    cp "$best_relaxed" "relaxed_wildtype.pdb"
    cd - > /dev/null
    
    success "Structure relaxation completed"
}

# Parse mutations and create mutation files
parse_mutations() {
    log "Parsing mutations..."
    
    local mutation_count=0
    
    while IFS= read -r line; do
        [[ -z "$line" || "$line" =~ ^# ]] && continue
        
        mutation_count=$((mutation_count + 1))
        local mut_file="${DDG_DIR}/mutation_${mutation_count}.txt"
        
        # Count mutations in this line
        local num_muts=$(echo "$line" | tr ',' '\n' | wc -l)
        
        echo "$num_muts" > "$mut_file"
        
        # Parse each mutation
        echo "$line" | tr ',' '\n' | while read -r mut; do
            local wt_aa="${mut:0:1}"
            local pos="${mut:1:-1}"
            local mut_aa="${mut: -1}"
            echo "$wt_aa $pos $mut_aa" >> "$mut_file"
        done
        
        echo "$line" > "${DDG_DIR}/mutation_${mutation_count}_label.txt"
        
    done < "$MUTATIONS_FILE"
    
    success "Parsed $mutation_count mutation sets"
}

# Calculate ΔΔG for mutations
calculate_ddg() {
    log "Calculating ΔΔG values..."
    
    local ddg_results="${DDG_DIR}/ddg_results.txt"
    echo -e "Mutation\tddG_REU" > "$ddg_results"
    
    for mut_file in "${DDG_DIR}"/mutation_*.txt; do
        [[ ! -f "$mut_file" ]] && continue
        
        local mut_num=$(basename "$mut_file" | sed 's/mutation_\([0-9]*\).txt/\1/')
        local mut_label=$(cat "${DDG_DIR}/mutation_${mut_num}_label.txt")
        
        log "Processing mutation: $mut_label"
        
        cd "$DDG_DIR"
        
        "${ROSETTA_BIN}/ddg_monomer.linuxgccrelease" \
            -database "$ROSETTA_DB" \
            -in:file:s "../${RELAX_DIR}/relaxed_wildtype.pdb" \
            -ddg:mut_file "$mut_file" \
            -ddg:iterations "$DDG_ITERATIONS" \
            -ddg:dump_pdbs true \
            -ddg:local_opt_only false \
            -ddg:min_cst true \
            -ddg:mean false \
            -ddg:min true \
            -ddg:sc_min_only false \
            -ddg:ramp_repulsive true \
            -out:prefix "mut_${mut_num}_" \
            > "ddg_${mut_num}.log" 2>&1
        
        # Extract ΔΔG value
        local ddg_value=$(grep "ddG" "ddg_${mut_num}.log" | tail -1 | awk '{print $2}')
        echo -e "${mut_label}\t${ddg_value}" >> "$ddg_results"
        
        cd - > /dev/null
    done
    
    success "ΔΔG calculations completed"
}

# Create PyMOL visualization script
create_pymol_script() {
    log "Creating PyMOL visualization script..."
    
    cat > "${IMAGES_DIR}/visualize_mutations.py" << 'EOF'
#!/usr/bin/env python3

import pymol
from pymol import cmd
import os
import sys

def visualize_mutations(wt_pdb, mut_dir, output_dir):
    """Create superposition images for all mutations"""
    
    pymol.pymol_argv = ['pymol', '-qc']  # Quiet and no GUI
    pymol.finish_launching()
    
    # Load wild-type structure
    cmd.load(wt_pdb, 'wt')
    cmd.show('cartoon', 'wt')
    cmd.color('cyan', 'wt')
    
    # Find all mutant PDB files
    mut_files = []
    for f in os.listdir(mut_dir):
        if f.endswith('.pdb') and f.startswith('mut_'):
            mut_files.append(f)
    
    for mut_file in mut_files:
        mut_path = os.path.join(mut_dir, mut_file)
        mut_name = mut_file.replace('.pdb', '')
        
        # Load mutant structure
        cmd.load(mut_path, mut_name)
        
        # Align to wild-type
        cmd.align(mut_name, 'wt')
        
        # Style mutant structure
        cmd.show('cartoon', mut_name)
        cmd.color('orange', mut_name)
        
        # Set up view
        cmd.show('cartoon')
        cmd.bg_color('white')
        cmd.set('cartoon_transparency', 0.3, 'wt')
        
        # Save image
        output_file = os.path.join(output_dir, f'{mut_name}_superposition.png')
        cmd.png(output_file, width=800, height=600, dpi=300, ray=1)
        
        # Remove mutant for next iteration
        cmd.delete(mut_name)
        
        print(f"Created: {output_file}")
    
    cmd.quit()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python visualize_mutations.py <wt_pdb> <mut_dir> <output_dir>")
        sys.exit(1)
    
    wt_pdb = sys.argv[1]
    mut_dir = sys.argv[2] 
    output_dir = sys.argv[3]
    
    visualize_mutations(wt_pdb, mut_dir, output_dir)
EOF

    chmod +x "${IMAGES_DIR}/visualize_mutations.py"
    success "PyMOL visualization script created"
}

# Generate structure visualizations
generate_visualizations() {
    if ! command -v pymol &> /dev/null; then
        warning "PyMOL not available. Skipping visualization."
        return 0
    fi
    
    log "Generating structure visualizations..."
    
    python3 "${IMAGES_DIR}/visualize_mutations.py" \
        "${RELAX_DIR}/relaxed_wildtype.pdb" \
        "$DDG_DIR" \
        "$IMAGES_DIR"
    
    success "Structure visualizations generated"
}

# Create summary report
create_summary() {
    log "Creating summary report..."
    
    local summary_file="${OUTPUT_DIR}/pipeline_summary.txt"
    local protein_name=$(get_protein_name)
    
    cat > "$summary_file" << EOF
Rosetta Protein Folding and Mutation Analysis Pipeline Summary
============================================================

Protein: $protein_name
Analysis Date: $(date)

Input Files:
- Protein FASTA: $PROTEIN_FASTA
- Mutations File: $MUTATIONS_FILE

Pipeline Steps Completed:
1. Fragment generation
2. Ab initio folding ($NSTRUCT_FOLDING structures)
3. Structure relaxation ($NSTRUCT_RELAX structures)
4. ΔΔG calculations ($DDG_ITERATIONS iterations each)
5. Structure visualization

Output Files:
- Wild-type relaxed structure: ${RELAX_DIR}/relaxed_wildtype.pdb
- ΔΔG results table: ${DDG_DIR}/ddg_results.txt
- Structure images: ${IMAGES_DIR}/*.png
- Detailed logs: */log files in each subdirectory

Results Summary:
EOF

    # Add ΔΔG results if available
    if [[ -f "${DDG_DIR}/ddg_results.txt" ]]; then
        echo "" >> "$summary_file"
        echo "ΔΔG Results:" >> "$summary_file"
        echo "============" >> "$summary_file"
        cat "${DDG_DIR}/ddg_results.txt" >> "$summary_file"
    fi
    
    success "Summary report created: $summary_file"
}

# Main pipeline execution
main() {
    log "Starting Rosetta Protein Folding and Mutation Analysis Pipeline"
    
    check_dependencies
    setup_directories
    
    generate_fragments
    run_ab_initio_folding
    select_best_model
    relax_structure
    
    parse_mutations
    calculate_ddg
    
    create_pymol_script
    generate_visualizations
    
    create_summary
    
    success "Pipeline completed successfully!"
    log "Results available in: $OUTPUT_DIR"
}

# Error handling
trap 'error "Pipeline failed at line $LINENO"' ERR

# Run main pipeline
main "$@"