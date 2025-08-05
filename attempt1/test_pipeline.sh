#!/bin/bash

# Test script for Rosetta Protein Folding and Mutation Analysis Pipeline
# Author: Claude Code Assistant

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
    exit 1
}

success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Test functions
test_python_scripts() {
    log "Testing Python scripts syntax..."
    
    local scripts=(
        "parse_mutations.py"
        "aggregate_ddg_results.py"
        "visualize_structures.py"
    )
    
    for script in "${scripts[@]}"; do
        if [[ -f "$script" ]]; then
            if python3 -m py_compile "$script"; then
                success "✓ $script syntax OK"
            else
                error "✗ $script syntax error"
            fi
        else
            warning "Script not found: $script"
        fi
    done
}

test_mutation_parsing() {
    log "Testing mutation parsing..."
    
    if [[ -f "parse_mutations.py" && -f "example_mutations.txt" ]]; then
        # Test dry run
        if ./parse_mutations.py --dry-run example_mutations.txt test_output/; then
            success "✓ Mutation parsing dry run successful"
        else
            error "✗ Mutation parsing failed"
        fi
        
        # Test with FASTA validation
        if [[ -f "example_protein.faa" ]]; then
            if ./parse_mutations.py --validate-fasta example_protein.faa --dry-run example_mutations.txt test_output/; then
                success "✓ Mutation parsing with FASTA validation successful"
            else
                warning "Mutation parsing with FASTA validation failed"
            fi
        fi
    else
        warning "Missing files for mutation parsing test"
    fi
}

test_file_permissions() {
    log "Testing file permissions..."
    
    local executables=(
        "rosetta_mutation_pipeline.sh"
        "parse_mutations.py"
        "aggregate_ddg_results.py"
        "visualize_structures.py"
        "test_pipeline.sh"
    )
    
    for exe in "${executables[@]}"; do
        if [[ -f "$exe" ]]; then
            if [[ -x "$exe" ]]; then
                success "✓ $exe is executable"
            else
                warning "✗ $exe is not executable"
                chmod +x "$exe"
                success "✓ Made $exe executable"
            fi
        else
            warning "File not found: $exe"
        fi
    done
}

test_dependencies() {
    log "Testing dependencies..."
    
    # Test Python
    if command -v python3 &> /dev/null; then
        local python_version=$(python3 --version 2>&1 | cut -d' ' -f2)
        success "✓ Python 3 found: $python_version"
        
        # Test Python modules
        local modules=("pandas" "numpy" "pathlib" "argparse")
        for module in "${modules[@]}"; do
            if python3 -c "import $module" 2>/dev/null; then
                success "✓ Python module $module available"
            else
                warning "✗ Python module $module not found"
            fi
        done
    else
        error "Python 3 not found"
    fi
    
    # Test PyMOL (optional)
    if command -v pymol &> /dev/null; then
        success "✓ PyMOL command line tool found"
    else
        warning "PyMOL command line tool not found (optional)"
    fi
    
    if python3 -c "import pymol" 2>/dev/null; then
        success "✓ PyMOL Python module available"
    else
        warning "PyMOL Python module not found (optional)"
    fi
    
    # Test Snakemake (optional)
    if command -v snakemake &> /dev/null; then
        success "✓ Snakemake found"
    else
        warning "Snakemake not found (optional)"
    fi
}

test_input_files() {
    log "Testing input file formats..."
    
    # Test FASTA file
    if [[ -f "example_protein.faa" ]]; then
        if grep -q "^>" example_protein.faa && grep -q "[ACDEFGHIKLMNPQRSTVWY]" example_protein.faa; then
            success "✓ Example FASTA file format OK"
        else
            warning "Example FASTA file format may be invalid"
        fi
    else
        warning "Example FASTA file not found"
    fi
    
    # Test mutations file  
    if [[ -f "example_mutations.txt" ]]; then
        if grep -q "^[A-Z][0-9][A-Z]" example_mutations.txt; then
            success "✓ Example mutations file format OK"
        else
            warning "Example mutations file format may be invalid"
        fi
    else
        warning "Example mutations file not found"
    fi
}

test_rosetta_path() {
    log "Testing Rosetta installation..."
    
    local rosetta_path="/usr/local/rosetta.source.release-371"
    
    if [[ -d "$rosetta_path" ]]; then
        success "✓ Rosetta directory found: $rosetta_path"
        
        local rosetta_bin="${rosetta_path}/main/source/bin"
        if [[ -d "$rosetta_bin" ]]; then
            success "✓ Rosetta bin directory found"
            
            # Check for key executables
            local required_exes=(
                "fragment_picker"
                "AbinitioRelax"
                "relax"
                "ddg_monomer"
                "extract_pdbs"
            )
            
            for exe in "${required_exes[@]}"; do
                local found=false
                for ext in ".linuxgccrelease" ".macosclangrelease" ".default.linuxgccrelease"; do
                    if [[ -f "${rosetta_bin}/${exe}${ext}" ]]; then
                        success "✓ Found ${exe}${ext}"
                        found=true
                        break
                    fi
                done
                
                if [[ "$found" == false ]]; then
                    warning "✗ Rosetta executable not found: $exe"
                fi
            done
            
        else
            warning "Rosetta bin directory not found: $rosetta_bin"
        fi
        
        local rosetta_db="${rosetta_path}/main/database"
        if [[ -d "$rosetta_db" ]]; then
            success "✓ Rosetta database found"
        else
            warning "Rosetta database not found: $rosetta_db"
        fi
        
    else
        warning "Rosetta installation not found at: $rosetta_path"
        warning "Please update ROSETTA_PATH in the scripts"
    fi
}

create_test_files() {
    log "Creating test output directory..."
    mkdir -p test_output
    success "✓ Test output directory created"
}

run_comprehensive_test() {
    log "Running comprehensive pipeline test..."
    
    # Create small test files if they don't exist
    if [[ ! -f "protein.faa" ]]; then
        cp example_protein.faa protein.faa
        log "Using example protein as test input"
    fi
    
    if [[ ! -f "mutations.txt" ]]; then
        cp example_mutations.txt mutations.txt  
        log "Using example mutations as test input"
    fi
    
    # Test mutation parsing with actual execution
    if ./parse_mutations.py mutations.txt test_output/ddg_input/; then
        success "✓ Mutation parsing completed successfully"
        
        # Check output files
        if [[ -d "test_output/ddg_input" ]]; then
            local mut_files=$(find test_output/ddg_input -name "mutation_*.txt" | wc -l)
            local label_files=$(find test_output/ddg_input -name "*_label.txt" | wc -l)
            success "✓ Created $mut_files mutation files and $label_files label files"
        fi
    else
        warning "Mutation parsing test failed"
    fi
    
    # Test visualization script preparation
    if ./visualize_structures.py --script-only --list-mutants example_protein.faa test_output/ddg_input/ test_output/images/ 2>/dev/null; then
        success "✓ Visualization script preparation successful"
    else
        warning "Visualization script test had issues"
    fi
}

cleanup_test_files() {
    log "Cleaning up test files..."
    rm -rf test_output/
    success "✓ Test files cleaned up"
}

# Main test execution
main() {
    log "Starting Rosetta Pipeline Test Suite"
    
    create_test_files
    test_file_permissions
    test_python_scripts
    test_dependencies
    test_input_files
    test_mutation_parsing
    test_rosetta_path
    run_comprehensive_test
    
    success "Pipeline test suite completed!"
    log "Review any warnings above. The pipeline should be ready to use."
    
    if [[ "${1:-}" != "--keep-files" ]]; then
        cleanup_test_files
    fi
}

# Run tests
main "$@"