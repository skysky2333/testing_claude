# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a computational biology pipeline for protein structure prediction and mutation analysis using Rosetta. The pipeline performs ab initio protein folding from sequence, structure relaxation, ΔΔG calculations for mutations, and structure visualization.

## Core Workflow Architecture

The pipeline follows this sequential workflow:

1. **Fragment Generation**: Creates 3mer and 9mer fragments using `fragment_picker`
2. **Ab Initio Folding**: Generates ~1000 folded structures using `AbinitioRelax`
3. **Model Selection**: Selects best model by lowest energy score
4. **Structure Relaxation**: Refines structure using `relax` protocol
5. **Mutation Processing**: Parses mutations and creates ddg_monomer input files
6. **ΔΔG Calculation**: Calculates stability changes using `ddg_monomer`
7. **Visualization**: Creates structure superpositions using PyMOL

## Common Commands

### Run Full Pipeline
```bash
# Bash script method (single execution)
./rosetta_mutation_pipeline.sh

# Snakemake method (parallel execution)
snakemake --cores 8
```

### Test Pipeline Components  
```bash
# Run comprehensive tests
./test_pipeline.sh

# Test with file cleanup
./test_pipeline.sh --keep-files
```

### Individual Script Usage
```bash
# Parse mutations into ddg_monomer format
./parse_mutations.py mutations.txt ddg_input/
./parse_mutations.py --validate-fasta protein.faa mutations.txt ddg_input/

# Aggregate ΔΔG results
./aggregate_ddg_results.py ddg_output/ -o results.tsv

# Create structure visualizations  
./visualize_structures.py wildtype.pdb ddg_output/ images/
```

### Development Testing
```bash
# Check Python syntax
python3 -m py_compile *.py

# Validate mutation parsing
./parse_mutations.py --dry-run example_mutations.txt test_output/

# List available mutant structures
./visualize_structures.py --list-mutants wildtype.pdb ddg_output/ images/
```

## Configuration

### Primary Configuration Files
- **`config.yaml`**: Snakemake workflow parameters and Rosetta settings
- **`rosetta_mutation_pipeline.sh`**: Direct script configuration (lines 10-29)

### Key Configuration Parameters
- `rosetta_path`: Path to Rosetta installation (default: `/usr/local/rosetta.source.release-371`)
- `nstruct_folding`: Number of ab initio structures (default: 1000)
- `nstruct_relax`: Number of relaxation structures (default: 10) 
- `ddg_iterations`: DDG sampling iterations (default: 50)

## Input File Formats

### Protein Sequence (`protein.faa`)
Standard FASTA format:
```
>protein_name Description
MKFLVLLFNILCLFPVLAADNHGVGPQGAS...
```

### Mutations (`mutations.txt`)
One mutation set per line, comma-separated for multiple simultaneous mutations:
```
A1G
S2M
A1G,S2M
```

## Project Structure

```
rosetta_results/          # Main output directory
├── fragments/           # Fragment files for folding
├── folding/            # Ab initio folding results  
├── relax/              # Relaxed structures
├── ddg/                # ΔΔG calculations and mutant PDBs
└── images/             # Structure visualization images
```

## Dependencies

### Required
- **Rosetta 3.14+** (compiled from source)
- **Python 3.7+** with pandas, numpy
- **Snakemake** (for parallel workflow execution)

### Optional  
- **PyMOL** (for structure visualization)

## Common Development Patterns

### Error Handling
- All Python scripts use proper argument parsing with `argparse`
- Bash scripts use `set -euo pipefail` for strict error handling
- Log files are generated for each major step

### File Validation
- FASTA sequences validated against standard amino acid codes
- Mutation positions checked against sequence length
- PDB file existence verified before processing

### Parallel Execution
- Snakemake handles parallel DDG calculations automatically
- Individual mutation files processed independently
- Thread counts configurable in `config.yaml`

## Testing Strategy

The test suite (`test_pipeline.sh`) validates:
- Python script syntax and executability
- Input file formats and dependencies
- Rosetta installation and executables
- Mutation parsing with dry-run validation
- File permissions and directory structure