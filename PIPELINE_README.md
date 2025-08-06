# Rosetta Pipeline - Self-Contained Python Script

A comprehensive, idempotent Python script for protein folding and mutation analysis using the Rosetta software suite.

## Overview

`rosetta_pipeline.py` is a self-contained script that consolidates the entire Rosetta protein analysis workflow into a single executable with resume capability, robust error handling, and flexible configuration.

## Features

- **Self-contained**: Single Python file with all functionality embedded
- **Idempotent**: Safe to re-run - automatically skips completed steps
- **Resumable**: Continue interrupted runs with `--resume` flag
- **Configurable**: Multiple configuration methods (defaults → file → CLI)
- **Comprehensive**: Full workflow from fragments to visualization
- **Production-ready**: Robust error handling and logging

## Workflow

The pipeline executes these steps sequentially:

1. **Dependency Validation** - Check Rosetta installation and input files
2. **Directory Setup** - Create organized output structure
3. **Fragment Generation** - Generate 3mer/9mer fragments for folding
4. **Ab Initio Folding** - Generate ~1000 protein conformations
5. **Best Model Selection** - Select lowest energy structure
6. **Structure Relaxation** - Refine structure with energy minimization
7. **Mutation Parsing** - Parse mutations into ddg_monomer format
8. **ΔΔG Calculation** - Calculate stability changes for each mutation
9. **Result Aggregation** - Collect and analyze all ΔΔG values
10. **Visualization** - Create PyMOL structure superpositions
11. **Summary Report** - Generate comprehensive analysis report

## Requirements

### Required Dependencies
- **Rosetta 3.14+** (compiled from source)
- **Python 3.7+** with packages:
  - `pandas` - Data manipulation and analysis
  - `numpy` - Numerical computations
  - `PyYAML` - YAML configuration support

### Optional Dependencies
- **PyMOL** - For structure visualization (auto-skipped if unavailable)

## Installation

1. **Install Rosetta** from https://downloads.rosettacommons.org/software/academic/
2. **Install Python dependencies**:
   ```bash
   pip install pandas numpy pyyaml
   ```
3. **Optional: Install PyMOL**:
   ```bash
   pip install pymol-open-source
   ```

## Quick Start

### Basic Usage

```bash
# Run complete pipeline with defaults
python rosetta_pipeline.py

# Validate setup without running
python rosetta_pipeline.py --validate-only
```

### Input Files

Create these input files in your working directory:

**`protein.faa`** - Protein sequence in FASTA format:
```
>my_protein Description
MKFLVLLFNILCLFPVLAADNHGVGPQGAS...
```

**`mutations.txt`** - One mutation set per line:
```
A1G
S2M
A1G,S2M
```

### Resume Interrupted Runs

```bash
# Resume from last completed step
python rosetta_pipeline.py --resume
```

## Configuration

### Configuration Priority
1. **Built-in defaults** (lowest priority)
2. **Configuration file** (medium priority)
3. **Command line arguments** (highest priority)

### Configuration File Example

**`my_config.yaml`:**
```yaml
rosetta_path: "/usr/local/rosetta"
protein_fasta: "my_protein.faa"
mutations_file: "my_mutations.txt"
output_dir: "results"

# Computational parameters
nstruct_folding: 5000
nstruct_relax: 20
ddg_iterations: 100

# Advanced options
relaxation:
  fast_mode: false
  constrain_coords: true

visualization:
  image_width: 1200
  image_height: 900
  ray_tracing: true
```

### Command Line Examples

```bash
# Use custom config file
python rosetta_pipeline.py --config my_config.yaml

# Override specific parameters
python rosetta_pipeline.py --nstruct-folding 5000 --ddg-iterations 100

# Custom input/output
python rosetta_pipeline.py --protein-fasta my_protein.faa --output-dir results

# Skip visualization
python rosetta_pipeline.py --skip-visualization

# Verbose logging
python rosetta_pipeline.py --verbose
```

## Output Structure

```
rosetta_results/
├── fragments/              # Fragment files for ab initio folding
│   ├── protein.200.3mers
│   ├── protein.200.9mers
│   └── fragment_picker.log
├── folding/                # Ab initio folding results
│   ├── folding_results.out
│   ├── best_model.pdb
│   └── folding.log
├── relax/                  # Relaxed structures
│   ├── relaxed_wildtype.pdb
│   └── relax.log
├── ddg/                    # ΔΔG calculations
│   ├── mutation_*.txt      # ddg_monomer input files
│   ├── mutation_*_label.txt
│   ├── ddg_*.log          # Calculation logs
│   ├── mut_*_ddg_mut_1.pdb # Mutant structures
│   └── ddg_results.txt    # Final results table
├── images/                 # Structure visualizations
│   └── *_superposition.png
├── pipeline_state.json     # Resume state tracking
└── pipeline_summary.txt    # Comprehensive report
```

## Key Output Files

- **`relaxed_wildtype.pdb`** - Final wild-type structure
- **`ddg_results.txt`** - ΔΔG values for all mutations (TSV format)
- **`*_superposition.png`** - Structure alignment images
- **`pipeline_summary.txt`** - Complete analysis report with configuration and results

## Advanced Usage

### Customizing Rosetta Parameters

The pipeline uses sensible defaults but allows customization via configuration:

```yaml
fragment_picking:
  n_candidates: 200
  n_frags: 25
  frag_sizes: [3, 9]

relaxation:
  fast_mode: true          # vs thorough mode
  constrain_coords: true   # constrain to starting coordinates
  ramp_constraints: false

ddg_calculation:
  local_opt_only: false    # optimize beyond local shell
  min_cst: true           # use minimization constraints
  dump_pdbs: true         # save mutant PDB files
  ramp_repulsive: true    # ramp repulsive weights
```

### High-Throughput Analysis

For large mutation sets:

```bash
# Use more folding structures for better sampling
python rosetta_pipeline.py --nstruct-folding 10000

# Increase DDG sampling
python rosetta_pipeline.py --ddg-iterations 200

# Use custom output directory
python rosetta_pipeline.py --output-dir large_analysis
```

### Development and Debugging

```bash
# Verbose logging
python rosetta_pipeline.py --verbose

# Validate without running
python rosetta_pipeline.py --validate-only

# Quick test run
python rosetta_pipeline.py --nstruct-folding 100 --ddg-iterations 5
```

## State Management

The pipeline automatically tracks completed steps in `pipeline_state.json`:

```json
{
  "completed_steps": [
    "setup_directories",
    "generate_fragments",
    "ab_initio_folding"
  ],
  "last_run": "2024-01-15T10:30:45",
  "protein_name": "my_protein"
}
```

This enables:
- **Resume capability** - Pick up where you left off
- **Idempotent execution** - Safe to re-run completed steps
- **Progress tracking** - Monitor long-running analyses

## Error Handling

The pipeline includes comprehensive error handling:

- **Dependency validation** - Check Rosetta installation and input files
- **Step-wise validation** - Verify outputs before proceeding
- **Graceful failures** - Continue pipeline even if some mutations fail
- **Detailed logging** - Clear error messages and log file references
- **State preservation** - Resume capability even after failures

## Performance Considerations

### Memory Requirements
- **Fragment generation**: ~1-2 GB
- **Ab initio folding**: ~2-4 GB (depends on protein size and nstruct)
- **DDG calculations**: ~500 MB per mutation

### Runtime Estimates
- **Fragment generation**: 5-30 minutes
- **Ab initio folding**: 1-12 hours (1000 structures)
- **Relaxation**: 10-60 minutes
- **DDG calculations**: 10-30 minutes per mutation

### Optimization Tips
- Use **SSDs** for better I/O performance
- Adjust **nstruct_folding** based on protein size and available time
- Consider **cluster computing** for large mutation sets
- Use **--skip-visualization** to save time if images not needed

## Troubleshooting

### Common Issues

**Rosetta executables not found**:
- Verify `rosetta_path` points to correct installation
- Check that executables have proper extensions for your platform

**Fragment generation fails**:
- Ensure internet connection (for database access)
- Verify protein sequence is valid
- Check Rosetta database path

**High memory usage**:
- Reduce `nstruct_folding`
- Process mutations in smaller batches
- Use machines with more RAM

**PyMOL visualization fails**:
- Install PyMOL: `pip install pymol-open-source`
- Or skip visualization: `--skip-visualization`

### Log Files

Each step creates detailed log files:
- `fragments/fragment_picker.log`
- `folding/folding.log`
- `relax/relax.log`
- `ddg/ddg_*.log`

Check these for detailed error messages and Rosetta output.

## Citation

If you use this pipeline, please cite:
- **Rosetta software suite**: Alford et al. (2017) J Chem Theory Comput
- **ddg_monomer protocol**: Kellogg et al. (2011) Proteins
- **PyMOL** (if used): The PyMOL Molecular Graphics System

## License

This pipeline script is provided under the same license terms as Rosetta. The script itself is open source, but Rosetta requires academic/commercial licensing.

## Support

- **Rosetta issues**: https://forum.rosettacommons.org/
- **Script issues**: Check log files and error messages
- **PyMOL issues**: https://pymol.org/support