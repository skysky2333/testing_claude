# Rosetta Protein Folding and Mutation Analysis Pipeline

A comprehensive pipeline for ab initio protein folding and ΔΔG mutation analysis using Rosetta.

## Overview

This pipeline performs:
1. **Ab initio folding** from protein sequence using fragment-based sampling
2. **Structure relaxation** to optimize the folded structure
3. **ΔΔG calculations** for mutations using ddg_monomer with cartesian sampling
4. **Structure visualization** and superposition using PyMOL
5. **Results aggregation** in tabular format

## Requirements

### Software Dependencies
- **Rosetta 3.14+** (compiled from source)
- **Python 3.7+** with pandas, numpy
- **PyMOL** (optional, for visualization)
- **Snakemake** (optional, for workflow execution)

### System Requirements
- Linux/macOS operating system
- 8+ GB RAM (more for larger proteins)
- Multiple CPU cores recommended

## Installation

1. **Install Rosetta**:
   ```bash
   # Download from https://downloads.rosettacommons.org/software/academic/
   # Compile according to Rosetta documentation
   # Set ROSETTA_PATH in scripts or config.yaml
   ```

2. **Install Python dependencies**:
   ```bash
   pip install pandas numpy snakemake
   ```

3. **Install PyMOL** (optional):
   ```bash
   # Via conda
   conda install -c conda-forge pymol-open-source
   
   # Or via pip
   pip install pymol-open-source
   ```

## Usage

### Input Files

1. **Protein sequence** (`protein.faa`):
   ```
   >protein_name Description
   MKFLVLLFNILCLFPVLAADNHGVGPQGAS...
   ```

2. **Mutations file** (`mutations.txt`):
   ```
   # Single mutations
   A1G
   S2M
   
   # Multiple simultaneous mutations
   A1G,S2M
   ```

### Execution Methods

#### Method 1: Bash Script (Recommended for single runs)

```bash
# Configure Rosetta path in script
vim rosetta_mutation_pipeline.sh

# Run pipeline
./rosetta_mutation_pipeline.sh
```

#### Method 2: Snakemake Workflow (Recommended for parallel execution)

```bash
# Configure settings
vim config.yaml

# Dry run to check workflow
snakemake -n

# Execute workflow
snakemake --cores 8
```

#### Method 3: Individual Scripts

```bash
# Parse mutations
./parse_mutations.py mutations.txt ddg_input/

# Aggregate results
./aggregate_ddg_results.py ddg_output/ -o results.tsv

# Create visualizations
./visualize_structures.py wildtype.pdb ddg_output/ images/
```

## Configuration

### Bash Script Configuration
Edit variables at the top of `rosetta_mutation_pipeline.sh`:
```bash
ROSETTA_PATH="/usr/local/rosetta.source.release-371"
NSTRUCT_FOLDING=1000
NSTRUCT_RELAX=10
DDG_ITERATIONS=50
```

### Snakemake Configuration
Edit `config.yaml`:
```yaml
rosetta_path: "/usr/local/rosetta.source.release-371"
protein_fasta: "protein.faa"
mutations_file: "mutations.txt"
nstruct_folding: 1000
nstruct_relax: 10
ddg_iterations: 50
```

## Output Files

### Directory Structure
```
rosetta_results/
├── fragments/              # Fragment files for ab initio folding
│   ├── protein.200.3mers
│   └── protein.200.9mers
├── folding/                # Ab initio folding results
│   ├── folding_results.out
│   └── best_model.pdb
├── relax/                  # Relaxed structures
│   └── relaxed_wildtype.pdb
├── ddg/                    # ΔΔG calculations
│   ├── mutation_*.txt
│   ├── ddg_*.log
│   ├── mut_*_ddg_mut_1.pdb
│   └── ddg_results.txt
├── images/                 # Structure visualizations
│   └── *_superposition.png
└── pipeline_summary.txt    # Summary report
```

### Key Output Files

1. **`relaxed_wildtype.pdb`**: Final wild-type structure
2. **`ddg_results.txt`**: ΔΔG values for all mutations
3. **`*_superposition.png`**: Structure superposition images
4. **`pipeline_summary.txt`**: Complete analysis summary

## Example Results

### ΔΔG Results Table
```
Mutation    ddG_REU
A1G         -0.845
S2M          1.234
A1G,S2M      0.389
```

- **Negative values**: Stabilizing mutations
- **Positive values**: Destabilizing mutations  
- **Units**: Rosetta Energy Units (REU)

## Advanced Usage

### Custom Fragment Generation
```bash
# Modify fragment picker settings in the pipeline
# Or use external fragment servers like ROBETTA
```

### Constraint-Based Relaxation
```bash
# Add constraints to relax protocol
-constraints:cst_fa_file constraints.cst
-relax:constrain_relax_to_start_coords
```

### High-Throughput Analysis
```bash
# Use Snakemake for parallel execution
snakemake --cluster "sbatch --cpus-per-task=4" --jobs 100
```

## Troubleshooting

### Common Issues

1. **Fragment generation fails**:
   - Check internet connection (for database access)
   - Verify Rosetta database path
   - Ensure protein sequence is valid

2. **Ab initio folding produces poor models**:
   - Increase `NSTRUCT_FOLDING` (try 5000-10000)
   - Add secondary structure predictions
   - Consider using comparative modeling instead

3. **ΔΔG calculations fail**:
   - Check that relaxed structure is valid
   - Verify mutation positions exist in sequence
   - Increase DDG_ITERATIONS for better sampling

4. **High memory usage**:
   - Reduce `NSTRUCT_FOLDING`
   - Run on cluster with more memory
   - Process mutations sequentially

### Performance Tips

- **Use SSDs** for better I/O performance
- **Parallel execution** with Snakemake
- **Cluster computing** for large-scale analysis
- **Fragment caching** to avoid regeneration

## Citation

If you use this pipeline, please cite:
- Rosetta software suite
- ddg_monomer protocol papers
- PyMOL (if used for visualization)

## Support

For issues with:
- **Rosetta**: https://forum.rosettacommons.org/
- **This pipeline**: Create an issue in the repository
- **PyMOL**: https://pymol.org/support

## License

This pipeline is provided under the same license terms as Rosetta. See Rosetta documentation for details.