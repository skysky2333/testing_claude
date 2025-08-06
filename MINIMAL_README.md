# Rosetta Pipeline - Minimal Version

A streamlined, idempotent Python script for Rosetta protein folding and mutation analysis.

## Features

- **Simple & Direct**: Single Python file, ~500 lines
- **Idempotent**: Safe to re-run - skips completed steps automatically
- **No Dependencies**: Only uses Python standard library (no pandas, numpy, yaml)
- **File-based Checks**: Uses output file existence for step completion detection

## Quick Start

### 1. Setup Input Files

Create these files in your working directory:

**`protein.faa`** - Your protein sequence:
```
>my_protein
MKFLVLLFNILCLFPVLAADNHGVGPQGAS...
```

**`mutations.txt`** - Mutations to test:
```
A1G
S2M
A1G,S2M
```

### 2. Configure Rosetta Path

Edit the script's top section or use command line:
```python
ROSETTA_PATH = "/path/to/your/rosetta.source.release-371"
```

### 3. Run Pipeline

```bash
# Basic run
python rosetta_pipeline_minimal.py

# With custom Rosetta path
python rosetta_pipeline_minimal.py --rosetta-path /path/to/rosetta

# With custom files
python rosetta_pipeline_minimal.py --protein-fasta my_protein.faa --mutations-file my_mutations.txt
```

## Workflow Steps

The pipeline runs these steps in sequence:

1. **Validate** - Check Rosetta installation and input files
2. **Setup** - Create output directories
3. **Fragments** - Generate 3mer/9mer fragments (idempotent: checks if `*.3mers` and `*.9mers` exist)
4. **Folding** - Ab initio folding ~1000 structures (idempotent: checks if `best_model.pdb` exists)
5. **Relax** - Energy minimize structure (idempotent: checks if `relaxed_wildtype.pdb` exists)
6. **Parse** - Convert mutations to ddg_monomer format (always runs, fast)
7. **DDG** - Calculate ΔΔG for each mutation (idempotent: checks individual log/pdb files)
8. **Results** - Aggregate final table (idempotent: checks if `ddg_results.txt` exists)

## Idempotent Behavior

Each step checks for its expected output files before running:

- **Fragments**: Skips if both `.3mers` and `.9mers` files exist
- **Folding**: Skips if `best_model.pdb` exists
- **Relax**: Skips if `relaxed_wildtype.pdb` exists  
- **DDG**: Skips individual mutations if both `.log` and `_ddg_mut_1.pdb` exist
- **Results**: Skips if `ddg_results.txt` exists

This means you can safely interrupt and re-run the pipeline - it will continue from where it left off.

## Command Line Options

```bash
python rosetta_pipeline_minimal.py [options]

Options:
  --rosetta-path PATH      Path to Rosetta installation
  --protein-fasta FILE     Protein FASTA file (default: protein.faa)
  --mutations-file FILE    Mutations file (default: mutations.txt)  
  --output-dir DIR         Output directory (default: rosetta_results)
  --nstruct-folding N      Number of folding structures (default: 1000)
  --ddg-iterations N       DDG iterations (default: 50)
```

## Output Structure

```
rosetta_results/
├── fragments/
│   ├── protein.200.3mers
│   ├── protein.200.9mers
│   └── fragment_picker_config
├── folding/
│   ├── folding_results.out
│   └── best_model.pdb
├── relax/
│   └── relaxed_wildtype.pdb
└── ddg/
    ├── mutation_1.txt
    ├── mutation_1_label.txt
    ├── ddg_1.log
    ├── mut_1_ddg_mut_1.pdb
    ├── ...
    └── ddg_results.txt
```

## Key Output Files

- **`relaxed_wildtype.pdb`** - Final wild-type structure
- **`ddg_results.txt`** - Final results table:
  ```
  Mutation    ddG_REU
  A1G         -0.845
  S2M          1.234
  A1G,S2M      0.389
  ```

## Configuration

Edit these variables at the top of the script:

```python
# Default paths and files
ROSETTA_PATH = "/usr/local/rosetta.source.release-371"
PROTEIN_FASTA = "protein.faa"
MUTATIONS_FILE = "mutations.txt"
OUTPUT_DIR = "rosetta_results"

# Computational parameters
NSTRUCT_FOLDING = 1000    # More = better sampling, longer runtime
NSTRUCT_RELAX = 10        # Usually sufficient
DDG_ITERATIONS = 50       # More = better DDG accuracy
```

## Requirements

- **Rosetta 3.14+** (compiled from source)
- **Python 3.6+** (standard library only)
- **8+ GB RAM** (for folding)

## Typical Runtime

- **Fragments**: 5-30 minutes
- **Folding**: 1-12 hours (1000 structures)
- **Relax**: 10-60 minutes  
- **DDG**: 10-30 minutes per mutation

## Troubleshooting

### Common Issues

**"Rosetta bin directory not found"**
- Check `ROSETTA_PATH` points to correct installation
- Verify directory structure: `$ROSETTA_PATH/main/source/bin/`

**"Required executable not found"**
- Check Rosetta compiled successfully
- Look for alternative extensions (`.macosclangrelease`, `.default.linuxgccrelease`)

**"Fragment generation failed"**  
- Ensure internet connection (downloads database files)
- Check protein sequence is valid amino acids only

**"No valid scores found"**
- Folding may have failed completely
- Try reducing `NSTRUCT_FOLDING` or check input sequence

### Resuming Interrupted Runs

Simply re-run the same command - the script will automatically:
- Skip completed fragments
- Skip folding if model exists  
- Skip relaxation if structure exists
- Skip individual DDG calculations that completed
- Regenerate final results table if needed

### Log Files

Check these for detailed error information:
- Rosetta output goes to respective directories
- Script prints timestamp and status messages
- Individual DDG logs: `ddg/ddg_*.log`

## Performance Tips

- Use **faster storage** (SSD) for better I/O
- Adjust **NSTRUCT_FOLDING** based on available time
- Run on **cluster** for large mutation sets
- **Monitor memory** usage during folding step

## License

This script is provided under the same terms as Rosetta. Rosetta requires separate academic/commercial licensing.