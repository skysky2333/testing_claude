# Snakemake workflow for Rosetta Protein Folding and Mutation Analysis
# Author: Claude Code Assistant

import os
import pandas as pd
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Rosetta paths
ROSETTA_PATH = config.get("rosetta_path", "/usr/local/rosetta.source.release-371")
ROSETTA_BIN = f"{ROSETTA_PATH}/main/source/bin"
ROSETTA_DB = f"{ROSETTA_PATH}/main/database"

# Input files
PROTEIN_FASTA = config.get("protein_fasta", "protein.faa")
MUTATIONS_FILE = config.get("mutations_file", "mutations.txt")

# Parameters
NSTRUCT_FOLDING = config.get("nstruct_folding", 1000)
NSTRUCT_RELAX = config.get("nstruct_relax", 10)
DDG_ITERATIONS = config.get("ddg_iterations", 50)

# Output directories
OUTPUT_DIR = "rosetta_results"
FRAGMENTS_DIR = f"{OUTPUT_DIR}/fragments"
FOLDING_DIR = f"{OUTPUT_DIR}/folding"
RELAX_DIR = f"{OUTPUT_DIR}/relax"
DDG_DIR = f"{OUTPUT_DIR}/ddg"
IMAGES_DIR = f"{OUTPUT_DIR}/images"

# Parse mutations to determine number of mutation sets
def parse_mutations():
    """Parse mutations file and return list of mutation IDs"""
    mutation_ids = []
    if os.path.exists(MUTATIONS_FILE):
        with open(MUTATIONS_FILE, 'r') as f:
            for i, line in enumerate(f, 1):
                line = line.strip()
                if line and not line.startswith('#'):
                    mutation_ids.append(str(i))
    return mutation_ids

MUTATION_IDS = parse_mutations()

# Helper function to get protein name
def get_protein_name():
    """Extract protein name from FASTA file"""
    if os.path.exists(PROTEIN_FASTA):
        with open(PROTEIN_FASTA, 'r') as f:
            header = f.readline().strip()
            if header.startswith('>'):
                name = header[1:].split()[0].replace('|', '_')
                return name if name else "protein"
    return "protein"

PROTEIN_NAME = get_protein_name()

# Target rule - defines final outputs
rule all:
    input:
        f"{RELAX_DIR}/relaxed_wildtype.pdb",
        f"{DDG_DIR}/ddg_results.txt",
        f"{OUTPUT_DIR}/pipeline_summary.txt"

# Generate fragment files
rule generate_fragments:
    input:
        fasta = PROTEIN_FASTA
    output:
        frag3 = f"{FRAGMENTS_DIR}/{PROTEIN_NAME}.200.3mers",
        frag9 = f"{FRAGMENTS_DIR}/{PROTEIN_NAME}.200.9mers",
        config = f"{FRAGMENTS_DIR}/fragment_picker_config"
    params:
        prefix = f"{FRAGMENTS_DIR}/{PROTEIN_NAME}",
        rosetta_bin = ROSETTA_BIN,
        rosetta_db = ROSETTA_DB
    log:
        f"{FRAGMENTS_DIR}/fragment_picker.log"
    shell:
        """
        mkdir -p {FRAGMENTS_DIR}
        
        cat > {output.config} << EOF
database {params.rosetta_db}
in:file:fasta {input.fasta}
frags:describe_fragments {FRAGMENTS_DIR}/{PROTEIN_NAME}_frags.fsc
out:file:frag_prefix {params.prefix}
frags:n_candidates 200
frags:n_frags 25
frags:frag_sizes 3 9
EOF

        cd {FRAGMENTS_DIR}
        {params.rosetta_bin}/fragment_picker.linuxgccrelease \
            @fragment_picker_config \
            > {log} 2>&1
        """

# Ab initio folding
rule ab_initio_folding:
    input:
        fasta = PROTEIN_FASTA,
        frag3 = f"{FRAGMENTS_DIR}/{PROTEIN_NAME}.200.3mers",
        frag9 = f"{FRAGMENTS_DIR}/{PROTEIN_NAME}.200.9mers"
    output:
        silent = f"{FOLDING_DIR}/folding_results.out"
    params:
        rosetta_bin = ROSETTA_BIN,
        rosetta_db = ROSETTA_DB,
        nstruct = NSTRUCT_FOLDING
    log:
        f"{FOLDING_DIR}/folding.log"
    threads: 4
    shell:
        """
        mkdir -p {FOLDING_DIR}
        cd {FOLDING_DIR}
        
        {params.rosetta_bin}/AbinitioRelax.linuxgccrelease \
            -database {params.rosetta_db} \
            -in:file:fasta ../{input.fasta} \
            -in:file:frag3 ../{input.frag3} \
            -in:file:frag9 ../{input.frag9} \
            -abinitio:relax \
            -relax:fast \
            -nstruct {params.nstruct} \
            -out:file:silent folding_results.out \
            -out:nooutput \
            > {log} 2>&1
        """

# Select best folding model
rule select_best_model:
    input:
        silent = f"{FOLDING_DIR}/folding_results.out"
    output:
        best_pdb = f"{FOLDING_DIR}/best_model.pdb",
        best_score = f"{FOLDING_DIR}/best_model_score.txt"
    params:
        rosetta_bin = ROSETTA_BIN
    shell:
        """
        cd {FOLDING_DIR}
        
        # Extract scores and find best model
        grep "^SCORE:" folding_results.out | sort -k2,2n | head -1 > best_model_score.txt
        best_tag=$(awk '{{print $NF}}' best_model_score.txt)
        
        # Extract best model as PDB
        echo "$best_tag" > best_tag.txt
        {params.rosetta_bin}/extract_pdbs.linuxgccrelease \
            -in:file:silent folding_results.out \
            -in:file:tagfile best_tag.txt \
            -out:prefix "best_"
        
        mv "best_${{best_tag}}.pdb" best_model.pdb
        """

# Relax structure
rule relax_structure:
    input:
        pdb = f"{FOLDING_DIR}/best_model.pdb"
    output:
        relaxed = f"{RELAX_DIR}/relaxed_wildtype.pdb"
    params:
        rosetta_bin = ROSETTA_BIN,
        rosetta_db = ROSETTA_DB,
        nstruct = NSTRUCT_RELAX
    log:
        f"{RELAX_DIR}/relax.log"
    shell:
        """
        mkdir -p {RELAX_DIR}
        cd {RELAX_DIR}
        
        {params.rosetta_bin}/relax.linuxgccrelease \
            -database {params.rosetta_db} \
            -in:file:s ../{input.pdb} \
            -relax:fast \
            -relax:constrain_relax_to_start_coords \
            -relax:ramp_constraints false \
            -ex1 \
            -ex2 \
            -use_input_sc \
            -flip_HNQ \
            -no_optH false \
            -nstruct {params.nstruct} \
            -out:prefix "relaxed_" \
            -out:suffix "_wt" \
            > {log} 2>&1
        
        # Select best relaxed structure
        best_relaxed=$(ls relaxed_*_wt_*.pdb | head -1)
        cp "$best_relaxed" relaxed_wildtype.pdb
        """

# Parse mutations into individual files
rule parse_mutations:
    input:
        mutations = MUTATIONS_FILE
    output:
        mutation_files = expand(f"{DDG_DIR}/mutation_{{mut_id}}.txt", mut_id=MUTATION_IDS),
        label_files = expand(f"{DDG_DIR}/mutation_{{mut_id}}_label.txt", mut_id=MUTATION_IDS)
    run:
        import os
        os.makedirs(DDG_DIR, exist_ok=True)
        
        mutation_count = 0
        with open(input.mutations, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    mutation_count += 1
                    
                    # Create mutation file
                    mut_file = f"{DDG_DIR}/mutation_{mutation_count}.txt"
                    label_file = f"{DDG_DIR}/mutation_{mutation_count}_label.txt"
                    
                    # Count mutations in this line
                    mutations = line.split(',')
                    num_muts = len(mutations)
                    
                    with open(mut_file, 'w') as mf:
                        mf.write(f"{num_muts}\n")
                        
                        for mut in mutations:
                            mut = mut.strip()
                            wt_aa = mut[0]
                            pos = mut[1:-1]
                            mut_aa = mut[-1]
                            mf.write(f"{wt_aa} {pos} {mut_aa}\n")
                    
                    # Save label
                    with open(label_file, 'w') as lf:
                        lf.write(line)

# Calculate ΔΔG for individual mutations
rule calculate_ddg:
    input:
        pdb = f"{RELAX_DIR}/relaxed_wildtype.pdb",
        mut_file = f"{DDG_DIR}/mutation_{{mut_id}}.txt"
    output:
        log = f"{DDG_DIR}/ddg_{{mut_id}}.log",
        pdb = f"{DDG_DIR}/mut_{{mut_id}}_ddg_mut_1.pdb"
    params:
        rosetta_bin = ROSETTA_BIN,
        rosetta_db = ROSETTA_DB,
        iterations = DDG_ITERATIONS,
        prefix = lambda wildcards: f"mut_{wildcards.mut_id}_"
    shell:
        """
        cd {DDG_DIR}
        
        {params.rosetta_bin}/ddg_monomer.linuxgccrelease \
            -database {params.rosetta_db} \
            -in:file:s ../{input.pdb} \
            -ddg:mut_file mutation_{wildcards.mut_id}.txt \
            -ddg:iterations {params.iterations} \
            -ddg:dump_pdbs true \
            -ddg:local_opt_only false \
            -ddg:min_cst true \
            -ddg:mean false \
            -ddg:min true \
            -ddg:sc_min_only false \
            -ddg:ramp_repulsive true \
            -out:prefix {params.prefix} \
            > ddg_{wildcards.mut_id}.log 2>&1
        """

# Aggregate ΔΔG results
rule aggregate_ddg_results:
    input:
        logs = expand(f"{DDG_DIR}/ddg_{{mut_id}}.log", mut_id=MUTATION_IDS),
        labels = expand(f"{DDG_DIR}/mutation_{{mut_id}}_label.txt", mut_id=MUTATION_IDS)
    output:
        results = f"{DDG_DIR}/ddg_results.txt"
    run:
        with open(output.results, 'w') as out_f:
            out_f.write("Mutation\tddG_REU\n")
            
            for mut_id in MUTATION_IDS:
                log_file = f"{DDG_DIR}/ddg_{mut_id}.log"
                label_file = f"{DDG_DIR}/mutation_{mut_id}_label.txt"
                
                # Read mutation label
                with open(label_file, 'r') as lf:
                    mut_label = lf.read().strip()
                
                # Extract ΔΔG value from log
                ddg_value = "N/A"
                if os.path.exists(log_file):
                    with open(log_file, 'r') as lf:
                        for line in lf:
                            if 'ddG' in line:
                                parts = line.strip().split()
                                if len(parts) >= 2:
                                    try:
                                        ddg_value = float(parts[1])
                                        break
                                    except (ValueError, IndexError):
                                        continue
                
                out_f.write(f"{mut_label}\t{ddg_value}\n")

# Generate PyMOL visualization script
rule create_pymol_script:
    output:
        script = f"{IMAGES_DIR}/visualize_mutations.py"
    shell:
        """
        mkdir -p {IMAGES_DIR}
        
        cat > {output.script} << 'EOF'
#!/usr/bin/env python3

import pymol
from pymol import cmd
import os
import sys

def visualize_mutations(wt_pdb, mut_dir, output_dir):
    \"\"\"Create superposition images for all mutations\"\"\"
    
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
        
        chmod +x {output.script}
        """

# Generate structure visualizations
rule generate_visualizations:
    input:
        script = f"{IMAGES_DIR}/visualize_mutations.py",
        wt_pdb = f"{RELAX_DIR}/relaxed_wildtype.pdb",
        mut_pdbs = expand(f"{DDG_DIR}/mut_{{mut_id}}_ddg_mut_1.pdb", mut_id=MUTATION_IDS) if MUTATION_IDS else []
    output:
        images = expand(f"{IMAGES_DIR}/mut_{{mut_id}}_superposition.png", mut_id=MUTATION_IDS) if MUTATION_IDS else [f"{IMAGES_DIR}/no_mutations.txt"]
    shell:
        """
        if command -v pymol &> /dev/null; then
            python3 {input.script} {input.wt_pdb} {DDG_DIR} {IMAGES_DIR}
        else
            echo "PyMOL not found. Skipping visualization."
            # Create dummy files to satisfy rule
            for mut_id in {MUTATION_IDS}; do
                touch {IMAGES_DIR}/mut_${{mut_id}}_superposition.png
            done
        fi
        """

# Create summary report
rule create_summary:
    input:
        relaxed = f"{RELAX_DIR}/relaxed_wildtype.pdb",
        ddg_results = f"{DDG_DIR}/ddg_results.txt"
    output:
        summary = f"{OUTPUT_DIR}/pipeline_summary.txt"
    run:
        with open(output.summary, 'w') as f:
            f.write("Rosetta Protein Folding and Mutation Analysis Pipeline Summary\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Protein: {PROTEIN_NAME}\n")
            f.write(f"Analysis Date: {os.popen('date').read().strip()}\n\n")
            
            f.write("Input Files:\n")
            f.write(f"- Protein FASTA: {PROTEIN_FASTA}\n")
            f.write(f"- Mutations File: {MUTATIONS_FILE}\n\n")
            
            f.write("Pipeline Steps Completed:\n")
            f.write("1. Fragment generation\n")
            f.write(f"2. Ab initio folding ({NSTRUCT_FOLDING} structures)\n")
            f.write(f"3. Structure relaxation ({NSTRUCT_RELAX} structures)\n")
            f.write(f"4. ΔΔG calculations ({DDG_ITERATIONS} iterations each)\n")
            f.write("5. Structure visualization\n\n")
            
            f.write("Output Files:\n")
            f.write(f"- Wild-type relaxed structure: {input.relaxed}\n")
            f.write(f"- ΔΔG results table: {input.ddg_results}\n")
            f.write(f"- Structure images: {IMAGES_DIR}/*.png\n")
            f.write("- Detailed logs: */log files in each subdirectory\n\n")
            
            f.write("Results Summary:\n")
            f.write("=" * 15 + "\n")
            
            # Add ΔΔG results
            if os.path.exists(input.ddg_results):
                with open(input.ddg_results, 'r') as ddg_f:
                    f.write(ddg_f.read())