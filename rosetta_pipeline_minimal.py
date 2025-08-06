#!/usr/bin/env python3
"""
Minimal idempotent Rosetta protein folding and mutation analysis pipeline.
Simple, direct approach with file-based idempotency checks.

Author: Claude Code Assistant
"""

import os
import sys
import re
import argparse
import subprocess
from datetime import datetime

# Default configuration - edit these as needed
ROSETTA_PATH = "/usr/local/rosetta.source.release-371"
PROTEIN_FASTA = "protein.faa"
MUTATIONS_FILE = "mutations.txt"
OUTPUT_DIR = "rosetta_results"

# Pipeline parameters
NSTRUCT_FOLDING = 1000
NSTRUCT_RELAX = 10
DDG_ITERATIONS = 50

def log(msg):
    """Simple logging with timestamp."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")

def run_command(cmd, cwd=None):
    """Run a command and check for errors."""
    log(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, cwd=cwd, check=True, 
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return True
    except subprocess.CalledProcessError as e:
        log(f"ERROR: Command failed with return code {e.returncode}")
        log(f"STDERR: {e.stderr}")
        return False

def validate_setup():
    """Validate Rosetta installation and input files."""
    log("Validating setup...")
    
    # Check Rosetta
    rosetta_bin = os.path.join(ROSETTA_PATH, "main", "source", "bin")
    if not os.path.exists(rosetta_bin):
        log(f"ERROR: Rosetta bin directory not found: {rosetta_bin}")
        return False
    
    # Check required executables
    required_exes = [
        "fragment_picker.linuxgccrelease",
        "AbinitioRelax.linuxgccrelease", 
        "relax.linuxgccrelease",
        "ddg_monomer.linuxgccrelease",
        "extract_pdbs.linuxgccrelease"
    ]
    
    for exe in required_exes:
        exe_path = os.path.join(rosetta_bin, exe)
        if not os.path.exists(exe_path):
            # Try alternative extensions
            alt_found = False
            for ext in ['.macosclangrelease', '.default.linuxgccrelease']:
                alt_exe = exe.replace('.linuxgccrelease', ext)
                if os.path.exists(os.path.join(rosetta_bin, alt_exe)):
                    alt_found = True
                    break
            if not alt_found:
                log(f"ERROR: Required executable not found: {exe}")
                return False
    
    # Check input files
    if not os.path.exists(PROTEIN_FASTA):
        log(f"ERROR: Protein FASTA not found: {PROTEIN_FASTA}")
        return False
    
    if not os.path.exists(MUTATIONS_FILE):
        log(f"ERROR: Mutations file not found: {MUTATIONS_FILE}")
        return False
    
    log("Validation passed")
    return True

def setup_directories():
    """Create output directory structure."""
    log("Setting up directories...")
    
    dirs = [OUTPUT_DIR, f"{OUTPUT_DIR}/fragments", f"{OUTPUT_DIR}/folding", 
            f"{OUTPUT_DIR}/relax", f"{OUTPUT_DIR}/ddg"]
    
    for d in dirs:
        os.makedirs(d, exist_ok=True)
    
    log("Directories created")

def get_protein_name():
    """Extract protein name from FASTA file."""
    try:
        with open(PROTEIN_FASTA, 'r') as f:
            header = f.readline().strip()
            if header.startswith('>'):
                name = header[1:].split()[0].replace('|', '_').replace('/', '_')
                return name if name else "protein"
    except:
        pass
    return "protein"

def generate_fragments():
    """Generate fragment files for ab initio folding."""
    log("Checking fragments...")
    
    protein_name = get_protein_name()
    frag3_file = f"{OUTPUT_DIR}/fragments/{protein_name}.200.3mers"
    frag9_file = f"{OUTPUT_DIR}/fragments/{protein_name}.200.9mers"
    
    # Check if fragments already exist (idempotent)
    if os.path.exists(frag3_file) and os.path.exists(frag9_file):
        log("Fragment files already exist, skipping generation")
        return True
    
    log("Generating fragments...")
    
    # Create fragment picker config
    config_file = f"{OUTPUT_DIR}/fragments/fragment_picker_config"
    rosetta_db = os.path.join(ROSETTA_PATH, "main", "database")
    
    config_content = f"""database {rosetta_db}
in:file:fasta {PROTEIN_FASTA}
frags:describe_fragments {OUTPUT_DIR}/fragments/{protein_name}_frags.fsc
out:file:frag_prefix {OUTPUT_DIR}/fragments/{protein_name}
frags:n_candidates 200
frags:n_frags 25
frags:frag_sizes 3 9
"""
    
    with open(config_file, 'w') as f:
        f.write(config_content)
    
    # Run fragment picker
    rosetta_bin = os.path.join(ROSETTA_PATH, "main", "source", "bin")
    cmd = [os.path.join(rosetta_bin, "fragment_picker.linuxgccrelease"), f"@{config_file}"]
    
    if not run_command(cmd, cwd=f"{OUTPUT_DIR}/fragments"):
        log("ERROR: Fragment generation failed")
        return False
    
    # Verify outputs
    if os.path.exists(frag3_file) and os.path.exists(frag9_file):
        log("Fragment generation completed")
        return True
    else:
        log("ERROR: Fragment files not created")
        return False

def run_folding():
    """Perform ab initio protein folding."""
    log("Checking folding...")
    
    silent_file = f"{OUTPUT_DIR}/folding/folding_results.out"
    best_pdb = f"{OUTPUT_DIR}/folding/best_model.pdb"
    
    # Check if folding already done (idempotent)
    if os.path.exists(best_pdb):
        log("Best folding model already exists, skipping")
        return True
    
    # Check if we need to run folding
    if not os.path.exists(silent_file):
        log("Running ab initio folding...")
        
        protein_name = get_protein_name()
        frag3_file = f"{OUTPUT_DIR}/fragments/{protein_name}.200.3mers"
        frag9_file = f"{OUTPUT_DIR}/fragments/{protein_name}.200.9mers"
        rosetta_db = os.path.join(ROSETTA_PATH, "main", "database")
        rosetta_bin = os.path.join(ROSETTA_PATH, "main", "source", "bin")
        
        cmd = [
            os.path.join(rosetta_bin, "AbinitioRelax.linuxgccrelease"),
            "-database", rosetta_db,
            "-in:file:fasta", os.path.abspath(PROTEIN_FASTA),
            "-in:file:frag3", os.path.abspath(frag3_file),
            "-in:file:frag9", os.path.abspath(frag9_file),
            "-abinitio:relax",
            "-relax:fast",
            "-nstruct", str(NSTRUCT_FOLDING),
            "-out:file:silent", "folding_results.out",
            "-out:nooutput"
        ]
        
        if not run_command(cmd, cwd=f"{OUTPUT_DIR}/folding"):
            log("ERROR: Ab initio folding failed")
            return False
    
    # Select best model
    log("Selecting best model...")
    
    scores = []
    with open(silent_file, 'r') as f:
        for line in f:
            if line.startswith('SCORE:') and 'description' not in line:
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        score = float(parts[1])
                        tag = parts[-1]
                        scores.append((score, tag))
                    except (ValueError, IndexError):
                        continue
    
    if not scores:
        log("ERROR: No valid scores found")
        return False
    
    # Find best model
    best_score, best_tag = min(scores, key=lambda x: x[0])
    log(f"Best model: {best_tag} (score: {best_score:.3f})")
    
    # Extract best model
    rosetta_bin = os.path.join(ROSETTA_PATH, "main", "source", "bin")
    tag_file = f"{OUTPUT_DIR}/folding/best_tag.txt"
    
    with open(tag_file, 'w') as f:
        f.write(best_tag)
    
    cmd = [
        os.path.join(rosetta_bin, "extract_pdbs.linuxgccrelease"),
        "-in:file:silent", "folding_results.out",
        "-in:file:tagfile", "best_tag.txt",
        "-out:prefix", "best_"
    ]
    
    if not run_command(cmd, cwd=f"{OUTPUT_DIR}/folding"):
        log("ERROR: Model extraction failed")
        return False
    
    # Rename to standard name
    extracted_pdb = f"{OUTPUT_DIR}/folding/best_{best_tag}.pdb"
    if os.path.exists(extracted_pdb):
        os.rename(extracted_pdb, best_pdb)
        log("Folding completed")
        return True
    else:
        log("ERROR: Extracted model not found")
        return False

def relax_structure():
    """Relax the best folding model."""
    log("Checking relaxation...")
    
    relaxed_pdb = f"{OUTPUT_DIR}/relax/relaxed_wildtype.pdb"
    
    # Check if relaxation already done (idempotent)
    if os.path.exists(relaxed_pdb):
        log("Relaxed structure already exists, skipping")
        return True
    
    log("Relaxing structure...")
    
    best_pdb = f"{OUTPUT_DIR}/folding/best_model.pdb"
    rosetta_db = os.path.join(ROSETTA_PATH, "main", "database")
    rosetta_bin = os.path.join(ROSETTA_PATH, "main", "source", "bin")
    
    cmd = [
        os.path.join(rosetta_bin, "relax.linuxgccrelease"),
        "-database", rosetta_db,
        "-in:file:s", os.path.abspath(best_pdb),
        "-relax:fast",
        "-relax:constrain_relax_to_start_coords",
        "-relax:ramp_constraints", "false",
        "-ex1", "-ex2",
        "-use_input_sc",
        "-flip_HNQ",
        "-no_optH", "false",
        "-nstruct", str(NSTRUCT_RELAX),
        "-out:prefix", "relaxed_",
        "-out:suffix", "_wt"
    ]
    
    if not run_command(cmd, cwd=f"{OUTPUT_DIR}/relax"):
        log("ERROR: Structure relaxation failed")
        return False
    
    # Find best relaxed structure
    relax_files = [f for f in os.listdir(f"{OUTPUT_DIR}/relax") 
                   if f.startswith("relaxed_") and f.endswith("_wt.pdb")]
    
    if relax_files:
        best_relaxed = f"{OUTPUT_DIR}/relax/{relax_files[0]}"
        os.rename(best_relaxed, relaxed_pdb)
        log("Structure relaxation completed")
        return True
    else:
        log("ERROR: No relaxed structures found")
        return False

def parse_mutations():
    """Parse mutations file and create ddg_monomer input files."""
    log("Parsing mutations...")
    
    mutation_files = []
    
    with open(MUTATIONS_FILE, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Parse mutations from this line
            mutation_strings = [m.strip() for m in line.split(',')]
            mutations = []
            
            try:
                for mut_str in mutation_strings:
                    if len(mut_str) < 3:
                        raise ValueError(f"Invalid mutation format: {mut_str}")
                    
                    wt_aa = mut_str[0]
                    mut_aa = mut_str[-1]
                    position = mut_str[1:-1]
                    
                    # Validate
                    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
                    if wt_aa not in valid_aa or mut_aa not in valid_aa:
                        raise ValueError(f"Invalid amino acids in: {mut_str}")
                    
                    int(position)  # Check position is numeric
                    mutations.append((wt_aa, position, mut_aa))
                
                # Create mutation file
                mut_file = f"{OUTPUT_DIR}/ddg/mutation_{line_num}.txt"
                label_file = f"{OUTPUT_DIR}/ddg/mutation_{line_num}_label.txt"
                
                with open(mut_file, 'w') as mf:
                    mf.write(f"{len(mutations)}\n")
                    for wt_aa, pos, mut_aa in mutations:
                        mf.write(f"{wt_aa} {pos} {mut_aa}\n")
                
                with open(label_file, 'w') as lf:
                    lf.write(line)
                
                mutation_files.append(line_num)
                log(f"Created mutation file {line_num}: {line}")
                
            except (ValueError, IndexError) as e:
                log(f"WARNING: Skipping invalid mutation on line {line_num}: {e}")
                continue
    
    log(f"Parsed {len(mutation_files)} mutation sets")
    return len(mutation_files) > 0

def calculate_ddg():
    """Calculate ΔΔG for all mutations."""
    log("Calculating ΔΔG values...")
    
    relaxed_pdb = f"{OUTPUT_DIR}/relax/relaxed_wildtype.pdb"
    rosetta_db = os.path.join(ROSETTA_PATH, "main", "database")
    rosetta_bin = os.path.join(ROSETTA_PATH, "main", "source", "bin")
    
    # Find mutation files
    mutation_files = [f for f in os.listdir(f"{OUTPUT_DIR}/ddg") 
                     if f.startswith("mutation_") and f.endswith(".txt") and "_label" not in f]
    
    failed_count = 0
    
    for mut_file in sorted(mutation_files):
        mut_num = mut_file.replace("mutation_", "").replace(".txt", "")
        log_file = f"{OUTPUT_DIR}/ddg/ddg_{mut_num}.log"
        mut_pdb = f"{OUTPUT_DIR}/ddg/mut_{mut_num}_ddg_mut_1.pdb"
        
        # Check if already done (idempotent)
        if os.path.exists(log_file) and os.path.exists(mut_pdb):
            log(f"ΔΔG for mutation {mut_num} already calculated, skipping")
            continue
        
        # Read mutation label
        label_file = f"{OUTPUT_DIR}/ddg/mutation_{mut_num}_label.txt"
        try:
            with open(label_file, 'r') as f:
                mut_label = f.read().strip()
        except:
            mut_label = f"mutation_{mut_num}"
        
        log(f"Processing mutation: {mut_label}")
        
        # Run ddg_monomer
        cmd = [
            os.path.join(rosetta_bin, "ddg_monomer.linuxgccrelease"),
            "-database", rosetta_db,
            "-in:file:s", os.path.abspath(relaxed_pdb),
            "-ddg:mut_file", mut_file,
            "-ddg:iterations", str(DDG_ITERATIONS),
            "-ddg:dump_pdbs", "true",
            "-ddg:local_opt_only", "false",
            "-ddg:min_cst", "true",
            "-ddg:mean", "false",
            "-ddg:min", "true",
            "-ddg:sc_min_only", "false",
            "-ddg:ramp_repulsive", "true",
            "-out:prefix", f"mut_{mut_num}_"
        ]
        
        # Run with output to log file
        try:
            with open(f"{OUTPUT_DIR}/ddg/ddg_{mut_num}.log", 'w') as log_f:
                result = subprocess.run(cmd, cwd=f"{OUTPUT_DIR}/ddg", 
                                      stdout=log_f, stderr=subprocess.STDOUT, text=True)
            if result.returncode != 0:
                log(f"WARNING: ΔΔG calculation failed for {mut_label}")
                failed_count += 1
        except Exception as e:
            log(f"WARNING: ΔΔG calculation failed for {mut_label}: {e}")
            failed_count += 1
    
    if failed_count > 0:
        log(f"WARNING: {failed_count} ΔΔG calculations failed")
    
    log("ΔΔG calculations completed")
    return True

def aggregate_results():
    """Aggregate ΔΔG results into final table."""
    log("Aggregating results...")
    
    results_file = f"{OUTPUT_DIR}/ddg/ddg_results.txt"
    
    # Check if results already aggregated (idempotent)
    if os.path.exists(results_file):
        log("Results already aggregated, skipping")
        return True
    
    # Find label files
    label_files = [f for f in os.listdir(f"{OUTPUT_DIR}/ddg") 
                   if f.startswith("mutation_") and f.endswith("_label.txt")]
    
    results = []
    
    for label_file in sorted(label_files):
        mut_num = label_file.replace("mutation_", "").replace("_label.txt", "")
        
        # Read mutation label
        with open(f"{OUTPUT_DIR}/ddg/{label_file}", 'r') as f:
            mutation_label = f.read().strip()
        
        # Extract ΔΔG from log file
        log_file = f"{OUTPUT_DIR}/ddg/ddg_{mut_num}.log"
        ddg_value = "N/A"
        
        if os.path.exists(log_file):
            try:
                with open(log_file, 'r') as f:
                    for line in f:
                        if 'ddG' in line and ':' in line:
                            match = re.search(r'ddG\s*:?\s*([+-]?\d+\.?\d*)', line)
                            if match:
                                ddg_value = float(match.group(1))
                                break
                        elif line.strip().startswith('ddG'):
                            parts = line.strip().split()
                            if len(parts) >= 2:
                                try:
                                    ddg_value = float(parts[1])
                                    break
                                except ValueError:
                                    continue
            except Exception:
                pass
        
        results.append((mutation_label, ddg_value))
    
    # Write results
    with open(results_file, 'w') as f:
        f.write("Mutation\tddG_REU\n")
        for mutation, ddg in results:
            f.write(f"{mutation}\t{ddg}\n")
    
    # Print summary
    valid_ddgs = [ddg for _, ddg in results if isinstance(ddg, (int, float))]
    if valid_ddgs:
        log(f"Results summary:")
        log(f"  Total mutations: {len(results)}")
        log(f"  Successful: {len(valid_ddgs)}")
        log(f"  Mean ΔΔG: {sum(valid_ddgs)/len(valid_ddgs):.3f} REU")
        stabilizing = sum(1 for ddg in valid_ddgs if ddg < 0)
        destabilizing = sum(1 for ddg in valid_ddgs if ddg > 0)
        log(f"  Stabilizing: {stabilizing}, Destabilizing: {destabilizing}")
    
    log(f"Results saved to: {results_file}")
    return True

def main():
    """Main pipeline execution."""
    # Need to declare global before using in defaults
    global ROSETTA_PATH, PROTEIN_FASTA, MUTATIONS_FILE, OUTPUT_DIR, NSTRUCT_FOLDING, DDG_ITERATIONS
    
    parser = argparse.ArgumentParser(description="Minimal Rosetta protein folding and mutation pipeline")
    parser.add_argument('--rosetta-path', help='Path to Rosetta installation', default=ROSETTA_PATH)
    parser.add_argument('--protein-fasta', help='Protein FASTA file', default=PROTEIN_FASTA)
    parser.add_argument('--mutations-file', help='Mutations file', default=MUTATIONS_FILE)
    parser.add_argument('--output-dir', help='Output directory', default=OUTPUT_DIR)
    parser.add_argument('--nstruct-folding', type=int, help='Number of folding structures', default=NSTRUCT_FOLDING)
    parser.add_argument('--ddg-iterations', type=int, help='DDG iterations', default=DDG_ITERATIONS)
    
    args = parser.parse_args()
    
    # Update globals with arguments
    ROSETTA_PATH = args.rosetta_path
    PROTEIN_FASTA = args.protein_fasta
    MUTATIONS_FILE = args.mutations_file
    OUTPUT_DIR = args.output_dir
    NSTRUCT_FOLDING = args.nstruct_folding
    DDG_ITERATIONS = args.ddg_iterations
    
    log("Starting Rosetta pipeline")
    log(f"Rosetta path: {ROSETTA_PATH}")
    log(f"Protein: {PROTEIN_FASTA}")
    log(f"Mutations: {MUTATIONS_FILE}")
    
    # Run pipeline steps
    steps = [
        ("Validation", validate_setup),
        ("Directory setup", setup_directories),
        ("Fragment generation", generate_fragments),
        ("Ab initio folding", run_folding),
        ("Structure relaxation", relax_structure),
        ("Mutation parsing", parse_mutations),
        ("ΔΔG calculation", calculate_ddg),
        ("Result aggregation", aggregate_results)
    ]
    
    start_time = datetime.now()
    
    for step_name, step_func in steps:
        log(f"=== {step_name} ===")
        if not step_func():
            log(f"ERROR: {step_name} failed")
            sys.exit(1)
    
    end_time = datetime.now()
    duration = end_time - start_time
    
    log("=" * 40)
    log("Pipeline completed successfully!")
    log(f"Total time: {duration}")
    log(f"Results in: {OUTPUT_DIR}")
    log("=" * 40)

if __name__ == "__main__":
    main()