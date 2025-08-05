#!/usr/bin/env python3

import subprocess
import logging
import sys
from pathlib import Path
import os
import csv
import shutil

# Configuration constants
ROSETTA = Path("/usr/local/rosetta.source.release-371/main/source/bin")
INPUT = Path("/path/to/input")
OUTPUT = Path("./output")
THREADS = 1

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('pipeline.log')
    ]
)
logger = logging.getLogger(__name__)

def run_rosetta_command(executable, args, cwd=None):
    """Execute a Rosetta command and return success status."""
    cmd = [str(ROSETTA / executable)] + args
    logger.info(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd, 
            cwd=cwd, 
            capture_output=True, 
            text=True, 
            check=True
        )
        logger.info(f"Command completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with exit code {e.returncode}")
        logger.error(f"stderr: {e.stderr}")
        return False

def check_skip_step(step_dir):
    """Check if step should be skipped and log accordingly."""
    if step_dir.exists():
        logger.info(f"[skip] {step_dir.name}")
        return True
    return False

def parse_mutations_line(line):
    """Parse a line of comma-separated mutations into list of (pos, old, new) tuples."""
    mutations = []
    for mutation in line.strip().split(','):
        mutation = mutation.strip()
        if len(mutation) >= 3:
            old = mutation[0]
            new = mutation[-1]
            pos = mutation[1:-1]
            mutations.append((int(pos), old, new))
    return mutations

def create_mutfile(mutations, output_path):
    """Create a Rosetta .mutfile from a list of mutations."""
    with open(output_path, 'w') as f:
        f.write(f"total {len(mutations)}\n")
        f.write(f"{len(mutations)}\n")
        for pos, old, new in mutations:
            f.write(f"{old} {pos} {new}\n")

def extract_ddg_from_output(output_file):
    """Extract ddG value from Rosetta ddg_monomer output."""
    try:
        with open(output_file, 'r') as f:
            content = f.read()
            # Look for ddG value in output
            for line in content.split('\n'):
                if 'ddG' in line and 'REU' in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == 'ddG' and i + 1 < len(parts):
                            return float(parts[i + 1])
        return None
    except (FileNotFoundError, ValueError):
        return None

def run_pymol_rendering(pdb_file, output_png):
    """Generate PyMOL rendering script and execute it."""
    script_content = f"""
load {pdb_file}
hide everything
show cartoon
color purple, ss h
color yellow, ss s
color white, ss l
set cartoon_fancy_helices, 1
set ray_opaque_background, on
ray 800, 800
png {output_png}
quit
"""
    
    script_file = pdb_file.parent / f"{pdb_file.stem}_render.pml"
    with open(script_file, 'w') as f:
        f.write(script_content)
    
    try:
        result = subprocess.run(['pymol', '-qxci', str(script_file)], 
                              capture_output=True, text=True, check=True)
        script_file.unlink()  # Clean up script file
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        if script_file.exists():
            script_file.unlink()
        return False

def step_01_fold(fasta_file):
    """Step 1: Fold the FASTA sequence using AbinitioRelax."""
    step_dir = OUTPUT / "01_fold"
    
    if check_skip_step(step_dir):
        expected_pdb = step_dir / "folded.pdb"
        if expected_pdb.exists():
            return expected_pdb
        else:
            logger.warning(f"Expected output {expected_pdb} not found, re-running step")
    
    step_dir.mkdir(parents=True)
    logger.info("Starting step 01_fold")
    
    # For simplicity, we'll use a basic folding approach without fragments
    # In production, you would need fragment files
    args = [
        "-in:file:fasta", str(fasta_file),
        "-database", str(ROSETTA.parent.parent / "database"),
        "-out:nstruct", "1",
        "-out:path:pdb", str(step_dir),
        "-out:prefix", "folded_",
        "-relax:fast"
    ]
    
    success = run_rosetta_command("AbinitioRelax.default.linuxgccrelease", args, cwd=step_dir)
    
    if success:
        # Find the generated PDB file
        pdb_files = list(step_dir.glob("folded_*.pdb"))
        if pdb_files:
            # Rename to standard name
            output_pdb = step_dir / "folded.pdb"
            pdb_files[0].rename(output_pdb)
            logger.info(f"Folding completed: {output_pdb}")
            return output_pdb
    
    logger.error("Folding failed")
    return None

def step_02_relax(input_pdb):
    """Step 2: Relax the folded structure with coordinate constraints."""
    step_dir = OUTPUT / "02_relax"
    
    if check_skip_step(step_dir):
        expected_pdb = step_dir / "relaxed.pdb"
        if expected_pdb.exists():
            return expected_pdb
        else:
            logger.warning(f"Expected output {expected_pdb} not found, re-running step")
    
    step_dir.mkdir(parents=True)
    logger.info("Starting step 02_relax")
    
    args = [
        "-s", str(input_pdb),
        "-database", str(ROSETTA.parent.parent / "database"),
        "-constrain_relax_to_start_coords",
        "-relax:coord_constrain_sidechains",
        "-relax:coord_cst_stdev", "0.5",
        "-relax:ramp_constraints", "false",
        "-nstruct", "1",
        "-relax:cartesian",
        "-score:weights", "ref2015_cart",
        "-relax:min_type", "lbfgs_armijo_nonmonotone",
        "-fa_max_dis", "9.0",
        "-out:path:pdb", str(step_dir),
        "-out:prefix", "relaxed_"
    ]
    
    success = run_rosetta_command("relax.default.linuxgccrelease", args, cwd=step_dir)
    
    if success:
        # Find the generated PDB file
        pdb_files = list(step_dir.glob("relaxed_*.pdb"))
        if pdb_files:
            # Rename to standard name
            output_pdb = step_dir / "relaxed.pdb"
            pdb_files[0].rename(output_pdb)
            logger.info(f"Relaxation completed: {output_pdb}")
            return output_pdb
    
    logger.error("Relaxation failed")
    return None

def step_03_mutfiles(mutations_file):
    """Step 3: Generate .mutfile for each line of mutations."""
    step_dir = OUTPUT / "03_mutfiles"
    
    if check_skip_step(step_dir):
        existing_mutfiles = list(step_dir.glob("*.mutfile"))
        if existing_mutfiles:
            return existing_mutfiles
        else:
            logger.warning(f"No mutfiles found in {step_dir}, re-running step")
    
    step_dir.mkdir(parents=True)
    logger.info("Starting step 03_mutfiles")
    
    mutfiles = []
    with open(mutations_file, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
                
            mutations = parse_mutations_line(line)
            if mutations:
                mutfile_path = step_dir / f"mutations_{i+1:03d}.mutfile"
                create_mutfile(mutations, mutfile_path)
                mutfiles.append(mutfile_path)
                logger.info(f"Created mutfile: {mutfile_path}")
    
    logger.info(f"Generated {len(mutfiles)} mutfiles")
    return mutfiles

def step_04_ddg(relaxed_pdb, mutfiles):
    """Step 4: Run ddg_monomer for each mutfile and aggregate results."""
    step_dir = OUTPUT / "04_ddg" 
    
    if check_skip_step(step_dir):
        expected_csv = step_dir / "ddg.csv"
        if expected_csv.exists():
            return expected_csv
        else:
            logger.warning(f"Expected output {expected_csv} not found, re-running step")
    
    step_dir.mkdir(parents=True)
    logger.info("Starting step 04_ddg")
    
    # Initialize CSV file
    csv_file = step_dir / "ddg.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["mutation", "ddg_REU"])
    
    for i, mutfile in enumerate(mutfiles):
        # Create subdirectory for this mutation set
        mut_dir = step_dir / f"mut_{i+1:03d}"
        mut_dir.mkdir(exist_ok=True)
        
        # Read mutation description for CSV
        mutations = []
        with open(mutfile, 'r') as f:
            lines = f.readlines()
            for line in lines[2:]:  # Skip 'total' and count lines
                parts = line.strip().split()
                if len(parts) >= 3:
                    mutations.append(f"{parts[0]}{parts[1]}{parts[2]}")
        
        mutation_desc = ",".join(mutations)
        logger.info(f"Processing mutations: {mutation_desc}")
        
        # Run ddg_monomer
        args = [
            "-in:file:s", str(relaxed_pdb),
            "-ddg::mut_file", str(mutfile),
            "-ddg::iterations", "3",  # As specified
            "-database", str(ROSETTA.parent.parent / "database"),
            "-score:weights", "ref2015_cart",
            "-relax:cartesian", "true",
            "-fa_max_dis", "9.0",
            "-out:path:pdb", str(mut_dir),
            "-out:prefix", f"mut_{i+1:03d}_"
        ]
        
        success = run_rosetta_command("ddg_monomer.default.linuxgccrelease", args, cwd=mut_dir)
        
        if success:
            # Extract ddG value from output
            score_file = mut_dir / "ddg_predictions.out"
            ddg_value = extract_ddg_from_output(score_file)
            
            if ddg_value is None:
                # Try alternative output files
                for out_file in mut_dir.glob("*.out"):
                    ddg_value = extract_ddg_from_output(out_file)
                    if ddg_value is not None:
                        break
            
            if ddg_value is None:
                ddg_value = "N/A"
                logger.warning(f"Could not extract ddG value for {mutation_desc}")
            
            # Append to CSV
            with open(csv_file, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([mutation_desc, ddg_value])
            
            logger.info(f"Completed ddG calculation for {mutation_desc}: {ddg_value}")
        else:
            logger.error(f"ddG calculation failed for {mutation_desc}")
            # Still add to CSV with error marker
            with open(csv_file, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([mutation_desc, "ERROR"])
    
    logger.info(f"ddG calculations completed, results in {csv_file}")
    return csv_file

def step_05_superpose(native_pdb):
    """Step 5: Superimpose all mutant PDBs onto the relaxed native."""
    step_dir = OUTPUT / "05_superpose"
    ddg_dir = OUTPUT / "04_ddg"
    
    if check_skip_step(step_dir):
        return step_dir
    
    step_dir.mkdir(parents=True)
    logger.info("Starting step 05_superpose")
    
    # Find all mutant PDB files
    mutant_pdbs = []
    for mut_subdir in ddg_dir.glob("mut_*"):
        if mut_subdir.is_dir():
            for pdb_file in mut_subdir.glob("*_relaxed*.pdb"):
                mutant_pdbs.append(pdb_file)
            # Also look for other output PDB patterns
            for pdb_file in mut_subdir.glob("mut_*.pdb"):
                if pdb_file not in mutant_pdbs:
                    mutant_pdbs.append(pdb_file)
    
    if not mutant_pdbs:
        logger.warning("No mutant PDB files found for superposition")
        return step_dir
    
    for mutant_pdb in mutant_pdbs:
        # Create corresponding output directory
        output_subdir = step_dir / mutant_pdb.parent.name
        output_subdir.mkdir(exist_ok=True)
        
        output_pdb = output_subdir / f"superposed_{mutant_pdb.name}"
        
        logger.info(f"Superposing {mutant_pdb.name} onto native structure")
        
        # Use Rosetta's structural alignment tool
        args = [
            "-in:file:s", str(mutant_pdb),
            "-in:file:native", str(native_pdb),
            "-database", str(ROSETTA.parent.parent / "database"),
            "-out:path:pdb", str(output_subdir),
            "-out:prefix", "superposed_",
            "-out:file:scorefile", str(output_subdir / "scores.sc")
        ]
        
        # Try using score_jd2 with superimpose mover
        success = run_rosetta_command("score_jd2.default.linuxgccrelease", args, cwd=output_subdir)
        
        if not success:
            # Fallback: simple copy if superposition fails
            logger.warning(f"Superposition failed for {mutant_pdb.name}, copying original")
            shutil.copy2(mutant_pdb, output_pdb)
    
    logger.info("Superposition step completed")
    return step_dir

def step_06_images(native_pdb):
    """Step 6: Generate PyMOL renderings if PyMOL is available."""
    step_dir = OUTPUT / "06_images"
    
    if check_skip_step(step_dir):
        return step_dir
    
    step_dir.mkdir(parents=True)
    logger.info("Starting step 06_images")
    
    # Check if PyMOL is available
    try:
        subprocess.run(['pymol', '-h'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.warning("PyMOL not available, skipping image generation")
        return step_dir
    
    # Render native structure
    native_png = step_dir / "native_relaxed.png"
    if run_pymol_rendering(native_pdb, native_png):
        logger.info(f"Rendered native structure: {native_png}")
    else:
        logger.error("Failed to render native structure")
    
    # Render all mutant structures
    superpose_dir = OUTPUT / "05_superpose"
    if superpose_dir.exists():
        for mut_subdir in superpose_dir.glob("mut_*"):
            if mut_subdir.is_dir():
                for pdb_file in mut_subdir.glob("*.pdb"):
                output_png = step_dir / f"{mut_subdir.name}_{pdb_file.stem}.png"
                if run_pymol_rendering(pdb_file, output_png):
                    logger.info(f"Rendered mutant structure: {output_png}")
                else:
                    logger.error(f"Failed to render {pdb_file}")
    
    # Also render any structures from 04_ddg if superposition wasn't done
    ddg_dir = OUTPUT / "04_ddg"
    if ddg_dir.exists():
        for mut_subdir in ddg_dir.glob("mut_*"):
            if mut_subdir.is_dir():
                for pdb_file in mut_subdir.glob("*_relaxed*.pdb"):
                output_png = step_dir / f"{mut_subdir.name}_{pdb_file.stem}.png"
                if not output_png.exists():  # Don't overwrite superposed versions
                    if run_pymol_rendering(pdb_file, output_png):
                        logger.info(f"Rendered mutant structure: {output_png}")
                    else:
                        logger.error(f"Failed to render {pdb_file}")
    
    logger.info("Image generation completed")
    return step_dir

def main():
    """Main pipeline execution function."""
    # Input file paths
    fasta_file = INPUT / "protein.faa"
    mutations_file = INPUT / "mutations.txt"
    
    # Validate input files
    if not fasta_file.exists():
        logger.error(f"Input FASTA file not found: {fasta_file}")
        sys.exit(1)
    
    if not mutations_file.exists():
        logger.error(f"Input mutations file not found: {mutations_file}")
        sys.exit(1)
    
    # Create output directory
    OUTPUT.mkdir(parents=True, exist_ok=True)
    
    logger.info("Starting Rosetta 3.14 monomer ΔΔG pipeline")
    logger.info(f"Input FASTA: {fasta_file}")
    logger.info(f"Input mutations: {mutations_file}")
    logger.info(f"Output directory: {OUTPUT}")
    
    # Step 1: Fold the FASTA sequence
    folded_pdb = step_01_fold(fasta_file)
    if not folded_pdb:
        logger.error("Pipeline failed at folding step")
        sys.exit(1)
    
    # Step 2: Relax the folded structure
    relaxed_pdb = step_02_relax(folded_pdb)
    if not relaxed_pdb:
        logger.error("Pipeline failed at relaxation step")
        sys.exit(1)
    
    # Step 3: Generate mutfiles
    mutfiles = step_03_mutfiles(mutations_file)
    if not mutfiles:
        logger.error("Pipeline failed at mutfile generation step")
        sys.exit(1)
    
    # Step 4: Run ddG calculations
    ddg_csv = step_04_ddg(relaxed_pdb, mutfiles)
    if not ddg_csv:
        logger.error("Pipeline failed at ddG calculation step")
        sys.exit(1)
    
    # Step 5: Superpose structures
    superpose_dir = step_05_superpose(relaxed_pdb)
    
    # Step 6: Generate images (optional)
    images_dir = step_06_images(relaxed_pdb)
    
    logger.info("done")

if __name__ == "__main__":
    main()