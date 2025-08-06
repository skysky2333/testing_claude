#!/usr/bin/env python3
"""
Self-contained idempotent Rosetta protein folding and mutation analysis pipeline.
Consolidates all functionality into a single script with resume capability.

Author: Claude Code Assistant
"""

import os
import sys
import json
import yaml
import argparse
import subprocess
import logging
import re
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any
import pandas as pd
import numpy as np

class RosettaPipelineConfig:
    """Configuration management for the Rosetta pipeline."""
    
    DEFAULT_CONFIG = {
        # Rosetta installation
        "rosetta_path": "/usr/local/rosetta.source.release-371",
        
        # Input files
        "protein_fasta": "protein.faa",
        "mutations_file": "mutations.txt",
        
        # Output directory
        "output_dir": "rosetta_results",
        
        # Folding parameters
        "nstruct_folding": 1000,
        "nstruct_relax": 10,
        "ddg_iterations": 50,
        
        # Computational resources
        "threads_folding": 4,
        "threads_ddg": 2,
        
        # Advanced options
        "fragment_picking": {
            "n_candidates": 200,
            "n_frags": 25,
            "frag_sizes": [3, 9]
        },
        
        "relaxation": {
            "fast_mode": True,
            "constrain_coords": True,
            "ramp_constraints": False
        },
        
        "ddg_calculation": {
            "local_opt_only": False,
            "min_cst": True,
            "dump_pdbs": True,
            "ramp_repulsive": True
        },
        
        "visualization": {
            "image_width": 800,
            "image_height": 600,
            "image_dpi": 300,
            "ray_tracing": True
        },
        
        # Pipeline options
        "save_intermediate": True,
        "compress_output": False,
        "skip_visualization": False
    }
    
    def __init__(self, config_file: Optional[str] = None, cli_overrides: Dict[str, Any] = None):
        """Initialize configuration with defaults, file, and CLI overrides."""
        self.config = self.DEFAULT_CONFIG.copy()
        
        # Load from config file if provided
        if config_file and os.path.exists(config_file):
            self._load_config_file(config_file)
        
        # Apply CLI overrides
        if cli_overrides:
            self._apply_overrides(cli_overrides)
        
        # Compute derived paths
        self._compute_paths()
    
    def _load_config_file(self, config_file: str):
        """Load configuration from YAML or JSON file."""
        try:
            with open(config_file, 'r') as f:
                if config_file.endswith('.yaml') or config_file.endswith('.yml'):
                    file_config = yaml.safe_load(f)
                else:
                    file_config = json.load(f)
            
            self._deep_update(self.config, file_config)
        except Exception as e:
            logging.warning(f"Failed to load config file {config_file}: {e}")
    
    def _apply_overrides(self, overrides: Dict[str, Any]):
        """Apply CLI argument overrides to configuration."""
        for key, value in overrides.items():
            if value is not None:
                self._set_nested_value(self.config, key, value)
    
    def _deep_update(self, base_dict: Dict, update_dict: Dict):
        """Deep update of nested dictionary."""
        for key, value in update_dict.items():
            if isinstance(value, dict) and key in base_dict and isinstance(base_dict[key], dict):
                self._deep_update(base_dict[key], value)
            else:
                base_dict[key] = value
    
    def _set_nested_value(self, config: Dict, key_path: str, value: Any):
        """Set nested configuration value using dot notation."""
        keys = key_path.split('.')
        current = config
        for key in keys[:-1]:
            if key not in current:
                current[key] = {}
            current = current[key]
        current[keys[-1]] = value
    
    def _compute_paths(self):
        """Compute derived paths from configuration."""
        self.rosetta_bin = os.path.join(self.config["rosetta_path"], "main", "source", "bin")
        self.rosetta_db = os.path.join(self.config["rosetta_path"], "main", "database")
        
        output_dir = self.config["output_dir"]
        self.fragments_dir = os.path.join(output_dir, "fragments")
        self.folding_dir = os.path.join(output_dir, "folding")
        self.relax_dir = os.path.join(output_dir, "relax")
        self.ddg_dir = os.path.join(output_dir, "ddg")
        self.images_dir = os.path.join(output_dir, "images")
        self.state_file = os.path.join(output_dir, "pipeline_state.json")
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get configuration value using dot notation."""
        keys = key.split('.')
        current = self.config
        for key in keys:
            if isinstance(current, dict) and key in current:
                current = current[key]
            else:
                return default
        return current
    
    def __getitem__(self, key: str) -> Any:
        """Allow dictionary-style access."""
        return self.get(key)

class RosettaPipeline:
    """Main pipeline class for Rosetta protein folding and mutation analysis."""
    
    def __init__(self, config: RosettaPipelineConfig):
        """Initialize pipeline with configuration."""
        self.config = config
        self.logger = self._setup_logging()
        self.state = self._load_state()
        
        # Track protein name
        self.protein_name = self._get_protein_name()
        
    def _setup_logging(self) -> logging.Logger:
        """Set up logging configuration."""
        logger = logging.getLogger('RosettaPipeline')
        logger.setLevel(logging.INFO)
        
        # Create console handler
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        
        return logger
    
    def _load_state(self) -> Dict[str, Any]:
        """Load pipeline state from file."""
        if os.path.exists(self.config.state_file):
            try:
                with open(self.config.state_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                self.logger.warning(f"Failed to load state file: {e}")
        
        return {
            "completed_steps": [],
            "last_run": None,
            "protein_name": None
        }
    
    def _save_state(self):
        """Save pipeline state to file."""
        self.state["last_run"] = datetime.now().isoformat()
        self.state["protein_name"] = self.protein_name
        
        os.makedirs(os.path.dirname(self.config.state_file), exist_ok=True)
        try:
            with open(self.config.state_file, 'w') as f:
                json.dump(self.state, f, indent=2)
        except Exception as e:
            self.logger.warning(f"Failed to save state file: {e}")
    
    def _mark_step_complete(self, step_name: str):
        """Mark a pipeline step as completed."""
        if step_name not in self.state["completed_steps"]:
            self.state["completed_steps"].append(step_name)
        self._save_state()
    
    def _is_step_complete(self, step_name: str) -> bool:
        """Check if a pipeline step is already completed."""
        return step_name in self.state["completed_steps"]
    
    def _get_protein_name(self) -> str:
        """Extract protein name from FASTA file."""
        fasta_file = self.config["protein_fasta"]
        if os.path.exists(fasta_file):
            try:
                with open(fasta_file, 'r') as f:
                    header = f.readline().strip()
                    if header.startswith('>'):
                        name = header[1:].split()[0].replace('|', '_').replace('/', '_')
                        return name if name else "protein"
            except Exception:
                pass
        return "protein"
    
    def _run_command(self, cmd: List[str], cwd: str = None, log_file: str = None) -> bool:
        """Run a command with proper error handling and logging."""
        self.logger.info(f"Running command: {' '.join(cmd)}")
        
        try:
            with open(log_file, 'w') if log_file else open(os.devnull, 'w') as log_handle:
                result = subprocess.run(
                    cmd,
                    cwd=cwd,
                    stdout=log_handle if log_file else subprocess.PIPE,
                    stderr=subprocess.STDOUT if log_file else subprocess.PIPE,
                    text=True
                )
            
            if result.returncode != 0:
                self.logger.error(f"Command failed with return code {result.returncode}")
                if log_file:
                    self.logger.error(f"Check log file: {log_file}")
                return False
            
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to run command: {e}")
            return False
    
    def validate_dependencies(self) -> bool:
        """Validate that all required dependencies are available."""
        self.logger.info("Validating dependencies...")
        
        # Check Rosetta installation
        if not os.path.exists(self.config["rosetta_path"]):
            self.logger.error(f"Rosetta installation not found: {self.config['rosetta_path']}")
            return False
        
        if not os.path.exists(self.config.rosetta_bin):
            self.logger.error(f"Rosetta bin directory not found: {self.config.rosetta_bin}")
            return False
        
        if not os.path.exists(self.config.rosetta_db):
            self.logger.error(f"Rosetta database not found: {self.config.rosetta_db}")
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
            exe_path = os.path.join(self.config.rosetta_bin, exe)
            if not os.path.exists(exe_path):
                # Try alternative extensions
                alt_exes = [exe.replace('.linuxgccrelease', ext) 
                           for ext in ['.macosclangrelease', '.default.linuxgccrelease']]
                
                found = False
                for alt_exe in alt_exes:
                    if os.path.exists(os.path.join(self.config.rosetta_bin, alt_exe)):
                        found = True
                        break
                
                if not found:
                    self.logger.error(f"Required Rosetta executable not found: {exe}")
                    return False
        
        # Check input files
        if not os.path.exists(self.config["protein_fasta"]):
            self.logger.error(f"Protein FASTA file not found: {self.config['protein_fasta']}")
            return False
        
        if not os.path.exists(self.config["mutations_file"]):
            self.logger.error(f"Mutations file not found: {self.config['mutations_file']}")
            return False
        
        # Check Python modules
        try:
            import pandas
            import numpy
        except ImportError as e:
            self.logger.error(f"Required Python module not found: {e}")
            return False
        
        # Check optional PyMOL
        if not self.config["skip_visualization"]:
            try:
                import pymol
                self.logger.info("PyMOL module available for visualization")
            except ImportError:
                self.logger.warning("PyMOL not available - visualization will be skipped")
                self.config["skip_visualization"] = True
        
        self.logger.info("All dependencies validated successfully")
        return True
    
    def setup_directories(self) -> bool:
        """Set up the directory structure for pipeline outputs."""
        if self._is_step_complete("setup_directories"):
            self.logger.info("Directory setup already completed, skipping...")
            return True
        
        self.logger.info("Setting up directory structure...")
        
        try:
            directories = [
                self.config["output_dir"],
                self.config.fragments_dir,
                self.config.folding_dir,
                self.config.relax_dir,
                self.config.ddg_dir,
                self.config.images_dir
            ]
            
            for directory in directories:
                os.makedirs(directory, exist_ok=True)
            
            self._mark_step_complete("setup_directories")
            self.logger.info("Directory structure created successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to create directories: {e}")
            return False
    
    def generate_fragments(self) -> bool:
        """Generate fragment files for ab initio folding."""
        if self._is_step_complete("generate_fragments"):
            self.logger.info("Fragment generation already completed, skipping...")
            return True
        
        self.logger.info("Generating fragment files...")
        
        frag3_file = os.path.join(self.config.fragments_dir, f"{self.protein_name}.200.3mers")
        frag9_file = os.path.join(self.config.fragments_dir, f"{self.protein_name}.200.9mers")
        config_file = os.path.join(self.config.fragments_dir, "fragment_picker_config")
        log_file = os.path.join(self.config.fragments_dir, "fragment_picker.log")
        
        # Check if fragment files already exist
        if os.path.exists(frag3_file) and os.path.exists(frag9_file):
            self.logger.info("Fragment files already exist, skipping generation...")
            self._mark_step_complete("generate_fragments")
            return True
        
        try:
            # Create fragment picker configuration
            frag_config = self.config.get("fragment_picking", {})
            config_content = f"""database {self.config.rosetta_db}
in:file:fasta {self.config["protein_fasta"]}
frags:describe_fragments {self.config.fragments_dir}/{self.protein_name}_frags.fsc
out:file:frag_prefix {self.config.fragments_dir}/{self.protein_name}
frags:n_candidates {frag_config.get('n_candidates', 200)}
frags:n_frags {frag_config.get('n_frags', 25)}
frags:frag_sizes {' '.join(map(str, frag_config.get('frag_sizes', [3, 9])))}
"""
            
            with open(config_file, 'w') as f:
                f.write(config_content)
            
            # Run fragment picker
            cmd = [
                os.path.join(self.config.rosetta_bin, "fragment_picker.linuxgccrelease"),
                f"@{config_file}"
            ]
            
            if self._run_command(cmd, cwd=self.config.fragments_dir, log_file=log_file):
                if os.path.exists(frag3_file) and os.path.exists(frag9_file):
                    self._mark_step_complete("generate_fragments")
                    self.logger.info("Fragment generation completed successfully")
                    return True
                else:
                    self.logger.error("Fragment files not created")
                    return False
            else:
                return False
                
        except Exception as e:
            self.logger.error(f"Fragment generation failed: {e}")
            return False
    
    def run_ab_initio_folding(self) -> bool:
        """Perform ab initio protein folding."""
        if self._is_step_complete("ab_initio_folding"):
            self.logger.info("Ab initio folding already completed, skipping...")
            return True
        
        self.logger.info("Running ab initio folding...")
        
        silent_file = os.path.join(self.config.folding_dir, "folding_results.out")
        log_file = os.path.join(self.config.folding_dir, "folding.log")
        
        # Check if silent file already exists
        if os.path.exists(silent_file):
            self.logger.info("Folding results already exist, skipping...")
            self._mark_step_complete("ab_initio_folding")
            return True
        
        try:
            frag3_file = os.path.join(self.config.fragments_dir, f"{self.protein_name}.200.3mers")
            frag9_file = os.path.join(self.config.fragments_dir, f"{self.protein_name}.200.9mers")
            
            cmd = [
                os.path.join(self.config.rosetta_bin, "AbinitioRelax.linuxgccrelease"),
                "-database", self.config.rosetta_db,
                "-in:file:fasta", os.path.abspath(self.config["protein_fasta"]),
                "-in:file:frag3", os.path.abspath(frag3_file),
                "-in:file:frag9", os.path.abspath(frag9_file),
                "-abinitio:relax",
                "-relax:fast",
                "-nstruct", str(self.config["nstruct_folding"]),
                "-out:file:silent", "folding_results.out",
                "-out:nooutput"
            ]
            
            if self._run_command(cmd, cwd=self.config.folding_dir, log_file=log_file):
                if os.path.exists(silent_file):
                    self._mark_step_complete("ab_initio_folding")
                    self.logger.info("Ab initio folding completed successfully")
                    return True
                else:
                    self.logger.error("Folding results file not created")
                    return False
            else:
                return False
                
        except Exception as e:
            self.logger.error(f"Ab initio folding failed: {e}")
            return False
    
    def select_best_model(self) -> bool:
        """Select the best folding model based on energy."""
        if self._is_step_complete("select_best_model"):
            self.logger.info("Best model selection already completed, skipping...")
            return True
        
        self.logger.info("Selecting best folding model...")
        
        silent_file = os.path.join(self.config.folding_dir, "folding_results.out")
        best_pdb = os.path.join(self.config.folding_dir, "best_model.pdb")
        best_score_file = os.path.join(self.config.folding_dir, "best_model_score.txt")
        
        # Check if best model already exists
        if os.path.exists(best_pdb):
            self.logger.info("Best model already exists, skipping selection...")
            self._mark_step_complete("select_best_model")
            return True
        
        try:
            # Extract scores and find best model
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
                self.logger.error("No valid scores found in silent file")
                return False
            
            # Find best (lowest energy) model
            best_score, best_tag = min(scores, key=lambda x: x[0])
            
            # Save best score info
            with open(best_score_file, 'w') as f:
                f.write(f"{best_score} {best_tag}\n")
            
            # Extract best model as PDB
            tag_file = os.path.join(self.config.folding_dir, "best_tag.txt")
            with open(tag_file, 'w') as f:
                f.write(best_tag)
            
            cmd = [
                os.path.join(self.config.rosetta_bin, "extract_pdbs.linuxgccrelease"),
                "-in:file:silent", silent_file,
                "-in:file:tagfile", tag_file,
                "-out:prefix", "best_"
            ]
            
            if self._run_command(cmd, cwd=self.config.folding_dir):
                # Rename to standard name
                extracted_pdb = os.path.join(self.config.folding_dir, f"best_{best_tag}.pdb")
                if os.path.exists(extracted_pdb):
                    shutil.move(extracted_pdb, best_pdb)
                    self._mark_step_complete("select_best_model")
                    self.logger.info(f"Best model selected: {best_tag} (score: {best_score:.3f})")
                    return True
                else:
                    self.logger.error("Failed to extract best model PDB")
                    return False
            else:
                return False
                
        except Exception as e:
            self.logger.error(f"Best model selection failed: {e}")
            return False
    
    def relax_structure(self) -> bool:
        """Relax the best folding model."""
        if self._is_step_complete("relax_structure"):
            self.logger.info("Structure relaxation already completed, skipping...")
            return True
        
        self.logger.info("Relaxing structure...")
        
        best_pdb = os.path.join(self.config.folding_dir, "best_model.pdb")
        relaxed_pdb = os.path.join(self.config.relax_dir, "relaxed_wildtype.pdb")
        log_file = os.path.join(self.config.relax_dir, "relax.log")
        
        # Check if relaxed structure already exists
        if os.path.exists(relaxed_pdb):
            self.logger.info("Relaxed structure already exists, skipping...")
            self._mark_step_complete("relax_structure")
            return True
        
        try:
            relax_config = self.config.get("relaxation", {})
            
            cmd = [
                os.path.join(self.config.rosetta_bin, "relax.linuxgccrelease"),
                "-database", self.config.rosetta_db,
                "-in:file:s", os.path.abspath(best_pdb),
                "-relax:fast" if relax_config.get("fast_mode", True) else "-relax:thorough",
                "-ex1", "-ex2",
                "-use_input_sc",
                "-flip_HNQ",
                "-no_optH", "false",
                "-nstruct", str(self.config["nstruct_relax"]),
                "-out:prefix", "relaxed_",
                "-out:suffix", "_wt"
            ]
            
            if relax_config.get("constrain_coords", True):
                cmd.extend(["-relax:constrain_relax_to_start_coords"])
            
            if not relax_config.get("ramp_constraints", False):
                cmd.extend(["-relax:ramp_constraints", "false"])
            
            if self._run_command(cmd, cwd=self.config.relax_dir, log_file=log_file):
                # Find and rename best relaxed structure
                relaxed_files = [f for f in os.listdir(self.config.relax_dir) 
                               if f.startswith("relaxed_") and f.endswith("_wt.pdb")]
                
                if relaxed_files:
                    best_relaxed = os.path.join(self.config.relax_dir, relaxed_files[0])
                    shutil.copy2(best_relaxed, relaxed_pdb)
                    self._mark_step_complete("relax_structure")
                    self.logger.info("Structure relaxation completed successfully")
                    return True
                else:
                    self.logger.error("No relaxed structures found")
                    return False
            else:
                return False
                
        except Exception as e:
            self.logger.error(f"Structure relaxation failed: {e}")
            return False
    
    def parse_mutation_string(self, mutation_str: str) -> Tuple[str, str, str]:
        """Parse a single mutation string like 'A1G' into components."""
        mutation_str = mutation_str.strip()
        if len(mutation_str) < 3:
            raise ValueError(f"Invalid mutation format: {mutation_str}")
        
        wild_type_aa = mutation_str[0]
        mutant_aa = mutation_str[-1]
        position = mutation_str[1:-1]
        
        # Validate amino acid codes
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        if wild_type_aa not in valid_aa:
            raise ValueError(f"Invalid wild-type amino acid: {wild_type_aa}")
        if mutant_aa not in valid_aa:
            raise ValueError(f"Invalid mutant amino acid: {mutant_aa}")
        
        # Validate position is numeric
        try:
            int(position)
        except ValueError:
            raise ValueError(f"Invalid position: {position}")
        
        return wild_type_aa, position, mutant_aa
    
    def create_ddg_mutation_file(self, mutations: List[Tuple[str, str, str]], output_file: str):
        """Create a mutation file in ddg_monomer format."""
        num_mutations = len(mutations)
        
        with open(output_file, 'w') as f:
            f.write(f"{num_mutations}\n")
            for wt_aa, pos, mut_aa in mutations:
                f.write(f"{wt_aa} {pos} {mut_aa}\n")
    
    def parse_mutations(self) -> bool:
        """Parse mutations file and create individual mutation files."""
        if self._is_step_complete("parse_mutations"):
            self.logger.info("Mutation parsing already completed, skipping...")
            return True
        
        self.logger.info("Parsing mutations...")
        
        mutations_file = self.config["mutations_file"]
        mutation_files = []
        labels = []
        
        try:
            with open(mutations_file, 'r') as f:
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
                            wt_aa, pos, mut_aa = self.parse_mutation_string(mut_str)
                            mutations.append((wt_aa, pos, mut_aa))
                        
                        # Create output files
                        mut_file = os.path.join(self.config.ddg_dir, f"mutation_{line_num}.txt")
                        label_file = os.path.join(self.config.ddg_dir, f"mutation_{line_num}_label.txt")
                        
                        self.create_ddg_mutation_file(mutations, mut_file)
                        
                        with open(label_file, 'w') as lf:
                            lf.write(line)
                        
                        mutation_files.append(mut_file)
                        labels.append(line)
                        
                        self.logger.debug(f"Created mutation file {line_num}: {line}")
                        
                    except ValueError as e:
                        self.logger.warning(f"Error processing line {line_num} '{line}': {e}")
                        continue
            
            if mutation_files:
                self._mark_step_complete("parse_mutations")
                self.logger.info(f"Parsed {len(mutation_files)} mutation sets successfully")
                return True
            else:
                self.logger.error("No valid mutations found")
                return False
                
        except Exception as e:
            self.logger.error(f"Mutation parsing failed: {e}")
            return False
    
    def calculate_ddg(self) -> bool:
        """Calculate ΔΔG for all mutations using ddg_monomer."""
        if self._is_step_complete("calculate_ddg"):
            self.logger.info("ΔΔG calculations already completed, skipping...")
            return True
        
        self.logger.info("Calculating ΔΔG values...")
        
        relaxed_pdb = os.path.join(self.config.relax_dir, "relaxed_wildtype.pdb")
        if not os.path.exists(relaxed_pdb):
            self.logger.error("Relaxed wildtype structure not found")
            return False
        
        # Find all mutation files
        mutation_files = [f for f in os.listdir(self.config.ddg_dir) 
                         if f.startswith("mutation_") and f.endswith(".txt") and "_label" not in f]
        
        if not mutation_files:
            self.logger.error("No mutation files found")
            return False
        
        ddg_config = self.config.get("ddg_calculation", {})
        failed_calculations = 0
        
        try:
            for mut_file in sorted(mutation_files):
                # Extract mutation number
                mut_num = mut_file.replace("mutation_", "").replace(".txt", "")
                label_file = os.path.join(self.config.ddg_dir, f"mutation_{mut_num}_label.txt")
                
                # Read mutation label
                try:
                    with open(label_file, 'r') as f:
                        mut_label = f.read().strip()
                except FileNotFoundError:
                    mut_label = f"mutation_{mut_num}"
                
                self.logger.info(f"Processing mutation: {mut_label}")
                
                # Check if already completed
                log_file = os.path.join(self.config.ddg_dir, f"ddg_{mut_num}.log")
                mut_pdb = os.path.join(self.config.ddg_dir, f"mut_{mut_num}_ddg_mut_1.pdb")
                
                if os.path.exists(log_file) and os.path.exists(mut_pdb):
                    self.logger.debug(f"ΔΔG calculation for {mut_label} already completed")
                    continue
                
                # Run ddg_monomer
                cmd = [
                    os.path.join(self.config.rosetta_bin, "ddg_monomer.linuxgccrelease"),
                    "-database", self.config.rosetta_db,
                    "-in:file:s", os.path.abspath(relaxed_pdb),
                    "-ddg:mut_file", os.path.join(self.config.ddg_dir, mut_file),
                    "-ddg:iterations", str(self.config["ddg_iterations"]),
                    "-ddg:dump_pdbs", str(ddg_config.get("dump_pdbs", True)).lower(),
                    "-ddg:local_opt_only", str(ddg_config.get("local_opt_only", False)).lower(),
                    "-ddg:min_cst", str(ddg_config.get("min_cst", True)).lower(),
                    "-ddg:mean", "false",
                    "-ddg:min", "true",
                    "-ddg:sc_min_only", "false",
                    "-ddg:ramp_repulsive", str(ddg_config.get("ramp_repulsive", True)).lower(),
                    "-out:prefix", f"mut_{mut_num}_"
                ]
                
                if not self._run_command(cmd, cwd=self.config.ddg_dir, log_file=log_file):
                    self.logger.warning(f"ΔΔG calculation failed for {mut_label}")
                    failed_calculations += 1
            
            if failed_calculations == 0:
                self._mark_step_complete("calculate_ddg")
                self.logger.info("All ΔΔG calculations completed successfully")
                return True
            else:
                self.logger.warning(f"{failed_calculations} ΔΔG calculations failed")
                # Mark as complete even with some failures to allow pipeline continuation
                self._mark_step_complete("calculate_ddg")
                return True
                
        except Exception as e:
            self.logger.error(f"ΔΔG calculation failed: {e}")
            return False
    
    def extract_ddg_from_log(self, log_file: str) -> Tuple[Optional[float], bool]:
        """Extract ΔΔG value from ddg_monomer log file."""
        if not os.path.exists(log_file):
            return None, False
        
        ddg_values = []
        
        try:
            with open(log_file, 'r') as f:
                for line in f:
                    # Look for ddG lines in the output
                    if 'ddG' in line and ':' in line:
                        # Try to extract numerical value
                        match = re.search(r'ddG\s*:?\s*([+-]?\d+\.?\d*)', line)
                        if match:
                            try:
                                value = float(match.group(1))
                                ddg_values.append(value)
                            except ValueError:
                                continue
                    
                    # Alternative patterns for ddG output
                    if line.strip().startswith('ddG'):
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            try:
                                value = float(parts[1])
                                ddg_values.append(value)
                            except (ValueError, IndexError):
                                continue
            
            if ddg_values:
                # Return the last/final ddG value found
                return ddg_values[-1], True
            else:
                return None, False
                
        except Exception as e:
            self.logger.warning(f"Error reading {log_file}: {e}")
            return None, False
    
    def aggregate_ddg_results(self) -> bool:
        """Aggregate ΔΔG results from all mutations."""
        if self._is_step_complete("aggregate_ddg_results"):
            self.logger.info("ΔΔG result aggregation already completed, skipping...")
            return True
        
        self.logger.info("Aggregating ΔΔG results...")
        
        results_file = os.path.join(self.config.ddg_dir, "ddg_results.txt")
        
        # Find all mutation label files
        label_files = [f for f in os.listdir(self.config.ddg_dir) 
                      if f.startswith("mutation_") and f.endswith("_label.txt")]
        
        if not label_files:
            self.logger.error("No mutation label files found")
            return False
        
        results = []
        
        try:
            for label_file in sorted(label_files):
                # Extract mutation number
                match = re.search(r'mutation_(\d+)_label\.txt', label_file)
                if not match:
                    continue
                
                mut_num = match.group(1)
                label_path = os.path.join(self.config.ddg_dir, label_file)
                
                # Read mutation label
                try:
                    with open(label_path, 'r') as f:
                        mutation_label = f.read().strip()
                except Exception:
                    mutation_label = f'mutation_{mut_num}'
                
                # Find corresponding log file
                log_file = os.path.join(self.config.ddg_dir, f"ddg_{mut_num}.log")
                ddg_value, success = self.extract_ddg_from_log(log_file)
                
                results.append({
                    'Mutation': mutation_label,
                    'ddG_REU': ddg_value if success else 'N/A',
                    'Status': 'Success' if success else 'Failed'
                })
            
            # Create DataFrame and save results
            df = pd.DataFrame(results)
            df = df.sort_values('Mutation').reset_index(drop=True)
            
            # Save as TSV
            df.to_csv(results_file, sep='\t', index=False)
            
            # Calculate summary statistics
            successful_results = df[df['Status'] == 'Success']
            if not successful_results.empty:
                ddg_values = pd.to_numeric(successful_results['ddG_REU'], errors='coerce')
                ddg_values = ddg_values.dropna()
                
                if not ddg_values.empty:
                    self.logger.info(f"Summary Statistics:")
                    self.logger.info(f"  Total mutations: {len(df)}")
                    self.logger.info(f"  Successful calculations: {len(successful_results)}")
                    self.logger.info(f"  Mean ΔΔG: {ddg_values.mean():.3f} REU")
                    self.logger.info(f"  Std ΔΔG: {ddg_values.std():.3f} REU")
                    self.logger.info(f"  Min ΔΔG: {ddg_values.min():.3f} REU")
                    self.logger.info(f"  Max ΔΔG: {ddg_values.max():.3f} REU")
                    
                    # Count stabilizing vs destabilizing
                    stabilizing = (ddg_values < 0).sum()
                    destabilizing = (ddg_values > 0).sum()
                    neutral = (ddg_values == 0).sum()
                    
                    self.logger.info(f"  Stabilizing (ΔΔG < 0): {stabilizing}")
                    self.logger.info(f"  Destabilizing (ΔΔG > 0): {destabilizing}")
                    self.logger.info(f"  Neutral (ΔΔG = 0): {neutral}")
            
            self._mark_step_complete("aggregate_ddg_results")
            self.logger.info(f"ΔΔG results aggregated and saved to: {results_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"ΔΔG result aggregation failed: {e}")
            return False
    
    def generate_visualizations(self) -> bool:
        """Generate structure visualizations using PyMOL."""
        if self._is_step_complete("generate_visualizations"):
            self.logger.info("Visualization generation already completed, skipping...")
            return True
        
        if self.config["skip_visualization"]:
            self.logger.info("Visualization skipped by configuration")
            self._mark_step_complete("generate_visualizations")
            return True
        
        self.logger.info("Generating structure visualizations...")
        
        relaxed_pdb = os.path.join(self.config.relax_dir, "relaxed_wildtype.pdb")
        if not os.path.exists(relaxed_pdb):
            self.logger.error("Relaxed wildtype structure not found")
            return False
        
        # Find mutant PDB files
        mut_pdbs = [f for f in os.listdir(self.config.ddg_dir) 
                   if f.endswith('.pdb') and f.startswith('mut_')]
        
        if not mut_pdbs:
            self.logger.warning("No mutant PDB files found for visualization")
            self._mark_step_complete("generate_visualizations")
            return True
        
        try:
            # Try to import PyMOL
            import pymol
            from pymol import cmd
            
            # Initialize PyMOL
            pymol.pymol_argv = ['pymol', '-qc']  # Quiet and no GUI
            pymol.finish_launching()
            
            vis_config = self.config.get("visualization", {})
            
            # Set up environment
            cmd.set('ray_opaque_background', 'off')
            cmd.set('ray_shadows', 'on')
            cmd.set('antialias', 2)
            cmd.set('orthoscopic', 'on')
            
            # Load wild-type structure
            cmd.load(relaxed_pdb, 'wildtype')
            cmd.show('cartoon', 'wildtype')
            cmd.color('cyan', 'wildtype')
            cmd.set('cartoon_transparency', 0.3, 'wildtype')
            
            for i, mut_pdb in enumerate(sorted(mut_pdbs)):
                mut_path = os.path.join(self.config.ddg_dir, mut_pdb)
                mut_name = f"mutant_{i+1}"
                mut_file = mut_pdb.replace('.pdb', '')
                
                self.logger.debug(f"Processing {mut_file}...")
                
                # Load mutant structure
                cmd.load(mut_path, mut_name)
                cmd.show('cartoon', mut_name)
                cmd.color('orange', mut_name)
                
                # Align to wild-type
                alignment_result = cmd.align(mut_name, 'wildtype')
                rmsd = alignment_result[0] if alignment_result else 0.0
                
                # Set up view
                cmd.bg_color('white')
                cmd.zoom()
                cmd.orient()
                
                # Save image
                output_file = os.path.join(self.config.images_dir, f'{mut_file}_superposition.png')
                cmd.png(output_file, 
                       width=vis_config.get('image_width', 800),
                       height=vis_config.get('image_height', 600),
                       dpi=vis_config.get('image_dpi', 300),
                       ray=1 if vis_config.get('ray_tracing', True) else 0)
                
                self.logger.debug(f"Created: {output_file} (RMSD: {rmsd:.3f} Å)")
                
                # Remove mutant for next iteration
                cmd.delete(mut_name)
            
            # Clean up PyMOL
            cmd.quit()
            
            self._mark_step_complete("generate_visualizations")
            self.logger.info(f"Generated {len(mut_pdbs)} structure superposition images")
            return True
            
        except ImportError:
            self.logger.warning("PyMOL not available - skipping visualization")
            self._mark_step_complete("generate_visualizations")
            return True
        except Exception as e:
            self.logger.error(f"Visualization generation failed: {e}")
            return False
    
    def create_summary_report(self) -> bool:
        """Create a comprehensive summary report of the pipeline results."""
        if self._is_step_complete("create_summary_report"):
            self.logger.info("Summary report already created, skipping...")
            return True
        
        self.logger.info("Creating summary report...")
        
        summary_file = os.path.join(self.config["output_dir"], "pipeline_summary.txt")
        relaxed_pdb = os.path.join(self.config.relax_dir, "relaxed_wildtype.pdb")
        ddg_results_file = os.path.join(self.config.ddg_dir, "ddg_results.txt")
        
        try:
            with open(summary_file, 'w') as f:
                f.write("Rosetta Protein Folding and Mutation Analysis Pipeline Summary\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"Protein: {self.protein_name}\n")
                f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write("Input Files:\n")
                f.write(f"- Protein FASTA: {self.config['protein_fasta']}\n")
                f.write(f"- Mutations File: {self.config['mutations_file']}\n\n")
                
                f.write("Pipeline Configuration:\n")
                f.write(f"- Rosetta Path: {self.config['rosetta_path']}\n")
                f.write(f"- Folding Structures: {self.config['nstruct_folding']}\n")
                f.write(f"- Relaxation Structures: {self.config['nstruct_relax']}\n")
                f.write(f"- DDG Iterations: {self.config['ddg_iterations']}\n\n")
                
                f.write("Pipeline Steps Completed:\n")
                completed_steps = self.state.get("completed_steps", [])
                step_descriptions = {
                    "setup_directories": "Directory structure setup",
                    "generate_fragments": "Fragment generation",
                    "ab_initio_folding": f"Ab initio folding ({self.config['nstruct_folding']} structures)",
                    "select_best_model": "Best model selection",
                    "relax_structure": f"Structure relaxation ({self.config['nstruct_relax']} structures)",
                    "parse_mutations": "Mutation parsing",
                    "calculate_ddg": f"ΔΔG calculations ({self.config['ddg_iterations']} iterations each)",
                    "aggregate_ddg_results": "Result aggregation",
                    "generate_visualizations": "Structure visualization",
                    "create_summary_report": "Summary report generation"
                }
                
                for i, (step, description) in enumerate(step_descriptions.items(), 1):
                    status = "✓" if step in completed_steps else "✗"
                    f.write(f"{i}. {status} {description}\n")
                f.write("\n")
                
                f.write("Output Files:\n")
                f.write(f"- Wild-type relaxed structure: {relaxed_pdb}\n")
                f.write(f"- ΔΔG results table: {ddg_results_file}\n")
                f.write(f"- Structure images: {self.config.images_dir}/*.png\n")
                f.write(f"- Detailed logs: */log files in each subdirectory\n\n")
                
                f.write("Results Summary:\n")
                f.write("=" * 15 + "\n")
                
                # Add ΔΔG results if available
                if os.path.exists(ddg_results_file):
                    f.write("\n")
                    with open(ddg_results_file, 'r') as ddg_f:
                        f.write(ddg_f.read())
                else:
                    f.write("ΔΔG results not yet available.\n")
            
            self._mark_step_complete("create_summary_report")
            self.logger.info(f"Summary report created: {summary_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"Summary report creation failed: {e}")
            return False
    
    def run_pipeline(self, resume: bool = False) -> bool:
        """Run the complete pipeline with optional resume capability."""
        if not resume:
            self.logger.info("Starting Rosetta Pipeline from the beginning")
            # Reset state for fresh run
            self.state["completed_steps"] = []
            self._save_state()
        else:
            completed_steps = len(self.state.get("completed_steps", []))
            self.logger.info(f"Resuming Rosetta Pipeline ({completed_steps} steps already completed)")
        
        pipeline_steps = [
            ("validate_dependencies", "Validating dependencies", self.validate_dependencies),
            ("setup_directories", "Setting up directories", self.setup_directories),
            ("generate_fragments", "Generating fragments", self.generate_fragments),
            ("ab_initio_folding", "Running ab initio folding", self.run_ab_initio_folding),
            ("select_best_model", "Selecting best model", self.select_best_model),
            ("relax_structure", "Relaxing structure", self.relax_structure),
            ("parse_mutations", "Parsing mutations", self.parse_mutations),
            ("calculate_ddg", "Calculating ΔΔG", self.calculate_ddg),
            ("aggregate_ddg_results", "Aggregating results", self.aggregate_ddg_results),
            ("generate_visualizations", "Generating visualizations", self.generate_visualizations),
            ("create_summary_report", "Creating summary report", self.create_summary_report)
        ]
        
        start_time = datetime.now()
        
        for step_name, step_desc, step_func in pipeline_steps:
            self.logger.info(f"Pipeline Step: {step_desc}")
            
            try:
                if not step_func():
                    self.logger.error(f"Pipeline failed at step: {step_desc}")
                    return False
            except KeyboardInterrupt:
                self.logger.info("Pipeline interrupted by user")
                return False
            except Exception as e:
                self.logger.error(f"Unexpected error in step '{step_desc}': {e}")
                return False
        
        end_time = datetime.now()
        duration = end_time - start_time
        
        self.logger.info("=" * 50)
        self.logger.info("Pipeline completed successfully!")
        self.logger.info(f"Total execution time: {duration}")
        self.logger.info(f"Results available in: {self.config['output_dir']}")
        self.logger.info("=" * 50)
        
        return True

def main():
    """Main entry point for the pipeline."""
    parser = argparse.ArgumentParser(
        description="Self-contained Rosetta protein folding and mutation analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline with default settings
  python rosetta_pipeline.py
  
  # Run with custom config file
  python rosetta_pipeline.py --config my_config.yaml
  
  # Resume interrupted pipeline
  python rosetta_pipeline.py --resume
  
  # Override specific parameters
  python rosetta_pipeline.py --nstruct-folding 5000 --ddg-iterations 100
  
  # Run with custom input files
  python rosetta_pipeline.py --protein-fasta my_protein.faa --mutations-file my_mutations.txt
  
  # Skip visualization
  python rosetta_pipeline.py --skip-visualization
  
  # Dry run (validate only)
  python rosetta_pipeline.py --validate-only

Configuration:
  The pipeline can be configured via:
  1. Built-in defaults
  2. Configuration file (--config)  
  3. Command line arguments (highest priority)
  
  Configuration files can be YAML or JSON format.
        """
    )
    
    # Input/Output options
    parser.add_argument('--config', '-c', 
                       help='Configuration file (YAML or JSON)')
    parser.add_argument('--protein-fasta',
                       help='Protein FASTA file (default: protein.faa)')
    parser.add_argument('--mutations-file', 
                       help='Mutations file (default: mutations.txt)')
    parser.add_argument('--output-dir', '-o',
                       help='Output directory (default: rosetta_results)')
    
    # Rosetta configuration
    parser.add_argument('--rosetta-path',
                       help='Path to Rosetta installation')
    
    # Pipeline parameters
    parser.add_argument('--nstruct-folding', type=int,
                       help='Number of folding structures (default: 1000)')
    parser.add_argument('--nstruct-relax', type=int,
                       help='Number of relaxation structures (default: 10)')
    parser.add_argument('--ddg-iterations', type=int,
                       help='Number of DDG iterations (default: 50)')
    
    # Pipeline control
    parser.add_argument('--resume', action='store_true',
                       help='Resume interrupted pipeline')
    parser.add_argument('--validate-only', action='store_true',
                       help='Only validate dependencies and inputs')
    parser.add_argument('--skip-visualization', action='store_true',
                       help='Skip PyMOL visualization step')
    
    # Logging
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Suppress non-error output')
    
    args = parser.parse_args()
    
    # Set up logging level
    if args.quiet:
        logging.getLogger('RosettaPipeline').setLevel(logging.ERROR)
    elif args.verbose:
        logging.getLogger('RosettaPipeline').setLevel(logging.DEBUG)
    
    # Prepare CLI overrides
    cli_overrides = {}
    if args.protein_fasta:
        cli_overrides['protein_fasta'] = args.protein_fasta
    if args.mutations_file:
        cli_overrides['mutations_file'] = args.mutations_file
    if args.output_dir:
        cli_overrides['output_dir'] = args.output_dir
    if args.rosetta_path:
        cli_overrides['rosetta_path'] = args.rosetta_path
    if args.nstruct_folding:
        cli_overrides['nstruct_folding'] = args.nstruct_folding
    if args.nstruct_relax:
        cli_overrides['nstruct_relax'] = args.nstruct_relax
    if args.ddg_iterations:
        cli_overrides['ddg_iterations'] = args.ddg_iterations
    if args.skip_visualization:
        cli_overrides['skip_visualization'] = True
    
    try:
        # Initialize configuration
        config = RosettaPipelineConfig(config_file=args.config, cli_overrides=cli_overrides)
        
        # Initialize pipeline
        pipeline = RosettaPipeline(config)
        
        # Validate dependencies
        if not pipeline.validate_dependencies():
            print("Dependency validation failed. Please check your Rosetta installation and input files.")
            sys.exit(1)
        
        if args.validate_only:
            print("Validation completed successfully.")
            sys.exit(0)
        
        # Run pipeline
        success = pipeline.run_pipeline(resume=args.resume)
        
        if success:
            print(f"Pipeline completed successfully! Results in: {config['output_dir']}")
            sys.exit(0)
        else:
            print("Pipeline failed. Check logs for details.")
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\nPipeline interrupted by user.")
        sys.exit(130)
    except Exception as e:
        print(f"Pipeline failed with error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()