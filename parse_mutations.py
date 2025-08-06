#!/usr/bin/env python3
"""
Helper script to parse mutations file and create individual mutation files for Rosetta ddg_monomer
Author: Claude Code Assistant
"""

import sys
import os
import argparse
from pathlib import Path

def parse_mutation_string(mutation_str):
    """
    Parse a single mutation string like 'A1G' into components
    Returns: (wild_type_aa, position, mutant_aa)
    """
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

def create_ddg_mutation_file(mutations, output_file):
    """
    Create a mutation file in ddg_monomer format
    """
    num_mutations = len(mutations)
    
    with open(output_file, 'w') as f:
        f.write(f"{num_mutations}\n")
        for wt_aa, pos, mut_aa in mutations:
            f.write(f"{wt_aa} {pos} {mut_aa}\n")

def parse_mutations_file(input_file, output_dir):
    """
    Parse the mutations file and create individual ddg_monomer input files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    mutation_files = []
    labels = []
    
    with open(input_file, 'r') as f:
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
                    wt_aa, pos, mut_aa = parse_mutation_string(mut_str)
                    mutations.append((wt_aa, pos, mut_aa))
                
                # Create output files
                mut_file = output_dir / f"mutation_{line_num}.txt"
                label_file = output_dir / f"mutation_{line_num}_label.txt"
                
                create_ddg_mutation_file(mutations, mut_file)
                
                with open(label_file, 'w') as lf:
                    lf.write(line)
                
                mutation_files.append(str(mut_file))
                labels.append(line)
                
                print(f"Created mutation file {line_num}: {line}")
                
            except ValueError as e:
                print(f"Error processing line {line_num} '{line}': {e}", file=sys.stderr)
                continue
    
    return mutation_files, labels

def validate_sequence_file(fasta_file):
    """
    Validate that the FASTA file exists and contains a protein sequence
    """
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    
    if not lines:
        raise ValueError("FASTA file is empty")
    
    if not lines[0].startswith('>'):
        raise ValueError("FASTA file does not start with header line")
    
    sequence = ''.join(line.strip() for line in lines[1:] if not line.startswith('>'))
    
    if not sequence:
        raise ValueError("No sequence found in FASTA file")
    
    # Check for valid amino acid codes
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    invalid_chars = set(sequence.upper()) - valid_aa
    if invalid_chars:
        print(f"Warning: Non-standard amino acids found: {invalid_chars}")
    
    return len(sequence)

def main():
    parser = argparse.ArgumentParser(
        description="Parse mutations file for Rosetta ddg_monomer analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python parse_mutations.py mutations.txt ddg_input/
  python parse_mutations.py --validate-fasta protein.faa mutations.txt ddg_input/
  
Mutation file format:
  Each line contains one or more mutations separated by commas.
  Mutation format: [WT_AA][Position][Mutant_AA]
  
  Example mutations.txt:
    A1G
    S2M  
    A1G,S2M
        """
    )
    
    parser.add_argument('mutations_file', 
                       help='Input file containing mutations')
    parser.add_argument('output_dir',
                       help='Output directory for ddg_monomer input files')
    parser.add_argument('--validate-fasta', 
                       help='FASTA file to validate mutations against')
    parser.add_argument('--dry-run', action='store_true',
                       help='Parse and validate without creating files')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.mutations_file):
        print(f"Error: Mutations file not found: {args.mutations_file}", file=sys.stderr)
        sys.exit(1)
    
    sequence_length = None
    if args.validate_fasta:
        try:
            sequence_length = validate_sequence_file(args.validate_fasta)
            print(f"Validated FASTA file: {sequence_length} amino acids")
        except (FileNotFoundError, ValueError) as e:
            print(f"Error validating FASTA file: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Parse mutations
    try:
        if args.dry_run:
            print("Dry run mode - parsing mutations without creating files")
            
        mutation_files, labels = parse_mutations_file(args.mutations_file, args.output_dir)
        
        print(f"\nSummary:")
        print(f"  Processed {len(labels)} mutation sets")
        print(f"  Output directory: {args.output_dir}")
        
        if not args.dry_run:
            print(f"  Created {len(mutation_files)} mutation files")
        
        # Validate mutations against sequence if provided
        if sequence_length:
            for i, label in enumerate(labels):
                mutations = [m.strip() for m in label.split(',')]
                for mut in mutations:
                    try:
                        wt_aa, pos, mut_aa = parse_mutation_string(mut)
                        pos_int = int(pos)
                        if pos_int > sequence_length:
                            print(f"Warning: Mutation {mut} position {pos_int} exceeds sequence length {sequence_length}")
                    except ValueError:
                        continue
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()