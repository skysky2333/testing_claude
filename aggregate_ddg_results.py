#!/usr/bin/env python3
"""
Helper script to aggregate ΔΔG results from multiple ddg_monomer runs
Author: Claude Code Assistant
"""

import sys
import os
import re
import argparse
import pandas as pd
from pathlib import Path
import numpy as np

def extract_ddg_from_log(log_file):
    """
    Extract ΔΔG value from ddg_monomer log file
    Returns: (ddg_value, success)
    """
    if not os.path.exists(log_file):
        return None, False
    
    ddg_values = []
    
    try:
        with open(log_file, 'r') as f:
            for line in f:
                # Look for ddG lines in the output
                if 'ddG' in line and ':' in line:
                    # Try to extract numerical value
                    # Format can be like "ddG: -1.234" or similar
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
        print(f"Error reading {log_file}: {e}", file=sys.stderr)
        return None, False

def parse_rosetta_scorefile(score_file):
    """
    Parse Rosetta score file to extract ddG values
    Returns: list of (tag, ddg_value) tuples
    """
    results = []
    
    if not os.path.exists(score_file):
        return results
    
    try:
        with open(score_file, 'r') as f:
            header_line = None
            for line in f:
                line = line.strip()
                if line.startswith('SCORE:') and 'description' in line:
                    # This is the header line
                    header_line = line.split()[1:]  # Remove 'SCORE:'
                elif line.startswith('SCORE:') and header_line:
                    # This is a data line
                    values = line.split()[1:]  # Remove 'SCORE:'
                    if len(values) == len(header_line):
                        score_dict = dict(zip(header_line, values))
                        
                        # Look for ddG-related columns
                        ddg_value = None
                        for col in ['ddG', 'total_score', 'score']:
                            if col in score_dict:
                                try:
                                    ddg_value = float(score_dict[col])
                                    break
                                except ValueError:
                                    continue
                        
                        tag = score_dict.get('description', 'unknown')
                        results.append((tag, ddg_value))
        
    except Exception as e:
        print(f"Error reading score file {score_file}: {e}", file=sys.stderr)
    
    return results

def aggregate_ddg_results(ddg_dir, output_file=None, format='tsv'):
    """
    Aggregate ΔΔG results from a directory containing ddg_monomer outputs
    """
    ddg_dir = Path(ddg_dir)
    
    if not ddg_dir.exists():
        raise FileNotFoundError(f"DDG directory not found: {ddg_dir}")
    
    results = []
    
    # Find all mutation label files
    label_files = list(ddg_dir.glob("mutation_*_label.txt"))
    
    if not label_files:
        print("No mutation label files found. Looking for log files...", file=sys.stderr)
        # Fallback: look for log files directly
        log_files = list(ddg_dir.glob("ddg_*.log"))
        for log_file in log_files:
            ddg_value, success = extract_ddg_from_log(log_file)
            mutation_id = log_file.stem.replace('ddg_', '')
            results.append({
                'Mutation': f'mutation_{mutation_id}',
                'ddG_REU': ddg_value if success else 'N/A',
                'Status': 'Success' if success else 'Failed'
            })
    else:
        # Process each mutation
        for label_file in sorted(label_files):
            # Extract mutation number
            match = re.search(r'mutation_(\d+)_label\.txt', label_file.name)
            if not match:
                continue
            
            mut_num = match.group(1)
            
            # Read mutation label
            try:
                with open(label_file, 'r') as f:
                    mutation_label = f.read().strip()
            except Exception:
                mutation_label = f'mutation_{mut_num}'
            
            # Find corresponding log file
            log_file = ddg_dir / f"ddg_{mut_num}.log"
            ddg_value, success = extract_ddg_from_log(log_file)
            
            # Try to find score file as alternative
            if not success:
                score_file = ddg_dir / f"score_{mut_num}.sc"
                if score_file.exists():
                    score_results = parse_rosetta_scorefile(score_file)
                    if score_results:
                        ddg_value = score_results[0][1]  # Take first result
                        success = ddg_value is not None
            
            results.append({
                'Mutation': mutation_label,
                'ddG_REU': ddg_value if success else 'N/A',
                'Status': 'Success' if success else 'Failed'
            })
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Sort by mutation if possible
    df = df.sort_values('Mutation').reset_index(drop=True)
    
    # Output results
    if output_file:
        if format.lower() == 'csv':
            df.to_csv(output_file, index=False)
        elif format.lower() == 'tsv':
            df.to_csv(output_file, sep='\t', index=False)
        elif format.lower() == 'excel':
            df.to_excel(output_file, index=False)
        else:
            raise ValueError(f"Unsupported format: {format}")
        
        print(f"Results written to: {output_file}")
    else:
        # Print to stdout
        if format.lower() == 'csv':
            print(df.to_csv(index=False))
        else:
            print(df.to_csv(sep='\t', index=False))
    
    # Summary statistics
    successful_results = df[df['Status'] == 'Success']
    if not successful_results.empty:
        ddg_values = pd.to_numeric(successful_results['ddG_REU'], errors='coerce')
        ddg_values = ddg_values.dropna()
        
        if not ddg_values.empty:
            print(f"\nSummary Statistics:", file=sys.stderr)
            print(f"  Total mutations: {len(df)}", file=sys.stderr)
            print(f"  Successful calculations: {len(successful_results)}", file=sys.stderr)
            print(f"  Mean ΔΔG: {ddg_values.mean():.3f} REU", file=sys.stderr)
            print(f"  Std ΔΔG: {ddg_values.std():.3f} REU", file=sys.stderr)
            print(f"  Min ΔΔG: {ddg_values.min():.3f} REU", file=sys.stderr)
            print(f"  Max ΔΔG: {ddg_values.max():.3f} REU", file=sys.stderr)
            
            # Count stabilizing vs destabilizing
            stabilizing = (ddg_values < 0).sum()
            destabilizing = (ddg_values > 0).sum()
            neutral = (ddg_values == 0).sum()
            
            print(f"  Stabilizing (ΔΔG < 0): {stabilizing}", file=sys.stderr)
            print(f"  Destabilizing (ΔΔG > 0): {destabilizing}", file=sys.stderr)
            print(f"  Neutral (ΔΔG = 0): {neutral}", file=sys.stderr)
    
    return df

def main():
    parser = argparse.ArgumentParser(
        description="Aggregate ΔΔG results from ddg_monomer calculations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python aggregate_ddg_results.py ddg_output/
  python aggregate_ddg_results.py ddg_output/ -o results.tsv
  python aggregate_ddg_results.py ddg_output/ -o results.csv -f csv
        """
    )
    
    parser.add_argument('ddg_dir',
                       help='Directory containing ddg_monomer output files')
    parser.add_argument('-o', '--output',
                       help='Output file (default: stdout)')
    parser.add_argument('-f', '--format', 
                       choices=['tsv', 'csv', 'excel'],
                       default='tsv',
                       help='Output format (default: tsv)')
    parser.add_argument('--summary-only', action='store_true',
                       help='Only print summary statistics')
    
    args = parser.parse_args()
    
    try:
        df = aggregate_ddg_results(args.ddg_dir, args.output, args.format)
        
        if args.summary_only:
            # Summary was already printed to stderr
            pass
        elif not args.output:
            # Results were already printed to stdout
            pass
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()