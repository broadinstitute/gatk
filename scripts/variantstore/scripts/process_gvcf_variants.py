#!/usr/bin/env python3
"""
Script to process variants by querying GVCF and reblocked GVCF files with bcftools.
Takes a JSON file with variant information and outputs enriched JSON with bcftools results.
"""

import json
import subprocess
import sys
import argparse


def run_bcftools(gvcf_path, chr_name, position, ref, alt):
    """
    Run bcftools view command and return the output line.
    
    Args:
        gvcf_path: Path to the GVCF file
        chr_name: Chromosome name (e.g., 'chr1')
        position: Position on chromosome
        ref: Reference allele
        alt: Alternate allele
    
    Returns:
        String containing the bcftools output line
    """
    cmd = [
        'bcftools', 'view', '--no-header',
        '--include', f'(ILEN={len(alt) - len(ref)})',
        '--regions', f'{chr_name}:{position}',
        gvcf_path
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Return the output, stripped of trailing whitespace
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running bcftools command: {' '.join(cmd)}", file=sys.stderr)
        print(f"Error message: {e.stderr}", file=sys.stderr)
        return f"ERROR: {e.stderr.strip()}"
    except FileNotFoundError:
        print("Error: bcftools not found. Please ensure bcftools is installed and in PATH.", file=sys.stderr)
        return "ERROR: bcftools not found"


def process_variant(variant_obj):
    """
    Process a single variant object by running bcftools on both GVCF files.
    
    Args:
        variant_obj: Dictionary containing variant information
    
    Returns:
        Dictionary with original fields plus bcftools results
    """
    # Extract required fields
    chr_name = variant_obj['chr']
    position = variant_obj['input_position']
    ref = variant_obj['input_ref']
    alt = variant_obj['input_alt']
    gvcf_path = variant_obj['gvcf_path']
    reblocked_gvcf = variant_obj['reblocked_gvcf']
    
    # Run bcftools on both files
    gvcf_line = run_bcftools(gvcf_path, chr_name, position, ref, alt)
    reblocked_gvcf_line = run_bcftools(reblocked_gvcf, chr_name, position, ref, alt)
    
    # Create result object with original fields plus new ones
    result = variant_obj.copy()
    result['gvcf_line'] = gvcf_line
    result['reblocked_gvcf_line'] = reblocked_gvcf_line
    
    return result


def main():
    parser = argparse.ArgumentParser(description='Process variants with bcftools queries')
    parser.add_argument('json_file', help='Input JSON file containing variant objects')
    parser.add_argument('--output', '-o', help='Output JSON file (default: stdout)')
    
    args = parser.parse_args()
    
    # Load input JSON
    try:
        with open(args.json_file, 'r') as f:
            variants = json.load(f)
    except FileNotFoundError:
        print(f"Error: Input file '{args.json_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in '{args.json_file}': {e}", file=sys.stderr)
        sys.exit(1)
    
    # Validate input is an array
    if not isinstance(variants, list):
        print("Error: Input JSON must be an array of objects.", file=sys.stderr)
        sys.exit(1)
    
    # Process each variant
    results = []
    for i, variant in enumerate(variants):
        print(f"Processing variant {i+1}/{len(variants)}...", file=sys.stderr)
        result = process_variant(variant)
        results.append(result)
        if i > 3:
            break
    
    # Output results
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Results written to {args.output}", file=sys.stderr)
    else:
        print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()