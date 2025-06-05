#!/usr/bin/env python3
"""
Script to generate bcftools invocations from VIDs (Variant IDs).

A VID consists of four components separated by dashes:
chromosome-position-reference_allele-variant_allele

Example usage:
python generate_bcftools_commands.py input_vids.txt
"""

import sys
import argparse


def parse_vid(vid):
    """Parse a VID into its components."""
    parts = vid.strip().split('-')
    if len(parts) != 4:
        raise ValueError(f"Invalid VID format: {vid}. Expected 4 components separated by dashes.")
    
    chromosome, position, ref_allele, var_allele = parts
    position = int(position)
    
    return chromosome, position, ref_allele, var_allele


def calculate_insert_size(ref_allele, var_allele):
    """Calculate insert size as length of variant allele minus length of reference allele."""
    return len(var_allele) - len(ref_allele)


def generate_bcftools_command(vid):
    """Generate bcftools command for a given VID."""
    chromosome, position, ref_allele, var_allele = parse_vid(vid)
    
    insert_size = calculate_insert_size(ref_allele, var_allele)
    abs_size = abs(insert_size)
    
    start = position - abs_size
    end = position + abs_size
    
    # Ensure start is at least 1 (genomic coordinates are 1-based)
    start = max(1, start)
    
    command = f"bcftools view --no-header -i '(ILEN = {insert_size})' --regions chr{chromosome}:{start}-{end} sites-only.vcf"
    
    return command


def main():
    parser = argparse.ArgumentParser(description='Generate bcftools commands from VIDs')
    parser.add_argument('input_file', help='File containing one VID per line')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    
    args = parser.parse_args()
    
    output_file = open(args.output, 'w') if args.output else sys.stdout
    
    try:
        with open(args.input_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):  # Skip empty lines and comments
                    continue
                
                try:
                    command = generate_bcftools_command(line)
                    print(command, file=output_file)
                except ValueError as e:
                    print(f"Error on line {line_num}: {e}", file=sys.stderr)
                    continue
    
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    finally:
        if args.output:
            output_file.close()


if __name__ == '__main__':
    main()