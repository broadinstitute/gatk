#!/usr/bin/env python3
"""
Script to extract QUAL and GQ values from VCF lines in TSV format.
Takes a TSV with gvcf_line and reblocked_gvcf_line columns and outputs
quality metrics for comparison.
"""

import sys
import csv
import argparse


def parse_vcf_line_for_qual(vcf_line):
    """
    Extract QUAL value from a VCF line.
    
    Args:
        vcf_line: Tab-delimited VCF line
    
    Returns:
        String containing the QUAL value (6th column, 0-based index 5)
    """
    if not vcf_line or vcf_line.strip() == "":
        return "."
    
    # Handle potential encoding issues or extra whitespace
    vcf_line = vcf_line.strip()
    fields = vcf_line.split('\\t')
    if len(fields) < 6:  # Need at least 6 fields to get QUAL (0-based index 5)
        return "."
    
    return fields[5]  # QUAL is at index 5 (0-based)


def parse_vcf_line_for_gq(vcf_line):
    """
    Extract GQ value from a VCF line.
    
    Args:
        vcf_line: Tab-delimited VCF line
    
    Returns:
        String containing the GQ value
    """
    if not vcf_line or vcf_line.strip() == "":
        return "."
    
    # Handle potential encoding issues or extra whitespace
    vcf_line = vcf_line.strip()
    fields = vcf_line.split('\\t')
    if len(fields) < 10:  # Need at least 10 fields for FORMAT and sample data
        return "."
    
    format_field = fields[8]  # FORMAT is at index 8
    sample_field = fields[9]  # First sample data is at index 9
    
    # Parse FORMAT field to find GQ index
    format_elements = format_field.split(':')
    try:
        gq_index = format_elements.index('GQ')
    except ValueError:
        # GQ not found in FORMAT field
        return "."
    
    # Parse sample field to extract GQ value
    sample_elements = sample_field.split(':')
    if gq_index >= len(sample_elements):
        return "."
    
    return sample_elements[gq_index]


def process_row(row):
    """
    Process a single row from the input TSV.
    
    Args:
        row: Dictionary containing row data
    
    Returns:
        Dictionary with extracted quality values
    """
    # Extract basic fields
    result = {
        'sample_id': row.get('sample_id', ''),
        'chr': row.get('chr', ''),
        'input_position': row.get('input_position', ''),
        'input_ref': row.get('input_ref', ''),
        'input_alt': row.get('input_alt', '')
    }
    
    # Extract QUAL and GQ from gvcf_line
    gvcf_line = row.get('gvcf_line', '')
    result['gvcf_qual'] = parse_vcf_line_for_qual(gvcf_line)
    result['gvcf_gq'] = parse_vcf_line_for_gq(gvcf_line)
    
    # Extract QUAL and GQ from reblocked_gvcf_line
    reblocked_line = row.get('reblocked_gvcf_line', '')
    result['reblocked_qual'] = parse_vcf_line_for_qual(reblocked_line)
    result['reblocked_gq'] = parse_vcf_line_for_gq(reblocked_line)
    
    return result


def main():
    parser = argparse.ArgumentParser(description='Extract VCF quality metrics from TSV')
    parser.add_argument('input_file', nargs='?', default='-', 
                        help='Input TSV file (default: stdin)')
    
    args = parser.parse_args()
    
    # Open input file
    if args.input_file == '-':
        input_file = sys.stdin
    else:
        try:
            input_file = open(args.input_file, 'r')
        except FileNotFoundError:
            print(f"Error: Input file '{args.input_file}' not found.", file=sys.stderr)
            sys.exit(1)
    
    try:
        # Read input TSV
        reader = csv.DictReader(input_file, delimiter='\t')
        
        # Write output header
        output_headers = ['sample_id', 'chr', 'input_position', 'input_ref', 'input_alt',
                         'gvcf_qual', 'gvcf_gq', 'reblocked_qual', 'reblocked_gq']
        writer = csv.DictWriter(sys.stdout, fieldnames=output_headers, delimiter='\t')
        writer.writeheader()
        
        # Process each row
        for row in reader:
            result = process_row(row)
            writer.writerow(result)
            
    finally:
        if args.input_file != '-':
            input_file.close()


if __name__ == '__main__':
    main()