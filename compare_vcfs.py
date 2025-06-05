#!/usr/bin/env python3
"""
Script to compare input variants VCF with left-aligned variants VCF.
Outputs tab-delimited comparison data to standard output.
"""

import sys
import argparse


def parse_vcf_line(line):
    """Parse a VCF line and return relevant fields."""
    if line.startswith('#'):
        return None
    
    fields = line.strip().split('\t')
    if len(fields) < 8:
        return None
    
    return {
        'CHROM': fields[0],
        'POS': fields[1],
        'REF': fields[3],
        'ALT': fields[4],
        'INFO': fields[7]
    }


def create_key(chrom, ref, alt, info):
    """Create dictionary key from CHROM, REF length, ALT length, and INFO."""
    return f"{chrom}:{len(ref)}:{len(alt)}:{info}"


def create_vid(chrom, pos, ref, alt):
    """Create variant ID in the format CHROM-POS-REF-ALT (without chr prefix)."""
    # Remove 'chr' prefix if present
    chrom_clean = chrom.replace('chr', '') if chrom.startswith('chr') else chrom
    return f"{chrom_clean}-{pos}-{ref}-{alt}"


def load_vcf_to_dict(vcf_file):
    """Load VCF file into dictionary with composite keys."""
    vcf_dict = {}
    
    with open(vcf_file, 'r') as f:
        for line in f:
            variant = parse_vcf_line(line)
            if variant is None:
                continue
            
            key = create_key(variant['CHROM'], variant['REF'], variant['ALT'], variant['INFO'])
            vcf_dict[key] = variant
    
    return vcf_dict


def main():
    parser = argparse.ArgumentParser(description='Compare input variants VCF with left-aligned variants VCF')
    parser.add_argument('input_vcf', help='Input variants VCF file')
    parser.add_argument('left_aligned_vcf', help='Left-aligned variants VCF file')
    
    args = parser.parse_args()
    
    # Load both VCF files into dictionaries
    input_variants = load_vcf_to_dict(args.input_vcf)
    left_aligned_variants = load_vcf_to_dict(args.left_aligned_vcf)
    
    # Print header
    print("vid\tinput_position\tinput_ref\tinput_alt\tleft_aligned_position\tleft_aligned_ref\tleft_aligned_alt\tinfo_field")
    
    # Iterate through left-aligned variants and find matches in input variants
    for key, left_variant in left_aligned_variants.items():
        if key not in input_variants:
            print(f"ERROR: No matching input variant found for key: {key}", file=sys.stderr)
            sys.exit(1)
        
        input_variant = input_variants[key]
        
        # Create VID using left-aligned variant data
        vid = create_vid(left_variant['CHROM'], left_variant['POS'], left_variant['REF'], left_variant['ALT'])
        
        # Output tab-delimited row
        print(f"{vid}\t{input_variant['POS']}\t{input_variant['REF']}\t{input_variant['ALT']}\t{left_variant['POS']}\t{left_variant['REF']}\t{left_variant['ALT']}\t{input_variant['INFO']}")


if __name__ == '__main__':
    main()