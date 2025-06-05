#!/usr/bin/env python3
"""
Script to filter VCF lines based on a set of VIDs (Variant IDs).

Reads VIDs from the first file into a set, then filters the VCF in the second file
to only output header lines and variant lines that match VIDs in the set.

Example usage:
python filter_vcf_by_vids.py vids.txt variants.vcf
"""

import sys
import argparse


def load_vids(filename):
    """Load VIDs from file into a set."""
    vids = set()
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):  # Skip empty lines and comments
                vids.add(line)
    return vids


def construct_vid_from_vcf_line(chrom, pos, ref, alt):
    """Construct VID string from VCF line components."""
    # Remove 'chr' prefix if present
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    
    return f"{chrom}-{pos}-{ref}-{alt}"


def process_vcf(vcf_filename, vids):
    """Process VCF file and output matching lines."""
    with open(vcf_filename, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            
            # Output all header lines
            if line.startswith('#'):
                print(line)
                continue
            
            # Process variant lines
            fields = line.split('\t')
            if len(fields) >= 5:
                chrom = fields[0]
                pos = fields[1]
                ref = fields[3]
                alt = fields[4]
                
                vid = construct_vid_from_vcf_line(chrom, pos, ref, alt)
                
                if vid in vids:
                    print(vid)


def main():
    parser = argparse.ArgumentParser(description='Filter VCF lines based on VIDs')
    parser.add_argument('vids_file', help='File containing VIDs, one per line')
    parser.add_argument('vcf_file', help='VCF file to filter')
    
    args = parser.parse_args()
    
    try:
        # Load VIDs into memory
        vids = load_vids(args.vids_file)
        
        # Process VCF file
        process_vcf(args.vcf_file, vids)
        
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()