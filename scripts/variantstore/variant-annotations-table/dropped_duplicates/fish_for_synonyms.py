#!/usr/bin/env python

import sys

def process_vcf_line(line):
    """
    Parses a single VCF line and generates a bcftools command.
    """
    if line.startswith('#'):
        return

    fields = line.strip().split('\t')
    if len(fields) < 5:
        return

    chrom = fields[0]
    try:
        pos = int(fields[1])
    except ValueError:
        # Handle cases where position is not a valid integer
        return
        
    ref = fields[3]
    alt = fields[4]

    # Skip lines with no REF or ALT
    if not ref or not alt:
        return

    insert_length = len(alt) - len(ref)
    end_pos = pos + 200

    command = f"bcftools view --no-header --regions {chrom}:{pos}-{end_pos} --include '(ILEN={insert_length})' sites-only.vcf.gz"
    print(command)

def main():
    """
    Main function to read from stdin and process each line.
    """
    print("bcftools head sites-only.vcf.gz")
    for line in sys.stdin:
        process_vcf_line(line)

if __name__ == "__main__":
    main()
