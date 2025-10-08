#!/usr/bin/env python

import sys
import csv

def process_vcf_line(vcf_line_str, sample_id, sample_name, dataset, vcf_type):
    """
    Parses a VCF line string and prints the formatted output.
    """
    if not vcf_line_str or vcf_line_str.isspace():
        return

    # The embedded VCF line uses the literal string '\t' as a separator.
    fields = vcf_line_str.strip().split('\\t')
    
    # A valid VCF record from a gVCF should have at least 10 fields
    # if it has genotype information.
    if len(fields) < 10:
        return

    format_field = fields[8]
    sample_field = fields[9]

    format_keys = format_field.split(':')
    format_values = sample_field.split(':')

    for key, value in zip(format_keys, format_values):
        # Construct the output string in TSV format
        output_line = f"{sample_id}\t{sample_name}\t{dataset}\t{vcf_type}\t{key}\t{value}"
        print(output_line)

def main(input_file):
    """
    Processes the input file and prints the transformed data to standard output.
    """
    # Print the header in TSV format
    print("sample_id\tsample_name\tdataset\tvcf\tformat\tinfo")

    with open(input_file, 'r') as f:
        # Use csv.reader to handle tab-delimited format properly
        reader = csv.reader(f, delimiter='\t')
        
        # Skip the header line
        try:
            next(reader)
        except StopIteration:
            # Handle empty file
            return

        for row in reader:
            if not row:
                continue

            # As per the header: sample_id, sample_name, gvcf_path, reblocked_gvcf, 
            # gvcf_line, reblocked_gvcf_line, dataset
            if len(row) < 7:
                continue

            sample_id = row[0]
            sample_name = row[1]
            gvcf_line = row[4]
            reblocked_gvcf_line = row[5]
            dataset = row[6]

            # Process the 'gvcf_line' as "raw"
            process_vcf_line(gvcf_line, sample_id, sample_name, dataset, "raw")

            # Process the 'reblocked_gvcf_line' as "reblocked"
            process_vcf_line(reblocked_gvcf_line, sample_id, sample_name, dataset, "reblocked")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python process_vcf_data.py <input_file>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])
