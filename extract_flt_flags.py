#!/usr/bin/env python3
"""

Usage: python extract_flt_flags.py <input_file> [output_file]

Script to extract research_id and flt_ flag columns from a tab-delimited "soft-filtered" AoU manifest file, producing an
output TSV of research_id and flt_ column names.  For each row, if any flt_ column has a value of True, the script will
output a line like:

<research_id><tab><flt_column_name>

This output should be suitable for loading to a BigQuery table for further analysis. This script is intended for a
"soft-filtered" manifest, which will have flt_ columns that contain True values; this is not the same as the "hard-filtered"
manifest which is usually distributed via email and whose flt_ columns are all False.

Note this script will produce multiple research_id/flt_column rows if the same combination of research_id and flt_ column with a
True value exists multple times in the input, a frequent occurrence in the Foxtrot "soft-filtered" manifest (see occurrences of
the flt_is_rid_duped flag).

Loading the output of this script into BigQuery would look like:

bq load --project_id aou-genomics-curation-prod --source_format=CSV --skip_leading_rows=1 --field_delimiter="\t" \
    foxtrot.foxtrot_flt_flags flt_flags.tsv \
    research_id:STRING,flt_flag:STRING

"""

import sys
import csv
import argparse

def extract_flt_flags(input_file, output_file=None):
    """
    Read input TSV file and extract research_id and flt_ flags where value is True.
    
    Args:
        input_file (str): Path to input tab-separated file
        output_file (str, optional): Path to output tab-separated file. If None, writes to stdout.
    """
    try:
        # Read the input file
        with open(input_file, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            
            # Get column names
            fieldnames = reader.fieldnames
            if not fieldnames:
                raise ValueError("Input file appears to be empty or has no headers")
            
            # Check if research_id column exists
            if 'research_id' not in fieldnames:
                raise ValueError("Input file must contain a 'research_id' column")
            
            # Find all columns that start with 'flt_'
            flt_columns = [col for col in fieldnames if col.startswith('flt_')]
            
            if not flt_columns:
                print("Warning: No columns starting with 'flt_' found in input file")
            # Helper function to convert string to boolean
            def str_to_bool(value):
                """Convert string representations of boolean to actual boolean"""
                if value is None:
                    return False
                value_str = str(value).strip().lower()
                return value_str in ('true', '1', 'yes')

            # Prepare output data
            output_rows = []

            # Iterate through each row
            for row in reader:
                research_id = row['research_id']

                # Check each flt_ column
                for flt_col in flt_columns:
                    if str_to_bool(row.get(flt_col, '')):
                        output_rows.append({
                            'research_id': research_id,
                            'flt_flag': flt_col
                        })

        # Write output file
        if output_file:
            # Write to file
            with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
                fieldnames = ['research_id', 'flt_flag']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')

                # Write header
                writer.writeheader()

                # Write data rows
                writer.writerows(output_rows)
        else:
            # Write to stdout
            fieldnames = ['research_id', 'flt_flag']
            writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter='\t')

            # Write header
            writer.writeheader()

            # Write data rows
            writer.writerows(output_rows)

        if output_rows:
            if output_file:
                print(f"Successfully wrote {len(output_rows)} flag records to {output_file}")
            else:
                print(f"Successfully wrote {len(output_rows)} flag records to stdout", file=sys.stderr)
        else:
            if output_file:
                print(f"No True flags found. Created empty output file: {output_file}")
            else:
                print(f"No True flags found. Wrote empty output to stdout", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Extract research_id and flt_ flag columns from TSV file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python extract_flt_flags.py input.tsv output.tsv
  python extract_flt_flags.py input.tsv  # writes to stdout

Input file format:
  Tab-separated file with 'research_id' column and columns starting with 'flt_'

Output file format:
  Two columns: 'research_id' and 'flt_flag'
  One row for each True value in any flt_ column
        """
    )

    parser.add_argument('input_file', 
                       help='Input tab-separated file')
    parser.add_argument('output_file', 
                       nargs='?',
                       help='Output tab-separated file (optional, writes to stdout if not specified)')
    parser.add_argument('--verbose', '-v', 
                       action='store_true',
                       help='Print verbose output')

    args = parser.parse_args()
    
    if args.verbose:
        print(f"Reading input file: {args.input_file}", file=sys.stderr)
        if args.output_file:
            print(f"Writing output file: {args.output_file}", file=sys.stderr)
        else:
            print("Writing output to stdout", file=sys.stderr)

    extract_flt_flags(args.input_file, args.output_file)
if __name__ == "__main__":
    main()
