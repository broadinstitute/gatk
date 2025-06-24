#!/usr/bin/env python3
"""
Script to extract research_id and flt_ flag columns from a tab-separated file.

For each row, if any flt_ column has a value of True, output a line with:
<research_id><tab><flt_column_name>

Usage: python extract_flt_flags.py <input_file> <output_file>
"""

import sys
import csv
import argparse

def extract_flt_flags(input_file, output_file):
    """
    Read input TSV file and extract research_id and flt_ flags where value is True.
    
    Args:
        input_file (str): Path to input tab-separated file
        output_file (str): Path to output tab-separated file
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
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['research_id', 'flt_flag']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')

            # Write header
            writer.writeheader()

            # Write data rows
            writer.writerows(output_rows)

        if output_rows:
            print(f"Successfully wrote {len(output_rows)} flag records to {output_file}")
        else:
            print(f"No True flags found. Created empty output file: {output_file}")

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
                       help='Output tab-separated file')
    parser.add_argument('--verbose', '-v', 
                       action='store_true',
                       help='Print verbose output')

    args = parser.parse_args()
    
    if args.verbose:
        print(f"Reading input file: {args.input_file}")
        print(f"Writing output file: {args.output_file}")

    extract_flt_flags(args.input_file, args.output_file)
if __name__ == "__main__":
    main()
