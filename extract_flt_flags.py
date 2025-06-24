#!/usr/bin/env python3
"""
Script to extract research_id and flt_ flag columns from a tab-separated file.

For each row, if any flt_ column has a value of True, output a line with:
<research_id><tab><flt_column_name>

Usage: python extract_flt_flags.py <input_file> <output_file>
"""

import sys
import pandas as pd
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
        df = pd.read_csv(input_file, sep='\t')
        
        # Check if research_id column exists
        if 'research_id' not in df.columns:
            raise ValueError("Input file must contain a 'research_id' column")
        
        # Find all columns that start with 'flt_'
        flt_columns = [col for col in df.columns if col.startswith('flt_')]
        
        if not flt_columns:
            print("Warning: No columns starting with 'flt_' found in input file")
        
        # Convert flt_ columns to boolean, handling string representations
        for col in flt_columns:
            # Convert string representations of boolean to actual boolean
            df[col] = df[col].astype(str).str.strip().str.lower().map({
                'true': True,
                'false': False,
                '1': True,
                '0': False,
                'yes': True,
                'no': False
            })
        
        # Prepare output data
        output_rows = []
        
        # Iterate through each row
        for idx, row in df.iterrows():
            research_id = row['research_id']
            
            # Check each flt_ column
            for flt_col in flt_columns:
                if row[flt_col] is True:  # Explicitly check for True
                    output_rows.append({
                        'research_id': research_id,
                        'flt_flag': flt_col
                    })
        
        # Create output DataFrame and write to file
        if output_rows:
            output_df = pd.DataFrame(output_rows)
            output_df.to_csv(output_file, sep='\t', index=False)
            print(f"Successfully wrote {len(output_rows)} flag records to {output_file}")
        else:
            # Create empty output file with headers
            output_df = pd.DataFrame(columns=['research_id', 'flt_flag'])
            output_df.to_csv(output_file, sep='\t', index=False)
            print(f"No True flags found. Created empty output file: {output_file}")
            
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: Input file '{input_file}' is empty")
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
