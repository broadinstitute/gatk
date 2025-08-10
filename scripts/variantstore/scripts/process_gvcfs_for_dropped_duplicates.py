#!/usr/bin/env python3
"""
Script to process variants by querying GVCF and reblocked GVCF files with bcftools.
Takes a JSON file with variant information and outputs enriched JSON with bcftools results.
"""

import argparse
import json
import re
import subprocess
import sys

# This may be called for data which has never attempted GVCF line fetching, or for data where this has been
# attempted but there were errors we want to retry.
RE_ERROR = re.compile(r"^ERROR:")


def run_bcftools(gvcf_path, chr_name, position):
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
    location = int(variant_obj['location'])
    location_chr = location // 1000000000000

    if location_chr < 23:
        chr_name = f'chr{location_chr}'
    elif location_chr == 23:
        chr_name = 'chrX'
    elif location_chr == 24:
        chr_name = 'chrY'
    else:
        raise Exception(f"location-based chromosome expected to be between 1 and 24: {location}" )

    position = location % 1000000000000
    gvcf_path = variant_obj['gvcf_path']
    reblocked_gvcf = variant_obj['reblocked_gvcf']

    def query_or_lookup(field, gvcf):
        if field not in variant_obj or RE_ERROR.search(variant_obj[field]):
            print(f"re-running bcftools for {variant_obj['sample_name']}:{location}", file=sys.stderr)
            value = run_bcftools(gvcf, chr_name, position)
        else:
            value = variant_obj[field]
        return value

    gvcf_line = query_or_lookup('gvcf_line', gvcf_path)
    reblocked_gvcf_line = query_or_lookup('reblocked_gvcf_line', reblocked_gvcf)

    # Create result object with original fields plus new ones
    result = variant_obj.copy()
    result['chr'] = chr_name
    result['position'] = position
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

    # Output results
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Results written to {args.output}", file=sys.stderr)
    else:
        print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()