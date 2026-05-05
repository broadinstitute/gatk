import hail as hl
import argparse
import sys
import os

def count_distinct_alleles(vds, chrom, pos):
    """
    Returns the number of distinct alleles (REF + all ALTs) at a specific locus
    by checking the length of the alleles array.
    """
    # 1. Safely define a strictly inclusive literal interval
    interval_str = f"[{chrom}:{pos}-{pos}]"
    interval = hl.parse_locus_interval(interval_str, reference_genome="GRCh38")

    # 2. We ONLY need the sparse variant_data for this!
    vd = vds.variant_data
    vd_filtered = hl.filter_intervals(vd, [interval])

    # 3. Collect the row keys (which automatically includes 'locus' and 'alleles')
    rows = vd_filtered.select_rows().rows().collect()

    # 4. Handle empty states
    if not rows:
        return 0, []

    # 5. Extract the alleles array from the first (and likely only) row at this locus
    alleles = rows[0].alleles

    return len(alleles), list(alleles)


def main():
    # Set up command line arguments
    parser = argparse.ArgumentParser(description="Count distinct alleles at specified loci in a Hail VDS.")
    parser.add_argument("--vds", required=True, help="Path to the Hail VDS")
    parser.add_argument("--loci_file", required=True, help="Path to a text file containing loci (Format: CHR-POS, e.g., 20-345678)")

    args = parser.parse_args()

    # Read the loci file
    if not os.path.exists(args.loci_file):
        print(f"Error: Loci file not found at {args.loci_file}")
        sys.exit(1)

    with open(args.loci_file, 'r') as f:
        loci_lines = [line.strip() for line in f if line.strip()]

    if not loci_lines:
        print("Error: Loci file is empty.")
        sys.exit(1)

    print(f"Found {len(loci_lines)} loci to process.")

    # Initialize Hail (quietly)
    hl.init(quiet=True)

    # Load the VDS
    print(f"Loading VDS from: {args.vds}...\n")
    vds = hl.vds.read_vds(args.vds)

    # Execute the loop
    for locus_str in loci_lines:
        try:
            # Parse the CHR-POS components
            parts = locus_str.split('-')
            if len(parts) < 2:
                raise ValueError("Format must be CHR-POS")

            chrom_raw = parts[0]
            pos_str = parts[1]

            # Format the arguments
            chrom = f"chr{chrom_raw}"
            pos = int(pos_str)

            # Invoke our much simpler function
            num_alleles, alleles_list = count_distinct_alleles(vds, chrom, pos)

            # Print the summary
            print(f"--- Locus: {chrom}:{pos} ---")
            if num_alleles == 0:
                print("  No variants found in variant_data.\n")
            else:
                print(f"  Distinct alleles count: {num_alleles}")
                print(f"  Alleles array: {alleles_list}\n")

        except ValueError as e:
            print(f"Locus: {locus_str:<15} | Error: {str(e)}")
        except Exception as e:
            print(f"Locus: {locus_str:<15} | Error: {str(e)}")

if __name__ == "__main__":
    main()