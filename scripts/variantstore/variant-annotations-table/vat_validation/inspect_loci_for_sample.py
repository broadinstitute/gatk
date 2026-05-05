import hail as hl
import argparse
import sys
import os

def get_sample_alleles(vds, chrom, pos, sample_id):
    """
    Finds the specific alleles for a given sample at a target position.
    """
    # 1. Define the locus
    target = hl.locus(chrom, pos, reference_genome="GRCh38")

    # 2. Access variant data
    vd = vds.variant_data

    # 3. Filter down to the exact locus AND the exact sample
    vd_filtered = vd.filter_rows(vd.locus == target)
    vd_filtered = vd_filtered.filter_cols(vd_filtered.s == sample_id)

    # 4. Flatten to an entries Table
    entries_ht = vd_filtered.entries()

    # 5. Translate the sparse Local Genotype (LGT) to actual string alleles.
    entries_ht = entries_ht.annotate(
        actual_alleles = hl.array([
            entries_ht.alleles[entries_ht.LA[entries_ht.LGT[0]]],
            entries_ht.alleles[entries_ht.LA[entries_ht.LGT[1]]]
        ])
    )

    # 6. Collect the result
    result = entries_ht.actual_alleles.collect()

    # 7. Return the alleles, or handle the sparse "empty" state
    if len(result) > 0:
        return list(result[0])
    else:
        return "No variant data (Likely Reference or Missing)"

def main():
    # Set up command line arguments
    parser = argparse.ArgumentParser(description="Query a Hail VDS for a specific sample's alleles at given VIDs.")
    parser.add_argument("--vds", required=True, help="Path to the Hail VDS (e.g., gs://my-bucket/data.vds or /local/path/data.vds)")
    parser.add_argument("--sample", required=True, help="Exact sample ID to search for")
    parser.add_argument("--vids-file", required=True, help="Path to a text file containing VIDs, one per line")

    args = parser.parse_args()

    # Read the VIDs file
    if not os.path.exists(args.vids_file):
        print(f"Error: VIDs file not found at {args.vids_file}")
        sys.exit(1)

    with open(args.vids_file, 'r') as f:
        # Read lines, strip whitespace, and ignore empty lines
        vids = [line.strip() for line in f if line.strip()]

    if not vids:
        print("Error: VIDs file is empty.")
        sys.exit(1)

    print(f"Found {len(vids)} VIDs to process.")

    # Initialize Hail (quietly)
    hl.init(quiet=True)

    # Load the VDS
    print(f"Loading VDS from: {args.vds}...")
    vds = hl.vds.read_vds(args.vds)

    print(f"\n--- Querying alleles for sample: {args.sample} ---")

    # Execute the loop
    for vid in vids:
        try:
            # Parse the VID components
            chrom_raw, pos_str, ref, alt = vid.split('-')

            # Format the arguments
            chrom = f"chr{chrom_raw}"
            pos = int(pos_str)

            # Invoke the function
            alleles = get_sample_alleles(vds, chrom, pos, args.sample)

            # Report
            print(f"VID: {vid:<20} | Locus: {chrom}:{pos:<10} | Alleles: {alleles}")

        except ValueError:
            print(f"VID: {vid:<20} | Error: Invalid VID format. Expected CHR-POS-REF-ALT")
        except Exception as e:
            print(f"VID: {vid:<20} | Error: {str(e)}")

if __name__ == "__main__":
    main()
