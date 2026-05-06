import hail as hl
import argparse
import sys

def main():
    # ==========================================
    # Command Line Arguments
    # ==========================================
    parser = argparse.ArgumentParser(description="Compare sample lists between a Hail VDS and a TSV.")
    parser.add_argument("--vds", required=True, help="Path to the Hail VDS (e.g., gs://your-bucket/dataset.vds)")
    parser.add_argument("--tsv", required=True, help="GCS Path to the TSV file (e.g., gs://your-bucket/samples.tsv)")
    parser.add_argument("--column", required=True, help="Exact name of the TSV column header containing sample IDs")

    args = parser.parse_args()

    # Initialize Hail (quietly)
    hl.init(quiet=True)

    # ==========================================
    # 1. Extract Sample Names from the VDS
    # ==========================================
    print(f"Loading VDS from: {args.vds}")
    try:
        vds = hl.vds.read_vds(args.vds)
    except Exception as e:
        print(f"Error loading VDS: {e}")
        sys.exit(1)

    print("Collecting sample IDs from VDS...")
    vds_samples = set(vds.variant_data.s.collect())
    print(f"Total samples in VDS: {len(vds_samples)}\n")

    # ==========================================
    # 2. Extract Sample Names from the TSV
    # ==========================================
    print(f"Loading TSV from: {args.tsv}")
    try:
        tsv_ht = hl.import_table(args.tsv)
    except Exception as e:
        print(f"Error loading TSV: {e}")
        sys.exit(1)

    # Verify the provided column name actually exists in the TSV
    available_columns = list(tsv_ht.row)
    if args.column not in available_columns:
        print(f"\nError: Column '{args.column}' not found in the TSV.")
        print(f"Available columns are: {available_columns}")
        sys.exit(1)

    print(f"Collecting sample IDs from TSV column '{args.column}'...")
    # Use bracket notation for dynamic column access.
    # We also strip whitespace to prevent silent mismatches (e.g. "Sample1 " vs "Sample1")
    tsv_samples = set([str(s).strip() for s in tsv_ht[args.column].collect() if s is not None])
    print(f"Total samples in TSV: {len(tsv_samples)}\n")

    # ==========================================
    # 3. Compare the Sets
    # ==========================================
    print("--- Comparison Results ---")

    overlapping_samples = vds_samples.intersection(tsv_samples)
    print(f"Samples in BOTH:      {len(overlapping_samples)}")

    only_in_vds = vds_samples - tsv_samples
    print(f"Samples ONLY in VDS:  {len(only_in_vds)}")

    only_in_tsv = tsv_samples - vds_samples
    print(f"Samples ONLY in TSV:  {len(only_in_tsv)}")

    # Optional: Print a few examples of missing samples for quick debugging
    if only_in_tsv:
        print(f"\nExample samples ONLY in TSV (missing from VDS):")
        print(list(only_in_tsv)[:5])

    if only_in_vds:
        print(f"\nExample samples ONLY in VDS (missing from TSV):")
        print(list(only_in_vds)[:5])

if __name__ == "__main__":
    main()