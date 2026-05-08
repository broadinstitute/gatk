import hail as hl
import argparse
import sys

def main():
    # ==========================================
    # Command Line Arguments
    # ==========================================
    parser = argparse.ArgumentParser(description='Compare the samples in a Foxtrot "minus" VDS, an Epsilon VDS, and the v9 r2 manifest file')
    parser.add_argument("--foxtrot-minus-vds", required=True, help='Path to the Foxtrot "minus" VDS (e.g., gs://your-bucket/foxtrot_minus.vds)')
    parser.add_argument("--epsilon-vds", required=True, help='Path to the Epsilon VDS (e.g., gs://your-bucket/epsilon.vds)')
    parser.add_argument("--foxtrot-r2-manifest", required=True, help="GCS Path to the Foxtrot r2 manifest file (e.g., gs://your-bucket/manifest.tsv)")
    parser.add_argument("--column", required=True, help="Exact name of the TSV column header containing sample IDs")

    args = parser.parse_args()

    # Initialize Hail (quietly)
    hl.init(quiet=True)

    # ==========================================
    # 1. Extract Sample Names from the VDS
    # ==========================================
    print(f"Loading Foxtrot minus VDS from: {args.foxtrot_minus_vds}")
    try:
        foxtrot_minus_vds = hl.vds.read_vds(args.foxtrot_minus_vds)
    except Exception as e:
        print(f"Error loading Foxtrot minus VDS: {e}")
        sys.exit(1)

    print(f"Loading Epsilon VDS from: {args.epsilon_vds}")
    try:
        epsilon_vds = hl.vds.read_vds(args.epsilon_vds)
    except Exception as e:
        print(f"Error loading Epsilon VDS: {e}")
        sys.exit(1)

    print("Collecting sample IDs from Foxtrot minus VDS...")
    foxtrot_minus_samples = set(foxtrot_minus_vds.variant_data.s.collect())
    print(f"Total samples in Foxtrot minus VDS: {len(foxtrot_minus_samples)}\n")

    print("Collecting sample IDs from Epsilon VDS...")
    epsilon_samples = set(epsilon_vds.variant_data.s.collect())
    print(f"Total samples in Epsilon VDS: {len(epsilon_samples)}\n")

    # ==========================================
    # 2. Extract Sample Names from the TSV
    # ==========================================
    print(f"Loading TSV from: {args.foxtrot_r2_manifest}")
    try:
        tsv_ht = hl.import_table(args.foxtrot_r2_manifest)
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
    vds_samples = foxtrot_minus_samples | epsilon_samples

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