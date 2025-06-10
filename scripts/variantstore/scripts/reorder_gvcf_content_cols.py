import argparse
import csv
import sys


def main():
    parser = argparse.ArgumentParser(description='Establish column order in GVCF content TSV file')
    parser.add_argument('tsv_file', help='Input TSV file containing variant objects')

    args = parser.parse_args()

    # Load input JSON
    try:
        # https://stackoverflow.com/a/33002011
        with open(args.tsv_file, 'r') as infile:
            fieldnames = ["sample_name", "sample_id", "chr", "input_position", "input_ref",
                          "input_alt", "gvcf_path", "reblocked_gvcf", "gvcf_line", "reblocked_gvcf_line"]
            writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=fieldnames)
            # reorder the header first
            writer.writeheader()
            for row in csv.DictReader(infile, delimiter="\t"):
                # writes the reordered rows to the new file
                writer.writerow(row)
    except FileNotFoundError:
        print(f"Error: Input file '{args.tsv_file}' not found.", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()