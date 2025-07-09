import argparse
import re
import sys


# Gather all the contigs from the sequence dictionary
def map_contigs(args):
    # Grab whatever is identified as SN:contig_name in the file
    contig_matcher = re.compile(r"SN:(\S+)")

    all_raw_contigs = []
    with open(args.sequence_dictionary) as sequence_dictionary:
        for line in sequence_dictionary:
            line = line.strip()
            line_match = contig_matcher.search(line)
            if line_match:
                all_raw_contigs.append(line_match.group(1))

    print(f'{len(all_raw_contigs)} raw contigs found', file=sys.stderr)

    # Make sure it's deduped and sorted list. Convert to set then back to a list, which we'll sort
    deduped = set(all_raw_contigs)
    sorted_contigs = list(deduped)
    sorted_contigs.sort()

    print(f'{len(sorted_contigs)} deduped and sorted contigs found', file=sys.stderr)

    contig_num = 1
    # assign each contig a basic, sequential number
    # write to file if output is specified, otherwise to stdout
    if args.output:
        with open(args.output, "w") as contig_mapping_file:
            for contig in sorted_contigs:
                contig_mapping_file.write(f'{contig}\t{contig_num}\n')
                contig_num += 1
    else:
        for contig in sorted_contigs:
            print(f'{contig}\t{contig_num}')
            contig_num += 1


def parse_args():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Generate custom reference mappings from a sequence dictionary.')
    parser.add_argument('sequence_dictionary', type=str, help='Path to the sequence dictionary file')
    parser.add_argument('--output', type=str, help='Output file name for mappings (if not specified, write to stdout)')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    map_contigs(args)
