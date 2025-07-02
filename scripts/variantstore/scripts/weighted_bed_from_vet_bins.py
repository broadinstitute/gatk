import argparse
import csv
import math
import pandas as pd
import re
import sys


def read_contig_sizes(args):
    contig_sizes = {}

    # Parse the sequence dictionary to get contig lengths
    with open(args.reference_dict, 'r') as interval_list_file:
        for header_line in interval_list_file:
            if header_line.startswith("@SQ"):
                # look for the contig and length fields
                contig_result = re.search(r'SN:(\S+)', header_line)
                length_result = re.search(r'LN:(\d+)', header_line)

                if contig_result and length_result:
                    contig = contig_result.group(1)
                    contig_length = length_result.group(1)
                    contig_sizes[contig] = int(contig_length)
    return contig_sizes


def main(args, bed_list):
    contig_sizes = read_contig_sizes(args)

    binsize_kb = args.binsize_kb

    last_contig = ""
    last_ending = 0
    interval_size = binsize_kb * 1000

    def maybe_say(s):
        if False:
            print(s, sys.stderr)


    output = open(args.outfile, 'w') if args.outfile else sys.stdout
    # logic is simple.  For every row, print it to the new file.  EXCEPT, we also want to
    # detect gaps in between intervals and fill them with 0 weight intervals.  That way,
    # when it's read in we'll have a solid block of data for an entire contig and can put
    # it in an array for fast looking
    for row in bed_list:
        contig = row[0]
        start = row[1]
        ending = row[2]
        # Logic to detect gaps in the intervals is only valid WITHIN a contig
        if contig == last_contig:
            # the intervals are read as [,) so the final entry in the previous one should
            # match the first entry in the next one if there are no gaps
            if start != last_ending:
                maybe_say(f"detected gap from {contig} {last_ending} to {start}")
                size_of_gap = start - last_ending
                num_empty_entries = size_of_gap / interval_size
                maybe_say(f"Will create {num_empty_entries} 0 weight entries here")
                for s in range(last_ending, start, interval_size):
                    maybe_say(f"\t{contig} {s} - {s + interval_size}")
                    new_row = [contig, str(s), str(s + interval_size), ".", "0"]
                    output.write("\t".join(new_row) + "\n")
        else:
            # see if there's anything on the end of the last contig to extend out to
            if last_contig in contig_sizes:
                # we'll want to round up the listed ending to the next block in case we
                # haven't reached it yet
                new_ending_block = int(math.ceil(contig_sizes[last_contig] / float(interval_size))) * interval_size
                maybe_say(f"Ending detected for contig: {last_contig} rounded up to block {new_ending_block}.  Last block written was {last_ending}")
                for s in range(last_ending, new_ending_block, interval_size):
                    maybe_say(f"\t{last_contig} {s} - {s + interval_size}")
                    new_row = [last_contig, str(s), str(s + interval_size), ".", "0"]
                    output.write("\t".join(new_row) + "\n")
            # make sure we always start at 0 on a new contig and search for gaps
            maybe_say(f"detected gap at beginning of {contig}")
            last_ending = 0
            size_of_gap = start
            num_empty_entries = size_of_gap / interval_size
            maybe_say(f"Would create {num_empty_entries} 0 weight entries here")
            for s in range(last_ending, start, interval_size):
                maybe_say(f"\t{contig} {s} - {s + interval_size}")
                new_row = [contig, str(s), str(s + interval_size), ".", "0"]
                output.write("\t".join(new_row) + "\n")

        last_contig = contig
        last_ending = ending
        # now that we've inserted any fake rows necessary to fill gaps, write the current one
        output.write("\t".join([str(i) for i in row]) + "\n")

    #end for loop
    # check for anything leftover at the end of the last contig
    if last_contig in contig_sizes:
        # we'll want to round up the listed ending to the next block in case we
        # haven't reached it yet
        new_ending_block = int(math.ceil(contig_sizes[last_contig] / float(interval_size))) * interval_size
        maybe_say(f"Ending detected for contig: {last_contig} rounded up to block {new_ending_block}.  Last block written was {last_ending}")
        for s in range(last_ending, new_ending_block, interval_size):
            maybe_say(f"\t{last_contig} {s} - {s + interval_size}")
            new_row = [last_contig, str(s), str(s + interval_size), ".", "0"]
            output.write("\t".join(new_row) + "\n")
    maybe_say("Padding of bed file completed")


def load_contig_mapping(mapping_file):
    contigs = {}
    with open(mapping_file, mode='r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) == 2:
                contigs[int(row[1])] = row[0]
            else:
                raise ValueError(f"Invalid line format: {row}")
    return contigs


def raw_bins_to_bed_list(args, contig_map):
    location_offset = 1000000000000

    w = pd.read_csv(args.input_bin_data, dtype={'bin': int, 'entries':int})

    w['contig'] = (w['bin'].astype(int) / location_offset).astype(int).map(contig_map)
    w['start_position'] = w['bin'].astype(int) - (w['bin'] / location_offset).astype(int) * location_offset
    w['end_position'] = w['start_position'] + args.binsize_kb * 1000
    w['name'] = "."
    o = w[['contig', 'start_position','end_position', 'name', 'entries' ]]

    return o.values.tolist()


def parse_args():
    parser = argparse.ArgumentParser(description='Create a padded BED file with gap intervals.')
    parser.add_argument('--binsize-kb', type=float, default=1,
                        help='Bin size in kilobases (default: 1)')
    parser.add_argument('--input-bin-data', type=str, required=True,
                        help='Raw bin data input')
    parser.add_argument('--outfile', type=str,
                        help='Output padded BED file (default: stdout)')
    parser.add_argument('--reference-dict', type=str, required=True,
                        help='Reference sequence dictionary file')
    parser.add_argument('--contig-mapping', type=str, required=True,
                        help='Contig mapping file')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Load the contig mapping
    contig_map = load_contig_mapping(args.contig_mapping)

    bed_list = raw_bins_to_bed_list(args, contig_map)
    main(args, bed_list)
