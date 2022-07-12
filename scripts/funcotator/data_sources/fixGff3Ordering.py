#!/usr/bin/env python3

# Reorders a GENCODE GTF file to be in ascending order by chromosome, start, end.

########################################################################
# Imports:

import sys
import argparse
import csv

########################################################################
# Constants:

DELIMITER = "\t"
TOTAL_EXPECTED_LINES = 2840283
PERCENT_MARKER = TOTAL_EXPECTED_LINES/20

########################################################################
# Set up parsing for arguments:

parser = argparse.ArgumentParser(description="Creates a new file with each entry in the given GFF3 file in "
                                             "increasing order by location.")
parser.add_argument("GFF3_FILE", help="GFF3 file to read and reorder")

########################################################################
# Functions:

########################################################################
# Main:


if __name__ == "__main__":

    args = parser.parse_args()
    if args.GFF3_FILE is None:
        parser.print_usage()
        sys.exit(1)

    out_writer = csv.writer(sys.stdout, delimiter=DELIMITER, quoting=csv.QUOTE_NONE, quotechar='', lineterminator='\n')

    sys.stderr.write("Processing file: " + args.GFF3_FILE + " ...\n")

    # Open our GTF file:
    with open(args.GFF3_FILE, 'r') as f:

        # Set up our CSV reader:
        gtf_csv_reader = csv.reader(f, delimiter=DELIMITER)

        current_data = []
        current_contig = None

        try:
            while True:
                row = next(gtf_csv_reader)

                if gtf_csv_reader.line_num % PERCENT_MARKER == 0:
                    percent_done = (float(gtf_csv_reader.line_num) / float(TOTAL_EXPECTED_LINES)) * 100.0
                    sys.stderr.write("\tRead " + "{0:1.0f}".format(percent_done) + "%\n")

                if row[0].startswith("#"):
                    out_writer.writerow(row)
                else:
                    if not current_contig:
                        current_contig = row[0]

                    if row[0] != current_contig:
                        sys.stderr.write(f"Writing {current_contig} ...    ")
                        for r in sorted(current_data, key=lambda x: int(x[3])):
                            out_writer.writerow(r)
                        sys.stderr.write(f"DONE!\n")
                        current_contig = row[0]
                        current_data = []
                    
                    current_data.append(row)

        except StopIteration as e:
            sys.stderr.write(f"Writing {current_contig} ...    ")
            for r in sorted(current_data, key=lambda x: int(x[3])):
                out_writer.writerow(r)
            sys.stderr.write(f"DONE!\n")

    sys.stderr.write("DONE!\n")
