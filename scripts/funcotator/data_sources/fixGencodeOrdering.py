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
CONTIG_PRINT_ORDER_LIST = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                           "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                           "chr22", "chrX", "chrY", "chrM"]

########################################################################
# Set up parsing for arguments:

parser = argparse.ArgumentParser(description="Creates a new file with each Gene in the given GENCODE_GTF_FILE in "
                                             "increasing order by location.")
parser.add_argument("GENCODE_GTF_FILE", help="Gencode GTF file to read and reorder")

########################################################################
# Functions:


def write_contig(contig, contig_dictionary, out_writer):
    sys.stderr.write("Contig: " + contig + ":\n")

    sys.stderr.write("\tSorting ...\n")

    gene_list = contig_dictionary[contig]

    # sort our contig by the start gene column(3)
    gene_list.sort(key=lambda x: int(x[0][3]))

    sys.stderr.write("\tWriting ...\n")
    # print what we have so far
    for g in gene_list:
        for element in g:
            out_writer.writerow(element)


########################################################################
# Main:


if __name__ == "__main__":

    args = parser.parse_args()
    if args.GENCODE_GTF_FILE is None:
        parser.print_usage()
        sys.exit(1)

    out_writer = csv.writer(sys.stdout, delimiter=DELIMITER, quoting=csv.QUOTE_NONE, quotechar='', lineterminator='\n')

    sys.stderr.write("Processing file: " + args.GENCODE_GTF_FILE + " ...\n")

    # Open our GTF file:
    with open(args.GENCODE_GTF_FILE, 'r') as f:

        # Set up our CSV reader:
        gtf_csv_reader = csv.reader(f, delimiter=DELIMITER)

        row = next(gtf_csv_reader)
        # Save the header:
        while row[0].startswith("#"):
            out_writer.writerow(row)
            row = next(gtf_csv_reader)

        gene = [row] 

        # contig -> gene_list
        contig_dictionary = {gene[0][0]: []}

        # Read until there's nothing left:
        try:
            while True:
                row = next(gtf_csv_reader) 

                if gtf_csv_reader.line_num % PERCENT_MARKER == 0:
                    percent_done = (float(gtf_csv_reader.line_num) / float(TOTAL_EXPECTED_LINES)) * 100.0
                    sys.stderr.write("\tRead " + "{0:1.0f}".format(percent_done) + "%\n")

                if row[2] == "gene":

                    # Add our gene to the contig:
                    if gene[0][0] not in contig_dictionary:
                        contig_dictionary[gene[0][0]] = [gene]
                    else:
                        contig_dictionary[gene[0][0]].append(gene)

                    # create new gene
                    gene = []

                # add this line to the new gene
                gene.append(row)

        except StopIteration as e:
            pass

        sys.stderr.write("File Read Complete.\n")

        # Get the contigs that we've found:
        contigs = list(contig_dictionary.keys())

        # Go through and print our prioritized contig list first:
        for contig in CONTIG_PRINT_ORDER_LIST:
            write_contig(contig, contig_dictionary, out_writer)
            contigs.remove(contig)

        # Print the remaining ones:
        for contig in sorted(contigs, key=str.lower):
            write_contig(contig, contig_dictionary, out_writer)

    sys.stderr.write("DONE!\n")
