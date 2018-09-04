#!/usr/bin/env python

# This script sorts the Clinvar/HGMD dataset to be in ascending order by chromosome, position.

########################################################################
# Imports:

import csv
import sys
from GenericTsvReader import GenericTsvReader

########################################################################
# Constants:

RAW_FILE_NAME = 'clinvar_hgmd.original.txt'
COLUMN_HEADERS = ['HGMD_ID','SYM', 'CHR', 'LOCATION', 'TYPE', 'ASSEMBLY', 'rs']
SORTED_FILE_NAME = 'clinvar_hgmd.sorted.txt'

ENABLE_STDOUT = True

########################################################################
# Functions:


FLUSH = sys.stdout.flush


def row_comparator(x,y):

    chrCmp = cmp(x['CHR'], y['CHR'])
    if chrCmp == 0:
        startCmp = cmp(int(x['LOCATION']), int(y['LOCATION']))
        if startCmp == 0:
            symCmp = cmp(x['SYM'], y['SYM'])
            if symCmp == 0:
                return cmp(x['HGMD_ID'], y['HGMD_ID'])
            else:
                return symCmp
        else:
            return startCmp
    else:
        return chrCmp


########################################################################
# Main:

if __name__ == '__main__':

    # Read the Clinvar / HGMD file:
    if ENABLE_STDOUT: print 'Reading input file:', RAW_FILE_NAME, "...",
    if ENABLE_STDOUT: FLUSH()
    tsvReader = GenericTsvReader(RAW_FILE_NAME)
    headers = tsvReader.getFieldNames()
    if ENABLE_STDOUT: print('Found headers (input): ' + str(headers))
    if ENABLE_STDOUT: FLUSH()

    clinvarHgmd = []

    if ENABLE_STDOUT: print 'Reading through input file ...',
    if ENABLE_STDOUT: FLUSH()

    # Go through our file:
    for rowNum, row in enumerate(tsvReader):
        clinvarHgmd.append(row)

    if ENABLE_STDOUT: print 'DONE!'
    if ENABLE_STDOUT: FLUSH()

    if ENABLE_STDOUT: print 'Creating output file ...',
    if ENABLE_STDOUT: FLUSH()
    # Set up the output files:
    out_tsv_writer = csv.DictWriter(file(SORTED_FILE_NAME, 'w'), COLUMN_HEADERS, delimiter='\t', lineterminator="\n")
    out_tsv_writer.fieldnames = COLUMN_HEADERS
    out_tsv_writer.writeheader()
    if ENABLE_STDOUT: print 'DONE!'
    if ENABLE_STDOUT: FLUSH()
    
    print 'Sorting dataset ...',
    FLUSH()
    clinvarHgmd.sort(row_comparator)
    print 'DONE!'
    FLUSH()

    print 'Writing dataset ...',
    FLUSH()
    for row in clinvarHgmd:
        row[2] = "chr" + row[2]
        out_tsv_writer.writerow(row)
    print 'DONE!'
    FLUSH()