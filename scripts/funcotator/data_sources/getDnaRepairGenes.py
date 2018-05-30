#!/usr/bin/env python

# Downloads and formats the list of known DNA Repair Genes from the Wood laborotory website.
# Relies on TableParser

########################################################################
# Imports:

from urllib2 import Request, urlopen, URLError
from TableParser import TableParser

import csv
import datetime

########################################################################
# Constants:

DNA_REPAIR_WEBSITE = 'https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html'

DNA_REPAIR_FILE = 'dnaRepairGenes.' + datetime.datetime.now().strftime('%Y%m%dT%H%M%S') + '.csv'

########################################################################
# Functions:


########################################################################
# Main:

if __name__ == '__main__':

    print 'Reading website: ' + DNA_REPAIR_WEBSITE + '...',

    # Create a web request:
    req = Request(DNA_REPAIR_WEBSITE)

    # Open the URL with our request:
    url = urlopen(req)

    # Create a TableParser:
    tp = TableParser()

    # Load the website into the table parser:
    tp.feed(url.read())

    print 'DONE!'

    # The table we want is the first table on the page:
    dna_repair_table = tp.get_tables()[0]

    print 'Writing outfile: ' + DNA_REPAIR_FILE + '...',

    with open(DNA_REPAIR_FILE, 'wb') as f:
        writer = csv.writer(f, delimiter='|', lineterminator="\n")

        isFirstRow = True
        headerRow = []
        prevRow = []

        for row in dna_repair_table:

            # Strip leading/trailing whitespace and replace newlines:
            strippedRow = map(lambda x: x.strip().replace('\n', ' '), row)

            # Need to fix the first header:
            if isFirstRow:
                strippedRow[0] = "Gene Name"
                headerRow = strippedRow
                isFirstRow = False

            # Skip empty rows:
            if len(filter(lambda x: x is not '', strippedRow)) == 0:
                continue

            # Skip section header rows:
            if len(filter(lambda x: x.lower() == 'top of page', strippedRow)) >= 1:
                continue

            # If we have fewer fields than our header row, we have to grab the
            # "Activity" from the previous row (due to how the table parser works).
            # "Activity" is the second column.
            if len(strippedRow) < len(headerRow):
                strippedRow.insert(1, prevRow[1])

            prevRow = strippedRow
            writer.writerow(strippedRow)

    print 'DONE!'
