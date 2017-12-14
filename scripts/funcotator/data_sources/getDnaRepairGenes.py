#!/usr/bin/env python


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
        writer = csv.writer(f)
        for row in dna_repair_table:
            # Skip empty rows:
            if len(filter(lambda x: x is not '', row)) == 0:
                continue
            # Skip section header rows:
            if len(filter(lambda x: x.strip().lower() == 'top of page', row)) >= 1:
                continue
            writer.writerow(row)

    print 'DONE!'
