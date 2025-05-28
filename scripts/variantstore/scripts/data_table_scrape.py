#!/usr/bin/env python

"""
Intended to be run in a terminal in the workspace from which data is being scraped.
"""

NUM_SAMPLES_TO_SCRAPE = 100000
NUM_SAMPLES_BETWEEN_LOG_MESSAGES = 1000
OUTPUT_FILE = "output.csv"

from terra_notebook_utils import table

samples = table.list_rows('sample')
keys = ['research_id', 'reblocked_gvcf', 'reblocked_gvcf_index']

i = 0

with open(OUTPUT_FILE, 'w') as f:
    for s in samples:
        if i == NUM_SAMPLES_TO_SCRAPE:
            break
        if i == 0:
            f.write('\t'.join(['entity:sample_id'] + keys))
            f.write('\n')

        values = [s.name] + [str(s.attributes[k]) for k in keys]
        f.write('\t'.join(values))
        f.write('\n')
        i = i + 1
        if i % NUM_SAMPLES_BETWEEN_LOG_MESSAGES == 0:
            print(f"{i}...")