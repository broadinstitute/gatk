#!/usr/bin/env python

"""
Script to scrape the data table in the reblocking workspace for reblocked and raw gVCF paths.
"""

NUM_SAMPLES_TO_SCRAPE = 1000000
NUM_SAMPLES_BETWEEN_LOG_MESSAGES = 1000
OUTPUT_FILE = "output.csv"

from terra_notebook_utils import table

samples = table.list_rows('sample', 'AoU_DRC_WGS_reblocking_dragen_378', 'allofus-drc-wgs-prod')
keys = ['research_id', 'reblocked_gvcf', 'gvcf_path']

i = 0

with open(OUTPUT_FILE, 'w') as f:
    for s in samples:
        if i == NUM_SAMPLES_TO_SCRAPE:
            break
        if i == 0:
            f.write('\t'.join(['entity:sample_id'] + keys))
            f.write('\n')

        values = [s.name] + [str(s.attributes.get(k, '')) for k in keys]
        f.write('\t'.join(values))
        f.write('\n')
        i = i + 1
        if i % NUM_SAMPLES_BETWEEN_LOG_MESSAGES == 0:
            print(f"{i}...")