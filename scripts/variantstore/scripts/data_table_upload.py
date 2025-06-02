#!/usr/bin/env python

"""
Can be run in any analysis terminal provided the associated proxy group has auth to write to the target table.
Note that this is currently hardcoded to write to the Foxtrot Batch / scale test workspace.
"""

import csv
from terra_notebook_utils import table

UPLOAD_BATCH_SIZE = 1000
WORKSPACE_NAMESPACE = 'gvs-dev'
WORKSPACE_NAME = 'Foxtrot_Batch_testing'
INPUT_FILE = "output.csv"


def put_items(batch):
    table.put_rows(
        table = "sample",
        items = batch,
        workspace = WORKSPACE_NAME,
        workspace_namespace = WORKSPACE_NAMESPACE)


i = 0
batch = []
with open(INPUT_FILE, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for csv_row in reader:
        name = csv_row['sample_id']
        attr_names = ['research_id', 'reblocked_gvcf', 'reblocked_gvcf_index']
        attrs = {k: csv_row[k] for k in attr_names}
        row = table.Row(name=name, attributes=attrs)
        batch.append(row)
        i = i + 1
        if i % UPLOAD_BATCH_SIZE == 0:
            put_items(batch)
            print(f"{i}...")
            batch = []

# if there's anything left to be uploaded, upload it
if batch:
    put_items(batch)