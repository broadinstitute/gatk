#!/usr/bin/env python

"""
Can be run in any terminal in any workspace provided the proxy group being used to access the workspace has adequate
auth to be able to write to the workspace containing the target table.
Note below that this is currently hardcoded to the Foxtrot scale test workspace.
"""

import csv
from terra_notebook_utils import table

UPLOAD_BATCH_SIZE = 1000
WORKSPACE_NAMESPACE = 'gvs-dev'
WORKSPACE_NAME = 'Foxtrot_Batch_testing'
INPUT_FILE = "output.csv"


def put_items(batch):
    table.put_rows(
        table = "samples", # should be 'sample' but that consistently gives cryptic 400 errors
        items = batch,
        workspace = WORKSPACE_NAME,
        workspace_namespace = WORKSPACE_NAMESPACE)


i = 0
batch = []
with open(INPUT_FILE, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        batch.append(dict(row))
        i = i + 1
        if i % UPLOAD_BATCH_SIZE == 0:
            put_items(batch)
            print(f"{i}...")
            batch = []

# if there's anything left to be uploaded, upload it
if batch:
    put_items(batch)