""" Creates a new table which maps patient ICD codes in UK Biobank to phecodes.

inputs in BQ
    bigquery hesin table: <HESIN_TABLE>
    previously loaded phecode data

outputs in BQ:
    <OUTPUT_TABLE> with fields <sampleid, datetime, ICD10, Phecode>.
"""
from collections import defaultdict
import csv
from google.cloud import bigquery, storage
import gzip
import json
import os

# CONSTANTS
CWD = os.getcwd()
PHECODES_DATASET = 'shared_data'
PHENO_DATASET = 'ukbb7089_201904'  # TODO: replace with an argument
PARAMS = {
    'phecodes_dataset': PHECODES_DATASET,
    'pheno_dataset': PHENO_DATASET,
    'phecode_dict': PHECODES_DATASET + '.phecode_dictionary',
    'phecode_icd10': PHECODES_DATASET + '.phecode_icd10',
    'gs_bucket': 'ml4h',
    'gs_location': 'data/tmp/',
    'output_schema': os.path.join(CWD, 'phecode_mapping.json'),
    'output_file': PHENO_DATASET + '_phecode_mapping.csv.gz',
    'output_table': PHENO_DATASET + '.phecode_mapping',
}
client = bigquery.Client() #should already be hooked up to your default project


def get_ranged_phecodes():
    """returns a sorted list of [(min_icd10, max_icd10, phecode), ]
    sorted by min_icd10, then max_icd10
    from raw phecode table"""

    # create list of ranged elements (note, same range may appear twice with
    # different phecode
    ranged_lookups = []  # [(low,high,phecode),...]
    query_job = client.query(
        """
        SELECT REPLACE(icd10,'.','') as icd10, phecode
        FROM %(phecode_icd10)s
        WHERE regexp_contains(icd10,'-')
        """ % PARAMS
    )
    for row in query_job:
        (low, high) = tuple(row['icd10'].split('-'))  # (start icd10,end icd10)
        ranged_lookups.append((low, high, row['phecode']))
    # convert to sorted list [(lowest, high, phecode), (...), (...),...]
    return sorted(ranged_lookups, key=lambda x: (x[0], x[1]))


def identify_phecode_from_ranged_list(sorted_ranged_list, icd10):
    """turns icd10 into a set of phecodes from sorted, ranged phecode list"""
    icd10 = icd10.replace('.', '')

    phecodes = set()
    for (low, high, phecode) in sorted_ranged_list:
        if icd10 < low:
            continue
        elif icd10 <= high:
            phecodes.add(phecode)
        elif icd10 > high:
            continue
        else:
            raise Exception(
                "not catching icd10 %s in range %s, %s" %
                (icd10, low, high),
            )
    return phecodes


def create_phecode_match_file():
    """function that produces output table of eid, record_id, admidate,
                diag_icd10, phecode

    the original phecode table is redundant in many ways:
    icd10 can be mapped to multiple phecodes exactly (ok)
    icd10 can be mapped via a range to multiple phecodes (ok)
    icd10 can be mapped to the same phecode once by range and once by exact
    match (problem that we want to avoid)

    output is written to local file
    """

    # grab ranged phecodes
    ranged_sorted_phecodes = get_ranged_phecodes()  # TODO: pass in client

    # now grab exact phecodes, excluding ones covered in the range already
    exact_phecodes = defaultdict(set)  # dict of exact icd -> set(phecode),
    query_job = client.query(
        """
        SELECT replace(icd10,'.','') as icd10, phecode
        from %(phecode_icd10)s
        where regexp_contains(icd10,'-') = false
        """ % PARAMS
    )

    for row in query_job:
        matched_phecodes = identify_phecode_from_ranged_list(
                           ranged_sorted_phecodes,
                           row['icd10'],
        )
        if row['phecode'] in matched_phecodes:
            print(
                '(%s,%s) icd10, phecode redundantly found in range'
                % (row['icd10'], row['phecode']),
            )
        else:
            exact_phecodes[row['icd10']].add(row['phecode'])

    #read row from hesin table, match to phecode, dump to file
    query_job = client.query(
        """
        SELECT eid, record_id, admidate, diag_icd10,
        replace(diag_icd10,'.','') as icd10
        FROM %(pheno_dataset)s.hesin
        where diag_icd10 is not null
        """ % PARAMS
    )
    count = 0

    #read schema from json to make sure we get the ordering right
    with open(PARAMS['output_schema']) as f:
        schema = json.load(f)
    #pass ordering into fieldnames for writing
    fieldnames = [fn['name'] for fn in schema]

    #find phecode matches, write to a file
    with gzip.open(PARAMS['output_file'], mode='wt') as csv_file:
        csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        for row in query_job:
            count += 1
            if count % 1000 == 0:
                print(count)
            #one at a time, should be slow as hell
            matched_phecodes = identify_phecode_from_ranged_list(
                            ranged_sorted_phecodes,
                            row['icd10'],
            )
            all_matches = matched_phecodes | exact_phecodes.get(
                row['icd10'],
                set(),
            )
            if len(all_matches):
                output_dict = {}
                for fn in fieldnames:
                    if fn == 'phecode':
                        continue
                    else:
                        output_dict[fn] = row[fn]
                for match in all_matches:
                    output_dict['phecode'] = match
                csv_writer.writerow(output_dict)

#not operational
#def upload_file_to_db():
#    with open(PARAMS['output_schema']) as f:
#        schema = json.load(f)
#    #pass ordering into fieldnames for writing
#    #pass ordering,types, mode into jobconfig schema for later bigquery writing
#    bq_schema = [bigquery.SchemaField(fn['name'],
#                                      fn['type'],
#                                      fn['mode']) for fn in schema]
#
#
#    table_ref=client.dataset(PHENO_DATASET).table(PARAMS['output_table'])
#    table = client.create_table(bigquery.Table('broad-ml4cvd.lubitz.phecode_mapping', schema=bq_schema))
#    with open(os.path.join(CWD, PARAMS['output_file']), 'rb') as f:
#        job=client.load_table_from_file(
#            f,table)


if __name__ == '__main__':
    test = {
        ('a', 'b'): 2, ('a', 'c'): 3, ('a1', 'b'): 2,
        ('a',  'b1'): 2, ('a1', 'b1'): 3,
    }
    correct = [
        ('a', 'b', 2), ('a', 'b1', 2), ('a', 'c', 3),
        ('a1', 'b', 2), ('a1', 'b1', 3),
    ]
    test_list = [(k[0], k[1], v) for k, v in test.items()]
    sorted_test_list = sorted(test_list, key=lambda x: (x[0], x[1]))
    print(sorted_test_list)
    assert sorted_test_list == correct
    print(identify_phecode_from_ranged_list(sorted_test_list, 'a'))
    create_phecode_match_file()
    #upload_file_to_db()
