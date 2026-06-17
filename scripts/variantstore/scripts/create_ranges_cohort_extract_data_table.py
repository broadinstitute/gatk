# -*- coding: utf-8 -*-
import uuid
import datetime
import argparse
import io
import pybedtools
import re

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig

import utils

JOBS = []

#
# CONSTANTS
#
REF_TABLE_PREFIX = "ref_ranges_"
VET_TABLE_PREFIX = "vet_"
SAMPLES_PER_PARTITION = 4000


# temp-table-uuid
output_table_prefix = str(uuid.uuid4()).split("-")[0]
print(f"running with prefix {output_table_prefix}")

INTERVAL_TABLE_NAME = f"{output_table_prefix}_intervals"
# Intervals up to this count are inlined directly into the SQL WHERE clause.
# Above this threshold the intervals are loaded into a temporary BigQuery table
# and filtered via an EXISTS subquery, removing the effective query-size limit.
INTERVAL_TEMP_TABLE_THRESHOLD = 5000

# Total bases in hg38 chromosomes 1-22, X (23), Y (24) — the chromosomes GVS encodes.
HG38_GENOME_SIZE = 3_088_269_832

REF_VET_TABLE_COUNT = -1
# noinspection PyTypeChecker
client: bigquery.Client = None
default_config = None

EXTRACT_SAMPLE_TABLE = f"{output_table_prefix}_sample_names"

COMPRESSED_REFS_HELPER_FUNCTION_DEFS="""CREATE TEMP FUNCTION intToState(state INT64)
RETURNS string
AS (
    CASE state
WHEN 7 THEN 'v'
WHEN 8 THEN '*'
WHEN 9 THEN 'm'
WHEN 10 THEN 'u'
ELSE CAST(state as string)
END
);

CREATE TEMP FUNCTION UnpackRefRangeInfo(superpackEntry int64)
RETURNS STRUCT<location INT64, len INT64, state string>
AS (
    STRUCT(
        1000000000000 * ((superpackEntry >> 48) & 0xFFFF) + ((superpackEntry >> 16) & 0xFFFFFFFF),
        (superpackEntry >> 4) & 0xFFF,
        intToState(superpackEntry & 0xF))
);
"""

CHROM_MAP = {'chr1': '1', 'chr2': '2', 'chr3': '3', 'chr4': '4', 'chr5': '5', 'chr6': '6', 'chr7': '7', 'chr8': '8', 'chr9': '9', 'chr10': '10', 'chr11': '11', 'chr12': '12', 'chr13': '13', 'chr14': '14', 'chr15': '15', 'chr16': '16', 'chr17': '17', 'chr18': '18', 'chr19': '19', 'chr20': '20', 'chr21': '21', 'chr22': '22', 'chrX': '23', 'chrY': '24', 'chrM': '25'}


def get_partition_range(i):
    if i < 1 or i > REF_VET_TABLE_COUNT:
        raise ValueError(f"out of partition range")

    return {'start': (i - 1) * SAMPLES_PER_PARTITION + 1, 'end': i * SAMPLES_PER_PARTITION}


def get_samples_for_partition(sample_ids, i):
    return [s for s in sample_ids if get_partition_range(i)['start'] <= s <= get_partition_range(i)['end']]


def split_lists(samples, n):
    return [samples[i * n:(i + 1) * n] for i in range((len(samples) + n - 1) // n)]


def load_sample_names(sample_names_to_extract, fq_temp_table_dataset):
    schema = [bigquery.SchemaField("sample_name", "STRING", mode="REQUIRED")]
    fq_sample_table = f"{fq_temp_table_dataset}.{EXTRACT_SAMPLE_TABLE}"

    job_labels = client._default_query_job_config.labels
    job_labels["gvs_query_name"] = "load-sample-names"

    job_config = bigquery.LoadJobConfig(source_format=bigquery.SourceFormat.CSV, skip_leading_rows=0, schema=schema,
                                        labels=job_labels)

    with open(sample_names_to_extract, "rb") as source_file:
        job = client.load_table_from_file(source_file, fq_sample_table, job_config=job_config)

    job.result()  # Waits for the job to complete.

    # setting the TTL needs to be done as a second API call
    table = bigquery.Table(fq_sample_table, schema=schema)
    expiration = datetime.datetime.utcnow() + datetime.timedelta(hours=TEMP_TABLE_TTL_HOURS)
    table.expires = expiration
    client.update_table(table, ["expires"])

    return fq_sample_table


def get_all_sample_ids(fq_destination_table_samples, only_output_vet_tables, fq_sample_mapping_table):
    if only_output_vet_tables:
        sql = f"select sample_id from `{fq_sample_mapping_table}` WHERE is_control = false AND withdrawn IS NULL"
        sample_table = fq_sample_mapping_table
    else:
        sql = f"select sample_id from `{fq_destination_table_samples}`"
        sample_table = fq_destination_table_samples

    query_return = utils.execute_with_retry(client, "read cohort sample table", sql)
    JOBS.append({'job': query_return['job'], 'label': query_return['label']})
    sample_ids = [row.sample_id for row in list(query_return['results'])]
    sample_ids.sort()
    print(f"Discovered {len(sample_ids)} samples in {sample_table}...")
    return sample_ids


def create_extract_samples_table(control_samples, fq_destination_table_samples, fq_sample_name_table,
                                 fq_sample_mapping_table, honor_withdrawn, enable_extract_table_ttl):
    ttl = ""
    if enable_extract_table_ttl:
        ttl = "OPTIONS( expiration_timestamp=TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL 14 DAY))"

    sql = f"""
    CREATE OR REPLACE TABLE `{fq_destination_table_samples}` 
    {ttl}
    AS (
        SELECT m.sample_id, m.sample_name, m.is_loaded, {"m.withdrawn," if honor_withdrawn else "NULL as withdrawn,"} m.is_control FROM `{fq_sample_name_table}` s JOIN
        `{fq_sample_mapping_table}` m ON (s.sample_name = m.sample_name) WHERE
             m.is_loaded IS TRUE AND m.is_control = {control_samples}
             {"AND m.withdrawn IS NULL" if honor_withdrawn else ""}
    )
    """
    print(sql)

    query_return = utils.execute_with_retry(client, "create extract sample table", sql)
    JOBS.append({'job': query_return['job'], 'label': query_return['label']})
    return query_return['results']


def load_intervals_to_temp_table(intervals, fq_temp_table_dataset, table_prefix):
    """Load a pybedtools BedTool of intervals into a temporary BigQuery table.

    Each interval is stored as (location_start INT64, location_end INT64, chrom_index INT64) using
    the GVS location encoding: chrom_index * 1_000_000_000_000 + position.

    chrom_index is stored as a separate column so it can be used as an equality join key,
    allowing BigQuery to hash-partition the join on chromosome before applying the BETWEEN
    range filter. Without the equality key, BQ falls back to a nested loop join against
    all interval rows, which times out on large datasets.

    Intervals are written row-by-row to a CSV buffer and uploaded via a BQ load
    job rather than the streaming insert API.  This avoids accumulating a large
    Python list in memory and is far faster for very large interval lists (e.g.
    tens of millions of intervals).

    Returns the fully-qualified name of the created table.
    """
    fq_interval_table = f"{fq_temp_table_dataset}.{table_prefix}__{INTERVAL_TABLE_NAME}"
    schema = [
        bigquery.SchemaField("location_start", "INT64", mode="REQUIRED"),
        bigquery.SchemaField("location_end", "INT64", mode="REQUIRED"),
        bigquery.SchemaField("chrom_index", "INT64", mode="REQUIRED"),
    ]

    # Write intervals as CSV directly into a buffer — avoids building a large
    # Python list and keeps peak memory proportional to the encoded CSV size
    # (~20 bytes/row) rather than the size of Python dict objects (~250 bytes/row).
    buf = io.BytesIO()
    rows_written = 0
    skipped = 0
    for interval in intervals:
        if interval.chrom not in CHROM_MAP:
            skipped += 1
            continue
        chrom_num = int(CHROM_MAP[interval.chrom])
        loc_start = chrom_num * 1_000_000_000_000 + int(interval.start)
        loc_end   = chrom_num * 1_000_000_000_000 + int(interval.end)
        buf.write(f"{loc_start},{loc_end},{chrom_num}\n".encode())
        rows_written += 1

    if skipped > 0:
        print(f"Warning: skipped {skipped} interval(s) with unrecognized chromosome names")

    buf.seek(0)
    job_config = bigquery.LoadJobConfig(
        schema=schema,
        source_format=bigquery.SourceFormat.CSV,
        skip_leading_rows=0,
        write_disposition=bigquery.WriteDisposition.WRITE_EMPTY,
    )
    job = client.load_table_from_file(buf, fq_interval_table, job_config=job_config)
    job.result()  # Wait for the load job to complete.

    # Set the TTL as a second API call (same pattern used by load_sample_names).
    table = bigquery.Table(fq_interval_table, schema=schema)
    expiration = datetime.datetime.utcnow() + datetime.timedelta(hours=TEMP_TABLE_TTL_HOURS)
    table.expires = expiration
    client.update_table(table, ["expires"])

    print(f"Loaded {rows_written} intervals into temporary table {fq_interval_table}")
    return fq_interval_table


def get_location_filters(interval_list, fq_temp_table_dataset, table_prefix, padding=1000, skip_filter_threshold=0.5):
    """Return a SQL WHERE clause fragment for filtering by genomic location.

    Expands every interval by `padding` bp on each side and merges overlaps, then chooses
    a strategy based on merged interval count and genome coverage:

    * merged count <= INTERVAL_TEMP_TABLE_THRESHOLD: inline SQL WHERE clause
    * merged count > threshold, coverage < skip_filter_threshold: temp-table EXISTS subquery
    * coverage >= skip_filter_threshold: return "" (no filtering; caller does a full prepare)

    Returns "" when no interval list is provided (no location filtering).
    """
    if not interval_list:
        return ""

    raw_intervals = pybedtools.BedTool(interval_list)
    raw_count = len(raw_intervals)

    if raw_count == 0:
        return ""

    # Expand every interval by `padding` bp on each side, then merge overlapping results.
    # This is the correct equivalent of GATK's IntervalListTools --PADDING: isolated intervals
    # get their outer boundaries expanded, not just the gaps between adjacent intervals.
    # Start is clamped to 0; end is not clamped since out-of-range BQ locations match nothing.
    if padding > 0:
        def _pad(interval):
            interval.start = max(0, interval.start - padding)
            interval.end = interval.end + padding
            return interval
        intervals = raw_intervals.each(_pad).saveas().merge()
    else:
        intervals = raw_intervals.merge()

    # Single pass: count merged intervals and sum bases covered (on known chromosomes only).
    interval_count = 0
    bases_covered = 0
    for i in intervals:
        if i.chrom not in CHROM_MAP:
            continue
        interval_count += 1
        bases_covered += int(i.end) - int(i.start)

    if interval_count < raw_count:
        print(f"Merged {raw_count:,} raw intervals into {interval_count:,} non-overlapping ranges "
              f"using {padding:,} bp padding.")
    else:
        print(f"Padding and merging with {padding:,} bp did not reduce the {interval_count:,} intervals.")

    pct_covered = bases_covered / HG38_GENOME_SIZE
    print(f"Merged intervals cover {pct_covered:.1%} of the genome "
          f"({bases_covered:,} / {HG38_GENOME_SIZE:,} bases, {interval_count:,} intervals).")

    if pct_covered >= skip_filter_threshold:
        print(f"Coverage {pct_covered:.1%} >= {skip_filter_threshold:.0%} threshold; "
              f"skipping location filtering.")
        return ""

    if interval_count <= INTERVAL_TEMP_TABLE_THRESHOLD:
        print(f"Using inline SQL filter for {interval_count:,} intervals.")
        location_clause_list = [
            f"""(location >= {CHROM_MAP[i.chrom]}{'0' * (12 - len(str(i.start)))}{i.start}
            AND location <= {CHROM_MAP[i.chrom]}{'0' * (12 - len(str(i.end)))}{i.end})"""
            for i in intervals
            if i.chrom in CHROM_MAP
        ]
        return "WHERE (" + " OR ".join(location_clause_list) + ")"
    else:
        print(f"Interval list has {interval_count:,} intervals (> {INTERVAL_TEMP_TABLE_THRESHOLD}); "
              f"loading into a temporary BigQuery table for efficient filtering.")
        fq_interval_table = load_intervals_to_temp_table(intervals, fq_temp_table_dataset, table_prefix)
        # BigQuery does not support WHERE EXISTS with a non-equality (range) correlated subquery —
        # it translates EXISTS to a LEFT SEMI JOIN and requires an equality condition.
        # An INNER JOIN with BETWEEN is equivalent here because the intervals are merged
        # (non-overlapping), so any location matches at most one interval row.
        # The equality condition on chrom_index lets BigQuery hash-partition the join by chromosome
        # before applying the BETWEEN range filter, avoiding a full nested loop across all intervals.
        return (f"INNER JOIN `{fq_interval_table}` "
                f"ON CAST(FLOOR(location / 1000000000000) AS INT64) = chrom_index "
                f"AND location BETWEEN location_start AND location_end")


def create_final_extract_vet_table(fq_destination_table_vet_data, enable_extract_table_ttl, vet_ranges_extract_table_version):
    ttl = ""
    if enable_extract_table_ttl:
        ttl = "OPTIONS( expiration_timestamp=TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL 14 DAY))"

    additional_columns = ""
    if vet_ranges_extract_table_version == "V2":
        additional_columns = f"""
            ,
            call_PGT      STRING,
            call_PID      STRING,
            call_PS       INT64
        """

    sql = f"""
        CREATE OR REPLACE TABLE `{fq_destination_table_vet_data}` 
        (
              location      INT64,
              sample_id	    INT64,
              ref           STRING,
              alt           STRING,
              call_GT       STRING,
              call_GQ       INT64,
              call_AD       STRING,
              AS_QUALapprox STRING,
              QUALapprox    STRING,
              call_PL       STRING
              {additional_columns}	
        )
          PARTITION BY RANGE_BUCKET(location, GENERATE_ARRAY(0, 26000000000000, 6500000000))
          CLUSTER BY location
          {ttl}        
        """
    print(sql)
    query_return = utils.execute_with_retry(client, "create final export vet table", sql)
    JOBS.append({'job': query_return['job'], 'label': query_return['label']})


def create_final_extract_ref_table(fq_destination_table_ref_data, enable_extract_table_ttl):
    ttl = ""
    if enable_extract_table_ttl:
        ttl = "OPTIONS( expiration_timestamp=TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL 14 DAY))"

    sql = f"""
        CREATE OR REPLACE TABLE `{fq_destination_table_ref_data}` 
        (
              location      INT64,
              sample_id	    INT64,
              length        INT64,
              state	        STRING	
        )
          PARTITION BY RANGE_BUCKET(location, GENERATE_ARRAY(0, 26000000000000, 6500000000))
          CLUSTER BY location
          {ttl}        
        """
    print(sql)
    query_return = utils.execute_with_retry(client, "create final export ref table", sql)
    JOBS.append({'job': query_return['job'], 'label': query_return['label']})

def populate_final_extract_table_with_ref(fq_ranges_dataset, fq_destination_table_data, sample_ids, use_compressed_references, location_string):
    # split file into files with x lines and then run
    def get_ref_subselect(fq_ref_table, samples, id):
        sample_stanza = ','.join([str(s) for s in samples])
        sql = f"    q_{id} AS (SELECT location, sample_id, length, state FROM \n" \
              f" `{fq_ref_table}` WHERE sample_id IN ({sample_stanza})), "
        return sql

    def get_compressed_ref_subselect(fq_ref_table, samples, id):
        sample_stanza = ','.join([str(s) for s in samples])
        sql = f"    q_{id} AS (SELECT UnpackRefRangeInfo(packed_ref_data).location as location, sample_id, UnpackRefRangeInfo(packed_ref_data).len as length, UnpackRefRangeInfo(packed_ref_data).state as state FROM \n" \
              f" `{fq_ref_table}` WHERE sample_id IN ({sample_stanza})), "
        return sql

    for i in range(1, REF_VET_TABLE_COUNT + 1):
        partition_samples = get_samples_for_partition(sample_ids, i)  # sample ids for the partition

        if len(partition_samples) > 0:
            subs = {}
            insert = f"\nINSERT INTO `{fq_destination_table_data}` (location, sample_id, length, state) \n WITH \n"
            fq_ref_table = f"{fq_ranges_dataset}.{REF_TABLE_PREFIX}{i:03}"
            j = 1

            for samples in split_lists(partition_samples, 1000):
                id = f"{i}_{j}"
                if use_compressed_references:
                    subs[id] = get_compressed_ref_subselect(fq_ref_table, samples, id)
                else:
                    subs[id] = get_ref_subselect(fq_ref_table, samples, id)
                j = j + 1

            helper_function_definitions = COMPRESSED_REFS_HELPER_FUNCTION_DEFS if use_compressed_references else ""

            sql = helper_function_definitions + insert + ("\n".join(subs.values())) + "\n" + \
                  "q_all AS (" + (" union all ".join([f"(SELECT * FROM q_{id})" for id in subs.keys()])) + ")\n" + \
                  f" (SELECT q_all.* FROM q_all {location_string})"
            print(sql)
            print(f"{fq_ref_table} query is {utils.utf8len(sql) / (1024 * 1024)} MB in length")
            query_return = utils.execute_with_retry(client, "populate destination table with reference data", sql)
            JOBS.append({'job': query_return['job'], 'label': query_return['label']})
    return


def populate_final_extract_table_with_vet(fq_ranges_dataset, fq_destination_table_data, sample_ids, vet_ranges_extract_table_version, location_string):
    additional_columns = ""
    if vet_ranges_extract_table_version == "V2":
        additional_columns = ", CALL_PGT, CALL_PID, CALL_PS"

        # split file into files with x lines and then run
    def get_ref_subselect(fq_vet_table, samples, id):
        sample_stanza = ','.join([str(s) for s in samples])
        sql = f"    q_{id} AS (SELECT location, sample_id, ref, alt, call_GT, call_GQ, call_AD, AS_QUALapprox, QUALapprox, CALL_PL{additional_columns} FROM \n" \
              f" `{fq_vet_table}` WHERE sample_id IN ({sample_stanza})), "
        return sql

    for i in range(1, REF_VET_TABLE_COUNT + 1):
        partition_samples = get_samples_for_partition(sample_ids, i)  # sample ids for the partition

        if len(partition_samples) > 0:
            subs = {}
            insert = f"\nINSERT INTO `{fq_destination_table_data}` (location, sample_id, ref, alt, call_GT, call_GQ, call_AD, AS_QUALapprox, QUALapprox, CALL_PL{additional_columns}) \n WITH \n"
            fq_vet_table = f"{fq_ranges_dataset}.{VET_TABLE_PREFIX}{i:03}"
            j = 1

            for samples in split_lists(partition_samples, 1000):
                id = f"{i}_{j}"
                subs[id] = get_ref_subselect(fq_vet_table, samples, id)
                j = j + 1

            sql = insert + ("\n".join(subs.values())) + "\n" + \
                  "q_all AS (" + (" union all ".join([f"(SELECT * FROM q_{id})" for id in subs.keys()])) + ")\n" + \
                  f" (SELECT q_all.* FROM q_all {location_string})"
            print(sql)
            print(f"{fq_vet_table} query is {utils.utf8len(sql) / (1024 * 1024)} MB in length")
            query_return = utils.execute_with_retry(client, "populate destination table with variant data", sql)
            JOBS.append({'job': query_return['job'], 'label': query_return['label']})
    return


def make_extract_table(call_set_identifier,
                       control_samples,
                       fq_ranges_dataset,
                       max_tables,
                       sample_names_to_extract,
                       fq_cohort_sample_names,
                       query_project,
                       query_labels,
                       fq_temp_table_dataset,
                       fq_destination_dataset,
                       destination_table_prefix,
                       fq_sample_mapping_table,
                       temp_table_ttl_hours,
                       only_output_vet_tables,
                       write_cost_to_db,
                       use_compressed_references,
                       vet_ranges_extract_table_version,
                       enable_extract_table_ttl,
                       interval_list,
                       interval_list_padding=1000,
                       skip_filter_threshold=0.5):
    try:
        fq_destination_table_ref_data = f"{fq_destination_dataset}.{destination_table_prefix}__REF_DATA"
        fq_destination_table_vet_data = f"{fq_destination_dataset}.{destination_table_prefix}__VET_DATA"
        fq_destination_table_samples = f"{fq_destination_dataset}.{destination_table_prefix}__SAMPLES"

        global client
        global default_config
        # this is where a set of labels are being created for the cohort extract
        query_labels_map = {
            "id": output_table_prefix,
            "gvs_tool_name": "gvs_prepare_ranges_callset"
        }

        # query_labels is string that looks like 'key1=val1, key2=val2'
        if query_labels is not None and len(query_labels) != 0:
            for query_label in query_labels:
                kv = query_label.split("=", 2)
                key = kv[0].strip().lower()
                value = kv[1].strip().lower()
                query_labels_map[key] = value

                if not (bool(re.match(r"[a-z0-9_-]+$", key)) & bool(re.match(r"[a-z0-9_-]+$", value))):
                    raise ValueError(
                        f"label key or value did not pass validation--format should be 'key1=val1, key2=val2'")

        # add labels for DSP Cloud Cost Control Labeling and Reporting
        query_labels_map.update({'service': 'gvs', 'team': 'variants', 'managedby': 'prepare_ranges_callset'})

        # Default QueryJobConfig will be merged into job configs passed in
        # but if a specific default config is being updated (eg labels), new config must be added
        # to the client._default_query_job_config that already exists
        default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE", use_query_cache=True, use_legacy_sql=False)

        client = bigquery.Client(project=query_project,
                                 default_query_job_config=default_config)

        # TODO -- provide a cmdline arg to override this (so we can simulate smaller datasets)
        global REF_VET_TABLE_COUNT
        REF_VET_TABLE_COUNT = max_tables

        global TEMP_TABLE_TTL_HOURS
        TEMP_TABLE_TTL_HOURS = temp_table_ttl_hours

        global TEMP_TABLE_TTL
        TEMP_TABLE_TTL = f" OPTIONS( expiration_timestamp=TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL {TEMP_TABLE_TTL_HOURS} HOUR)) "

        print(f"Using {REF_VET_TABLE_COUNT} tables in {fq_ranges_dataset}...")

        # Compute the location filter string once; for large interval lists this will also
        # create and populate a temporary BigQuery table with the encoded interval bounds.
        location_string = get_location_filters(interval_list, fq_temp_table_dataset, destination_table_prefix, interval_list_padding, skip_filter_threshold)

        # if we have a file of sample names, load it into a temporary table
        if sample_names_to_extract:
            fq_sample_name_table = load_sample_names(sample_names_to_extract, fq_temp_table_dataset)
        else:
            fq_sample_name_table = fq_cohort_sample_names

        # At this point one way or the other we have a table of sample names in BQ, join it to the sample_info table to
        # drive the extract. If this script was explicitly given a list of sample names then it should create the
        # cohort from those samples without regard to `withdrawn` on the `sample_info` table, otherwise only include
        # samples with a null `withdrawn` date in the cohort.
        if not only_output_vet_tables:
            create_extract_samples_table(control_samples, fq_destination_table_samples, fq_sample_name_table,
                                     fq_sample_mapping_table, not sample_names_to_extract, enable_extract_table_ttl)

        # pull the sample ids back down
        sample_ids = get_all_sample_ids(fq_destination_table_samples, only_output_vet_tables, fq_sample_mapping_table)

        # create and populate the tables for extract data
        if not only_output_vet_tables:
            create_final_extract_ref_table(fq_destination_table_ref_data, enable_extract_table_ttl)
            populate_final_extract_table_with_ref(fq_ranges_dataset, fq_destination_table_ref_data, sample_ids, use_compressed_references, location_string)

        create_final_extract_vet_table(fq_destination_table_vet_data, enable_extract_table_ttl, vet_ranges_extract_table_version)
        populate_final_extract_table_with_vet(fq_ranges_dataset, fq_destination_table_vet_data, sample_ids, vet_ranges_extract_table_version, location_string)

    finally:
        utils.write_job_stats(JOBS, client, f"{fq_destination_dataset}", call_set_identifier, 'GvsPrepareRanges',
                              'PrepareRangesCallsetTask', output_table_prefix, write_cost_to_db)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Extract a cohort from BigQuery Variant Store ')
    parser.add_argument('--call_set_identifier', type=str,
                        help='callset identifier used to track costs in cost_observability table', default='false')
    parser.add_argument('--control_samples', type=str,
                        help='true for control samples only, false for participant samples only', default='false')
    parser.add_argument('--fq_ranges_dataset', type=str, help='project.dataset location of ranges/vet data',
                        required=True)
    parser.add_argument('--fq_temp_table_dataset', type=str,
                        help='project.dataset location where results should be stored', required=True)
    parser.add_argument('--fq_destination_dataset', type=str,
                        help='project.dataset location where results should be stored', required=True)
    parser.add_argument('--destination_cohort_table_prefix', type=str,
                        help='prefix used for destination cohort extract tables (e.g. my_fantastic_cohort)',
                        required=True)
    parser.add_argument('--query_project', type=str, help='Google project where query should be executed',
                        required=True)
    parser.add_argument('--query_labels', type=str, action='append',
                        help='Labels to put on the BQ query that will show up in the billing. Ex: --query_labels key1=value1 --query_labels key2=value2',
                        required=False)
    parser.add_argument('--fq_sample_mapping_table', type=str,
                        help='Mapping table from sample_id to sample_name', required=True)
    parser.add_argument('--max_tables',type=int,
                        help='Maximum number of vet/ref ranges tables to consider', required=False, default=250)
    parser.add_argument('--ttl', type=int,
                        help='Temp table TTL in hours', required=False, default=72)
    parser.add_argument('--only_output_vet_tables', type=bool,
                        help='Only create __VET_DATA table, skip __REF_DATA and __SAMPLES tables', required=False, default=False)
    parser.add_argument('--write_cost_to_db', type=bool,
                        help='Populate cost_observability table with BigQuery query bytes scanned', required=False, default=True)
    parser.add_argument('--use_compressed_references', type=bool,
                        help='Expect compressed reference data and expand the fields', required=False, default=False)
    parser.add_argument('--vet-ranges-extract-table-version', type=str,
                       help='Version of the vet ranges extract table - for maintaining backwards-compatibility (V2 or V1)', required=False)
    parser.add_argument('--enable_extract_table_ttl', type=bool,
                        help='Add a TTL to the extract tables', required=False, default=False)
    parser.add_argument('--interval_list', type=str,
                        help='interval list or BAM file to limit the locations', required=False)
    parser.add_argument('--interval_list_padding', type=int,
                        help='Padding in base pairs applied to each side of every interval before merging. '
                             'Equivalent to merging intervals within a gap of 2 * padding bp. Default is 1000.',
                        required=False, default=1000)
    parser.add_argument('--skip_filter_coverage_threshold', type=float,
                        help='If the padded and merged intervals cover at least this fraction of the hg38 genome '
                             '(chromosomes 1-22, X, Y), location filtering is skipped and a full prepare is performed. '
                             'Default is 0.5 (50%%).',
                        required=False, default=0.5)

    sample_args = parser.add_mutually_exclusive_group(required=True)
    sample_args.add_argument('--sample_names_to_extract', type=str,
                             help='File containing list of samples to extract, 1 per line. ' +
                                  'All samples in this file will be included in the cohort regardless of `withdrawn` status in the `sample_info` table.')
    sample_args.add_argument('--fq_cohort_sample_names', type=str,
                             help='Fully qualified name of cohort table to extract, contains "sample_name" column. ' +
                                  'Only samples with null `withdrawn` fields in the `sample_info` table will be included in the cohort.')

    args = parser.parse_args()

    if not args.vet_ranges_extract_table_version:
        args.vet_ranges_extract_table_version = "V2"

    if args.vet_ranges_extract_table_version != "V1" and args.vet_ranges_extract_table_version != "V2":
        raise ValueError(f"Error: illegal value ({args.vet_ranges_extract_table_version}) for `vet-ranges-extract-table-version` must be either 'V1' or 'V2'")

    make_extract_table(args.call_set_identifier,
                       args.control_samples,
                       args.fq_ranges_dataset,
                       args.max_tables,
                       args.sample_names_to_extract,
                       args.fq_cohort_sample_names,
                       args.query_project,
                       args.query_labels,
                       args.fq_temp_table_dataset,
                       args.fq_destination_dataset,
                       args.destination_cohort_table_prefix,
                       args.fq_sample_mapping_table,
                       args.ttl,
                       args.only_output_vet_tables,
                       args.write_cost_to_db,
                       args.use_compressed_references,
                       args.vet_ranges_extract_table_version,
                       args.enable_extract_table_ttl,
                       args.interval_list,
                       args.interval_list_padding,
                       args.skip_filter_coverage_threshold)
