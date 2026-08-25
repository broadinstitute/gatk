-- Stratified sample selection for VDS dropout screening (VS-1998).
--
-- YOU PROBABLY DO NOT NEED TO RUN THIS. `vds_dropout_scan.py --action materialize` already
-- writes the sample list as a side effect, and that copy is the authoritative one because
-- it is drawn from the samples actually present in the VDS. The normal sequence is:
--
--   1. run vds_dropout_sample_map.sql            -> sample_map.tsv
--   2. --action probe       (source VDS, --sample-map-path)   -> cost estimate
--   3. --action materialize (source VDS, --sample-map-path)   -> narrow VDS + sample_list.tsv
--   4. --action scan        (narrow VDS, --sample-list-path from step 3)
--
-- This file is for the cases below, none of which are on that path.
--
-- Reproduces `stratified_sample` from vds_dropout_scan.py in SQL: for each GVS
-- superpartition, take the SAMPLES_PER_SUPERPARTITION samples whose hash sorts first.
--
-- Why it exists at all
-- --------------------
-- Three reasons, the third being the real one:
--
--   1. Inspect the selection before spending any cluster time -- how many samples land in
--      each superpartition, which superpartitions are short.
--   2. Fix the sample set without a cluster, so the set the r3 comparison will use can be
--      agreed before there is an r3 VDS to run against.
--   3. Cross-check the Python. The recorded sample list is the thing that makes figures
--      comparable across r2 and r3, so "the selection is reproducible" is load-bearing
--      rather than a nicety. Two independent implementations agreeing is real evidence;
--      one implementation agreeing with itself is not.
--
-- The output format and row order are deliberately identical to what
-- `vds_dropout_scan.py --action materialize` writes, so the two can be compared directly:
--
--   diff <(tail -n +2 sample_list_from_sql.tsv) <(tail -n +2 sample_list_from_hail.tsv)
--
-- Equivalence argument
-- --------------------
-- The Python sorts by `sha256(f'{seed}:{name}').hexdigest()` then by name; this sorts by
-- `TO_HEX(SHA256(CONCAT(seed, ':', sample_name)))` then by sample_name. Those agree
-- because:
--
--   * BigQuery SHA256 hashes the UTF-8 bytes of a STRING, as Python's `.encode('utf-8')`
--     does, so the digests are the same bytes.
--   * TO_HEX and Python's hexdigest() both emit lowercase, fixed-width 64-character hex,
--     so lexicographic string ordering is the same total order in both.
--   * Both tie-break on the sample name, so neither depends on input ordering. The Python
--     used to rely on sort stability here, which would have been a real divergence.
--
-- Superpartition membership is `CAST(CEIL(sample_id / 4000.0) AS INT64)`, matching
-- GvsExtractAvroFilesForHail.wdl so that the numbering lines up with the vet_NNN and
-- ref_ranges_NNN table names. 4000 is the BigQuery partition limit, not a tunable, so it
-- is written inline rather than declared.
--
-- One difference that is not a bug: this selects from sample_info, which does not know
-- what any particular VDS contains. `--action materialize` selects from samples actually
-- present in the VDS. So a list generated here can name samples a given VDS lacks;
-- vds_dropout_scan.py intersects and reports the shortfall, and the per-superpartition
-- sample counts it writes reflect what was actually screened. For a single authoritative
-- list, prefer the one materialize writes.
--
-- Substitute PROJECT_ID and DATASET_NAME before running, e.g.
--
--   sed -e 's/PROJECT_ID/aou-genomics-curation-prod/' -e 's/DATASET_NAME/foxtrot/' \
--       -e 's|BUCKET|fc-secure-...|' \
--       vds_dropout_sample_list.sql > /tmp/sample_list.sql
--   bq query --project_id=aou-genomics-curation-prod --use_legacy_sql=false \
--       --format=csv --max_rows=10000000 < /tmp/sample_list.sql > sample_list.tsv
--
-- There are two output paths at the bottom of the file: returning the rows (the default,
-- for inspecting and diffing) and EXPORT DATA straight to GCS (for feeding the WDL). The
-- selection itself is built once into a temp table, so the two cannot disagree.
--
-- Keep SEED and SAMPLES_PER_SUPERPARTITION in step with the defaults in
-- vds_dropout_scan.py (DEFAULT_SEED, DEFAULT_SAMPLES_PER_SUPERPARTITION); a unit test
-- asserts these literals match, because silent drift here would break exactly the
-- cross-VDS comparability the list exists to guarantee.

DECLARE SEED STRING DEFAULT 'vs-1998';
DECLARE SAMPLES_PER_SUPERPARTITION INT64 DEFAULT 100;

-- The selection is defined once, here, so the two output paths below cannot drift apart.
CREATE TEMP TABLE stratified_selection AS
SELECT
  sample_name,
  sample_id,
  superpartition
FROM (
  SELECT
    sample_name,
    sample_id,
    CAST(CEIL(sample_id / 4000.0) AS INT64) AS superpartition,
    ROW_NUMBER() OVER (
      PARTITION BY CAST(CEIL(sample_id / 4000.0) AS INT64)
      ORDER BY TO_HEX(SHA256(CONCAT(SEED, ':', sample_name))), sample_name
    ) AS selection_rank
  FROM `PROJECT_ID.DATASET_NAME.sample_info`
  WHERE withdrawn IS NULL
    AND is_control = false
)
WHERE selection_rank <= SAMPLES_PER_SUPERPARTITION;

-- Path 1: return the rows, for eyeballing the selection or diffing it against the list
-- that `--action materialize` writes.
SELECT
  sample_name,
  sample_id,
  superpartition
FROM stratified_selection
ORDER BY superpartition, sample_name;

-- ---------------------------------------------------------------------------
-- Path 2: write straight to GCS, tab-delimited with a header, ready to pass to
-- --sample-list-path without a local round trip. Comment out the SELECT above and
-- uncomment this block; the DECLAREs and the temp table above serve either path.
--
-- EXPORT DATA requires exactly one wildcard in the uri, so the result lands at
-- .../sample_list_000000000000.tsv. At a sampling depth of 100 across ~134
-- superpartitions this is ~13,400 rows, comfortably one file; pass that full path.
--
-- ORDER BY is preserved across exported files, so the diff against materialize's output
-- still holds even in the multi-file case -- the same property the Avro export relies on
-- to give each shard a contiguous location range.
-- ---------------------------------------------------------------------------
--
-- EXPORT DATA OPTIONS(
--   uri='gs://BUCKET/vs-1998/sample_list_*.tsv',
--   format='CSV',
--   field_delimiter='\t',
--   header=true,
--   overwrite=true) AS
-- SELECT
--   sample_name,
--   sample_id,
--   superpartition
-- FROM stratified_selection
-- ORDER BY superpartition, sample_name;

-- ---------------------------------------------------------------------------
-- Companion: superpartition sizes and how many samples each contributes.
--
-- Useful before a run to see which superpartitions are short of the sampling depth -- the
-- last one in a callset always is -- and to confirm the superpartition count matches
-- expectations (~134 for a 535K-sample Foxtrot).
-- ---------------------------------------------------------------------------
--
-- SELECT
--   superpartition,
--   COUNT(*) AS samples_to_be_screened
-- FROM stratified_selection
-- GROUP BY superpartition
-- ORDER BY superpartition;
