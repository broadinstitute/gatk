-- Sample name to sample ID map for VDS dropout screening (VS-1998).
--
-- YOU PROBABLY DO NOT NEED TO RUN THIS. GvsValidateVdsCompleteness.wdl generates the map
-- itself when given bq_project_id and bq_dataset_name. This file is the hand-run copy, for
-- building a map
-- outside the workflow or inspecting the sample universe before committing to a run. The
-- workflow's GenerateSampleMap task inlines the same projection and filter; a unit test
-- asserts the two agree, because a divergent sample universe changes what gets screened
-- without erroring.
--
-- Produces the `--sample-map-path` input to vds_dropout_scan.py, which needs it because a
-- VDS knows only sample names while superpartition membership is a function of sample_id.
--
-- The filter matches what GvsExtractAvroFilesForHail.wdl applies when exporting Avro, so
-- the sample universe here is the same one the VDS was built from. Samples withdrawn since
-- the VDS was written are therefore absent from this map even though the VDS still holds
-- them; vds_dropout_scan.py treats a VDS sample missing from the map as a hard error rather
-- than dropping it silently, because a silently excluded sample would bias the
-- superpartition totals the whole comparison rests on. If that fires, regenerate the map
-- without the `withdrawn` predicate and re-run.
--
-- Substitute PROJECT_ID and DATASET_NAME before running, e.g.
--
--   sed -e 's/PROJECT_ID/aou-genomics-curation-prod/' -e 's/DATASET_NAME/foxtrot/' \
--       vds_dropout_sample_map.sql > /tmp/sample_map.sql
--   bq query --project_id=aou-genomics-curation-prod --use_legacy_sql=false \
--       --format=csv --max_rows=10000000 < /tmp/sample_map.sql > sample_map.tsv
--   gsutil cp sample_map.tsv gs://<bucket>/vs-1998/sample_map.tsv
--
-- vds_dropout_scan.py accepts either a tab- or comma-separated map, so `--format=csv`
-- output can be used as-is. For a callset large enough that `bq query` is unwieldy, see
-- the EXPORT DATA variant at the bottom of this file, which writes straight to GCS.

SELECT
  sample_name,
  sample_id
FROM `PROJECT_ID.DATASET_NAME.sample_info`
WHERE withdrawn IS NULL
  AND is_control = false
ORDER BY sample_id;

-- ---------------------------------------------------------------------------
-- Alternative: write directly to GCS, tab-delimited, no local round trip.
--
-- EXPORT DATA requires exactly one wildcard in the uri, so the result lands at
-- .../sample_map_000000000000.tsv (one file at AoU scale). Pass that full path to
-- --sample-map-path.
-- ---------------------------------------------------------------------------
--
-- EXPORT DATA OPTIONS(
--   uri='gs://BUCKET/vs-1998/sample_map_*.tsv',
--   format='CSV',
--   field_delimiter='\t',
--   header=true,
--   overwrite=true) AS
-- SELECT
--   sample_name,
--   sample_id
-- FROM `PROJECT_ID.DATASET_NAME.sample_info`
-- WHERE withdrawn IS NULL
--   AND is_control = false
-- ORDER BY sample_id;
