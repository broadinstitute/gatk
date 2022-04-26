CREATE TEMPORARY FUNCTION titv(ref STRING, allele STRING)
RETURNS STRING
  LANGUAGE js AS """
     if ( ref.length > 1 || allele.length > 1) {
         return "other";
     } else if ( (ref == "A" && allele == "G") ||
                 (ref == "G" && allele == "A") ||
                 (ref == "C" && allele == "T") ||
                 (ref == "T" && allele == "C") ) {
          return "ti";
     } else {
          return "tv";
     }
    """;

CREATE TEMPORARY FUNCTION type(ref STRING, allele STRING, gt_str STRING)
RETURNS STRING
  LANGUAGE js AS """

alts = allele.split(",")

// get the the non-reference allele indexes
ai = gt_str.replace("|","/").split("/").filter(i => i != "0");

// the the distinct set of lengths of the alternates
alt_lengths = new Set(ai.map(i => alts[parseInt(i)-1].length))

if (alt_lengths.size > 1) {
  return "complex"
} else {
  // get first (only) element
  al = alt_lengths.keys().next().value

  if ( ref.length == al && al == 1) {
         return "snp"
     } else if (ref.length > al) {
         return "del"
     } else if (ref.length < al) {
         return "ins"
     } else {
         return "other"
     }
}
    """;

-- TODO Do we eventually want to support multiple filter sets in this table?
-- If so have separate 'create table if not exists' statement and make this an insert
-- Note that as the data gets larger and larger we will need to split up this query.
--   - Splitting by location is one simple way to do this so that the query only scans parts of the genome at a time.
--   - The agg / group by still needs to get run though outside of the location splitting since it's by sample across all positions.
--   - So likely this will need to be 2 groupings of queries----one group of 2+ queries that scans the VETs by location, and one that aggregates all of that data together by sample.

CREATE OR REPLACE TABLE `$FQ_PREFIX_sample_metrics` AS
SELECT "$NAME_OF_FILTER_SET" filter_set_name,
       sample_id,
       count(1) variant_entries,
       SUM(CASE WHEN type = "del" THEN 1 ELSE 0 END) del_count,
       SUM(CASE WHEN type = "ins" THEN 1 ELSE 0 END) ins_count,
       SUM(CASE WHEN type = "snp" THEN 1 ELSE 0 END) snp_count,
       SUM(CASE WHEN type = "snp" AND titv = "ti" THEN 1 ELSE 0 END) ti_count, # TODO: minimize alleles
       SUM(CASE WHEN type = "snp" AND titv = "tv" THEN 1 ELSE 0 END) tv_count, # TODO: minimize alleles
       SUM(CASE WHEN type = "snp" AND gt_type = "het" THEN 1 ELSE 0 END) snp_het_count,
       SUM(CASE WHEN type = "snp" AND gt_type = "homvar" THEN 1 ELSE 0 END) snp_homvar_count,
       SUM(CASE WHEN type IN ("ins","del") AND gt_type = "het" THEN 1 ELSE 0 END) indel_het_count,
       SUM(CASE WHEN type IN ("ins","del") AND gt_type = "homvar" THEN 1 ELSE 0 END) indel_homvar_count,
       COUNTIF(not in_gnomad) singleton,
       null AS pass_qc  -- placeholder
    FROM (
        SELECT sample_id,
               type(ref, alt, call_GT) as type,
               CASE WHEN INSTR(call_GT, "0") > 0 THEN "het" ELSE "homvar" END as gt_type, -- if GT contains a zero, its a het
               titv(ref, alt) as titv,
               CASE WHEN gnomad.location IS NULL THEN false ELSE true END in_gnomad
        FROM `$FQ_PREFIX__VET_DATA` v
        LEFT JOIN `spec-ops-aou.gvs_public_reference_data.gnomad_v3_sites` gnomad ON (v.location = gnomad.location)
        WHERE call_GT != "./." -- for safety
        AND v.location < 23000000000000 -- autosomal only
) GROUP BY 1,2


-- BELOW are rori's not so great queries that did their job for Charlie, but are not ideal
-- part of why they are not ideal is there has to be a better way of doing the group by---this is sloppy


-- Group 1
-- create table & do ~ first half

CREATE OR REPLACE TABLE `$FQ_PREFIX_sample_metrics` AS
SELECT "{NAME_OF_FILTER_SET}" filter_set_name,
       sample_id,
       count(1) variant_entries,
       SUM(CASE WHEN type = "del" THEN 1 ELSE 0 END) del_count,
       SUM(CASE WHEN type = "ins" THEN 1 ELSE 0 END) ins_count,
       SUM(CASE WHEN type = "snp" THEN 1 ELSE 0 END) snp_count,
       SUM(CASE WHEN type = "snp" AND titv = "ti" THEN 1 ELSE 0 END) ti_count, # TODO: minimize alleles
       SUM(CASE WHEN type = "snp" AND titv = "tv" THEN 1 ELSE 0 END) tv_count, # TODO: minimize alleles
       SUM(CASE WHEN type = "snp" AND gt_type = "het" THEN 1 ELSE 0 END) snp_het_count,
       SUM(CASE WHEN type = "snp" AND gt_type = "homvar" THEN 1 ELSE 0 END) snp_homvar_count,
       SUM(CASE WHEN type IN ("ins","del") AND gt_type = "het" THEN 1 ELSE 0 END) indel_het_count,
       SUM(CASE WHEN type IN ("ins","del") AND gt_type = "homvar" THEN 1 ELSE 0 END) indel_homvar_count,
       COUNTIF(not in_gnomad) singleton,
       null AS pass_qc
    FROM (
        SELECT sample_id,
               type(ref, alt, call_GT) as type,
               CASE WHEN INSTR(call_GT, "0") > 0 THEN "het" ELSE "homvar" END as gt_type,
               titv(ref, alt) as titv,
               CASE WHEN gnomad.location IS NULL THEN false ELSE true END in_gnomad
        FROM `{FQ_PREFIX}__VET_DATA` v
        LEFT JOIN `spec-ops-aou.gvs_public_reference_data.gnomad_v3_sites` gnomad ON (v.location = gnomad.location)
        WHERE call_GT != "./." -- for safety
        AND v.location < 10000000000000 -- split by location, ~ first half
) GROUP BY 1,2

-- insert ~ second half into table

INSERT `{FQ_PREFIX}_sample_metrics` (
filter_set_name,
sample_id,
variant_entries,
del_count,
ins_count,
snp_count,
ti_count,
tv_count,
snp_het_count,
snp_homvar_count,
indel_het_count,
indel_homvar_count,
singleton,
pass_qc
)
SELECT "{NAME_OF_FILTER_SET}" filter_set_name,
       sample_id,
       count(1) variant_entries,
       SUM(CASE WHEN type = "del" THEN 1 ELSE 0 END) del_count,
       SUM(CASE WHEN type = "ins" THEN 1 ELSE 0 END) ins_count,
       SUM(CASE WHEN type = "snp" THEN 1 ELSE 0 END) snp_count,
       SUM(CASE WHEN type = "snp" AND titv = "ti" THEN 1 ELSE 0 END) ti_count, # TODO: minimize alleles
    SUM(CASE WHEN type = "snp" AND titv = "tv" THEN 1 ELSE 0 END) tv_count, # TODO: minimize alleles
       SUM(CASE WHEN type = "snp" AND gt_type = "het" THEN 1 ELSE 0 END) snp_het_count,
       SUM(CASE WHEN type = "snp" AND gt_type = "homvar" THEN 1 ELSE 0 END) snp_homvar_count,
       SUM(CASE WHEN type IN ("ins","del") AND gt_type = "het" THEN 1 ELSE 0 END) indel_het_count,
       SUM(CASE WHEN type IN ("ins","del") AND gt_type = "homvar" THEN 1 ELSE 0 END) indel_homvar_count,
       COUNTIF(not in_gnomad) singleton,
       null AS pass_qc
FROM (
    SELECT sample_id,
    type(ref, alt, call_GT) as type,
    CASE WHEN INSTR(call_GT, "0") > 0 THEN "het" ELSE "homvar" END as gt_type,
    titv(ref, alt) as titv,
    CASE WHEN gnomad.location IS NULL THEN false ELSE true END in_gnomad
    FROM `{FQ_PREFIX}__VET_DATA` v
    LEFT JOIN `spec-ops-aou.gvs_public_reference_data.gnomad_v3_sites` gnomad ON (v.location = gnomad.location)
    WHERE call_GT != "./."
    AND v.location >= 10000000000000  -- split by location, ~ second half
    AND v.location < 23000000000000 -- autosomal only
) GROUP BY 1,2


-- Group 2
-- create aggregated table

CREATE OR REPLACE TABLE `{FQ_PREFIX}_sample_metrics_agg` AS -- note that now a different table name must be used in the next query
SELECT "{NAME_OF_FILTER_SET}" filter_set_name,
       sample_id,
       SUM(variant_entries) variant_entries,
       SUM(CASE WHEN type = "del" THEN 1 ELSE 0 END) del_count,
       SUM(CASE WHEN type = "ins" THEN 1 ELSE 0 END) ins_count,
       SUM(CASE WHEN type = "snp" THEN 1 ELSE 0 END) snp_count,
       SUM(CASE WHEN type = "snp" AND titv = "ti" THEN 1 ELSE 0 END) ti_count, # TODO: minimize alleles
    SUM(CASE WHEN type = "snp" AND titv = "tv" THEN 1 ELSE 0 END) tv_count, # TODO: minimize alleles
       SUM(CASE WHEN type = "snp" AND gt_type = "het" THEN 1 ELSE 0 END) snp_het_count,
       SUM(CASE WHEN type = "snp" AND gt_type = "homvar" THEN 1 ELSE 0 END) snp_homvar_count,
       SUM(CASE WHEN type IN ("ins","del") AND gt_type = "het" THEN 1 ELSE 0 END) indel_het_count,
       SUM(CASE WHEN type IN ("ins","del") AND gt_type = "homvar" THEN 1 ELSE 0 END) indel_homvar_count,
       SUM(singleton) singleton,
       null AS pass_qc
FROM `{FQ_PREFIX}_sample_metrics`
GROUP BY 1,2


