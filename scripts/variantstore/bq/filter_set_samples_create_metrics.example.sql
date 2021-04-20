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
CREATE OR REPLACE TABLE `aou-genomics-curation-prod.alpha1_1000.filter_set_samples` AS    
SELECT "kc-testing-1" filter_set_name,
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
        FROM `aou-genomics-curation-prod.alpha1_1000.vet_001` v
        LEFT JOIN `spec-ops-aou.gvs_public_reference_data.gnomad_v3_sites` gnomad ON (v.location = gnomad.location)
        WHERE call_GT != "./." -- for safety
        AND v.location < 23000000000000 -- autosomal only
) GROUP BY 1,2
