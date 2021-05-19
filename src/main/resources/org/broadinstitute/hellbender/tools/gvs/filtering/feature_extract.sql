WITH
    -- sum of REF AD for any sample where there is > 1 (ie >= 2) alternate alleles
    aa_ref_ad_info AS (
        SELECT location,
               IFNULL(sum(ref_ad),0) as ref_ad

        FROM (
               SELECT location, sample_id, MAX(ref_ad) ref_ad
               FROM `@altAllele`
               WHERE sample_id IN (SELECT sample_id FROM `@sample`)
               @locationStanza
               @trainingSitesStanza
               GROUP BY 1,2 HAVING SUM(ad) > 1
        )
        GROUP BY location),
    -- sum of ref sb_table info at the site level (and don't double count samples)
    aa_ref_sb_info AS (
        SELECT location,
               IFNULL(sum(sb_ref_plus),0) as sb_ref_plus,
               IFNULL(sum(sb_ref_minus),0) as sb_ref_minus
        FROM (
               SELECT DISTINCT location, sample_id, sb_ref_plus, sb_ref_minus
               FROM `@altAllele`
               WHERE sample_id IN (SELECT sample_id FROM `@sample`)
               @locationStanza
               @trainingSitesStanza
        )
        GROUP BY location),
     -- sum of qualapprox, count of het,homvar at the site level
   aa_site_info AS (
        SELECT location,
               IFNULL(sum(SAFE_CAST(qualapprox as FLOAT64)),0) as sum_qualapprox,
               COUNTIF(call_GT IN ('0/1', '0|1', '1/0', '1|0', '0/2', '0|2', '2/0', '2|0')) num_het_samples,
               COUNTIF(call_GT IN ('1/1', '1|1', '1/2', '1|2', '2/1', '2|1')) num_homvar_samples
        FROM (
               SELECT DISTINCT location, sample_id, qualapprox, call_GT
               FROM `@altAllele`
               WHERE sample_id IN (SELECT sample_id FROM `@sample`)
               @locationStanza
               @trainingSitesStanza
        )
        GROUP BY location),
   aa_site_allele_info AS (
        SELECT location,
               COUNT(DISTINCT ref || "," || allele) distinct_alleles,
               COUNTIF(length(ref) = length(allele)) num_snp_alleles
        FROM `@altAllele`
        WHERE sample_id IN (SELECT sample_id FROM `@sample`)
        @locationStanza
        @trainingSitesStanza
        GROUP BY location),
   hq_genotype_qc AS (
        SELECT location, COUNT(DISTINCT sample_id) hq_genotype_samples
        FROM `@altAllele`
        WHERE sample_id IN (SELECT sample_id FROM `@sample`)
        @locationStanza
        @trainingSitesStanza
        AND   call_GQ >= @hqGenotypeGQThreshold
        AND   (SELECT SUM(CAST(x AS INT64)) FROM UNNEST(SPLIT(call_AD,",")) x) >= @hqGenotypeDepthThreshold
        AND   ad / (SELECT SUM(CAST(x AS INT64)) FROM UNNEST(SPLIT(call_AD,",")) x) >= @hqGenotypeABThreshold
        GROUP BY location)
    SELECT
           ai.location,
           ai.ref,
           ai.allele,
           RAW_QUAL,
           aaradi.ref_ad as ref_ad,
           REF_AD_GT_1,
           AS_MQRankSum,
           AS_MQRankSum_ft,
           AS_ReadPosRankSum,
           AS_ReadPosRankSum_ft,
           RAW_MQ,
           RAW_AD,
           RAW_AD_GT_1,
           aarsbi.SB_REF_PLUS as SB_REF_PLUS,
           aarsbi.SB_REF_MINUS as SB_REF_MINUS,
           SB_ALT_PLUS,
           SB_ALT_MINUS,
           aasi.num_het_samples,
           aasi.num_homvar_samples,
           IFNULL(aasi.sum_qualapprox,0) as sum_qualapprox,
           IFNULL(aasai.distinct_alleles,0) distinct_alleles,
           IFNULL(aasai.num_snp_alleles,0) num_snp_alleles,
           IFNULL(hc_gt.hq_genotype_samples,0) hq_genotype_samples
    FROM (
    SELECT aa.location,
           ref,
           allele,
           IFNULL(SUM(qual),0) as RAW_QUAL,
           `bqutil`.fn.median(ARRAY_AGG( raw_mqranksum_x_10 IGNORE NULLS)) / 10.0 as AS_MQRankSum,
           freq_table(ARRAY_AGG(raw_mqranksum_x_10 IGNORE NULLS)) AS_MQRankSum_ft,
           `bqutil`.fn.median(ARRAY_AGG(raw_readposranksum_x_10 IGNORE NULLS)) / 10.0 as AS_ReadPosRankSum,
           freq_table(ARRAY_AGG(raw_readposranksum_x_10 IGNORE NULLS)) as AS_ReadPosRankSum_ft,
           IFNULL(SUM(RAW_MQ),0) as RAW_MQ,
           IFNULL(SUM(AD),0) as RAW_AD,
           IFNULL(SUM(CASE WHEN AD > 1 THEN AD ELSE 0 END),0) as RAW_AD_GT_1, # to match GATK implementation
           IFNULL(SUM(REF_AD),0) as REF_AD,
           IFNULL(SUM(CASE WHEN REF_AD > 1 THEN REF_AD ELSE 0 END),0) as REF_AD_GT_1, # to match GATK implementation
           IFNULL(SUM(SB_REF_PLUS),0)  as SB_REF_PLUS,
           IFNULL(SUM(SB_REF_MINUS),0) as SB_REF_MINUS,
           IFNULL(SUM(SB_ALT_PLUS),0)  as SB_ALT_PLUS,
           IFNULL(SUM(SB_ALT_MINUS),0) as SB_ALT_MINUS
    FROM `@altAllele` as aa
    WHERE allele != '*'
    AND sample_id in (SELECT sample_id from `@sample` )
    @locationStanza
    @trainingSitesStanza
    GROUP BY 1,2,3
    ) ai
    LEFT JOIN aa_ref_ad_info aaradi ON (ai.location = aaradi.location)
    LEFT JOIN aa_ref_sb_info aarsbi ON (ai.location = aarsbi.location)
    LEFT JOIN aa_site_info aasi ON (ai.location = aasi.location)
    LEFT JOIN aa_site_allele_info aasai ON (ai.location = aasai.location)
    LEFT JOIN hq_genotype_qc hc_gt ON (ai.location = hc_gt.location)

