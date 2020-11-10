package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;

public class ExtractFeaturesBQ {

    public static String getVQSRFeatureExtractQueryString(final TableReference vet, final TableReference altAllele, final TableReference sampleList,
                                                          final SimpleInterval interval, final boolean trainingSitesOnly) {

        //TODO try to remove the dependency on the vet tables since there are multiple tables that contain the data we need
        String trainingSitesStanza =
                !trainingSitesOnly?"":
                        "AND location IN (SELECT location FROM `broad-dsp-spec-ops.joint_genotyping_ref.vqsr_training_sites_*` WHERE chrom='chr20')\n";
        String query =
                "WITH ref_ad_info AS (\n" +
                        "SELECT \n" +
                        "  location,\n" +
                        "  IFNULL(sum(cast(SPLIT(call_AD,\",\")[OFFSET(0)] as int64)),0) as ref_ad\n" +
                        "FROM `@vet` \n" +
                        "WHERE sample_id IN (SELECT sample_id FROM `@sample`)\n" +
                        "AND (location >= @start AND location <= @end) \n" +
                        trainingSitesStanza +
                        "AND (\n" +
                        "  SELECT SUM(CAST(part AS int64)) FROM UNNEST(SPLIT(call_AD, ',')) part WITH OFFSET index WHERE index >= 1 ) \n" +
                        "  > 1\n" +
                        "GROUP BY location),\n" +
                        "ref_sb_info AS (\n" +
                        "SELECT \n" +
                        "  location,\n" +
                        "  IFNULL(sum(cast(SPLIT(SPLIT(as_sb_table,\"|\")[OFFSET(0)],\",\")[OFFSET(0)] as int64)),0) as sb_ref_plus, \n" +
                        "  IFNULL(sum(cast(SPLIT(SPLIT(as_sb_table,\"|\")[OFFSET(0)],\",\")[OFFSET(1)] as int64)),0) as sb_ref_minus  \n" +
                        "FROM `@vet` \n" +
                        "WHERE sample_id IN (SELECT sample_id FROM `@sample`)\n" +
                        "AND (location >= @start AND location <= @end) \n" +
                        trainingSitesStanza +
                        "GROUP BY location)\n" +
                        "SELECT \n" +
                        "       ai.location, \n" +
                        "       ai.ref, \n" +
                        "       ai.allele, \n" +
                        "       RAW_QUAL,\n" +
                        "       radi.ref_ad,\n" +
                        "       AS_MQRankSum,\n" +
                        "       AS_MQRankSum_ft,\n" +
                        "       AS_ReadPosRankSum,\n" +
                        "       AS_ReadPosRankSum_ft,\n" +
                        "       RAW_MQ,\n" +
                        "       RAW_AD,\n" +
                        "       RAW_AD_GT_1,\n" +
                        "       rsbi.SB_REF_PLUS,\n" +
                        "       rsbi.SB_REF_MINUS,\n" +
                        "       SB_ALT_PLUS, \n" +
                        "       SB_ALT_MINUS\n" +
                        "FROM (\n" +
                        "SELECT aa.location, \n" +
                        "       ref, \n" +
                        "       allele, \n" +
                        "       IFNULL(SUM(qual),0) as RAW_QUAL,\n" +
                        "       `bqutil`.fn.median(ARRAY_AGG( raw_mqranksum_x_10 IGNORE NULLS)) / 10.0 as AS_MQRankSum,\n" +
                        "       `broad-dsp-spec-ops`.joint_genotyping_ref.freq_table(ARRAY_AGG(raw_mqranksum_x_10 IGNORE NULLS)) AS_MQRankSum_ft,\n" +
                        "       `bqutil`.fn.median(ARRAY_AGG(raw_readposranksum_x_10 IGNORE NULLS)) / 10.0 as AS_ReadPosRankSum,\n" +
                        "       `broad-dsp-spec-ops`.joint_genotyping_ref.freq_table(ARRAY_AGG(raw_readposranksum_x_10 IGNORE NULLS)) as AS_ReadPosRankSum_ft,\n" +
                        "       IFNULL(SUM(RAW_MQ),0) as RAW_MQ,\n" +
                        "       IFNULL(SUM(AD),0) as RAW_AD, \n" +
                        "       IFNULL(SUM(CASE WHEN AD > 1 THEN AD ELSE 0 END),0) as RAW_AD_GT_1, # to match GATK implementation\n" +
                        "       IFNULL(SUM(SB_ALT_PLUS),0)  as SB_ALT_PLUS, \n" +
                        "       IFNULL(SUM(SB_ALT_MINUS),0) as SB_ALT_MINUS\n" +
                        "FROM `@altAllele` as aa\n" +
                        "WHERE (location >= @start AND location <= @end) \n" +
                        trainingSitesStanza +
                        "AND sample_id in (SELECT sample_id from `@sample` )\n" +
                        "AND allele != '*'\n" +
                        "GROUP BY 1,2,3\n" +
                        ") ai\n" +
                        "LEFT JOIN ref_ad_info radi ON (ai.location = radi.location)\n" +
                        "LEFT JOIN ref_sb_info rsbi ON (ai.location = rsbi.location)\n";

        return query
                .replaceAll("@vet", vet.getFQTableName())
                .replaceAll("@sample", sampleList.getFQTableName())
                .replaceAll("@start", String.format("%d", SchemaUtils.encodeLocation(interval.getContig(), interval.getStart())))
                .replaceAll( "@end", String.format("%d",SchemaUtils.encodeLocation(interval.getContig(),interval.getEnd())))
                .replaceAll( "@altAllele", altAllele.getFQTableName());
    }}
