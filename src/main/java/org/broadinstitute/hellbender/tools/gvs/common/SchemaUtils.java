package org.broadinstitute.hellbender.tools.gvs.common;

import org.apache.avro.generic.GenericRecord;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class SchemaUtils {

    public static final String LOCATION_FIELD_NAME = "location";

    public static final String SAMPLE_NAME_FIELD_NAME = "sample_name";
    public static final String SAMPLE_ID_FIELD_NAME = "sample_id";

    public static final String STATE_FIELD_NAME = "state";
    public static final String LENGTH_FIELD_NAME = "length";

    public static final String REF_ALLELE_FIELD_NAME = "ref";
    public static final String ALT_ALLELE_FIELD_NAME = "alt";

    public static final String GENOTYPE_FIELD_PREFIX = "call_";
    public static final String AS_FIELD_PREFIX = "AS_";
    public static final String MULTIVALUE_FIELD_DELIMITER = ",";

    public static final String CALL_GT = GENOTYPE_FIELD_PREFIX + "GT";
    public static final String CALL_GQ = GENOTYPE_FIELD_PREFIX + "GQ";
    public static final String CALL_RGQ = GENOTYPE_FIELD_PREFIX + "RGQ";

    // exomes & genomes
    public static final String RAW_QUAL = "RAW_QUAL";
    public static final String RAW_MQ = "RAW_MQ";
    public static final String AS_RAW_MQ = AS_FIELD_PREFIX + RAW_MQ;
    public static final String AS_MQRankSum = AS_FIELD_PREFIX + "MQRankSum";
    public static final String AS_RAW_MQRankSum = AS_FIELD_PREFIX + "RAW_MQRankSum";
    public static final String QUALapprox = "QUALapprox";
    public static final String AS_QUALapprox = AS_FIELD_PREFIX + "QUALapprox";
    public static final String AS_RAW_ReadPosRankSum = AS_FIELD_PREFIX + "RAW_ReadPosRankSum";
    public static final String AS_ReadPosRankSum = AS_FIELD_PREFIX + "ReadPosRankSum";
    public static final String AS_SB_TABLE = AS_FIELD_PREFIX + "SB_TABLE";
    public static final String AS_VarDP = AS_FIELD_PREFIX + "VarDP";
    public static final String CALL_AD = GENOTYPE_FIELD_PREFIX + "AD";
    public static final String SUM_AD = "SUM_AD";
    public static final String RAW_AD = "RAW_AD";
    public static final String CALL_PL = GENOTYPE_FIELD_PREFIX + "PL";

    //Filtering table
    public static final String FILTER_SET_NAME = "filter_set_name";
    public static final String VQSLOD = "vqslod";
    public static final String YNG_STATUS = "yng_status";

    //Tranches table
    public static final String TARGET_TRUTH_SENSITIVITY = "target_truth_sensitivity";
    public static final String MIN_VQSLOD = "min_vqslod";
    public static final String TRANCHE_FILTER_NAME = "filter_name";
    public static final String TRANCHE_MODEL = "model";

    // Site filtering table
    public static final String FILTERS = "filters";
    public static final List<String> FILTER_SET_SITE_FIELDS = Arrays.asList(FILTER_SET_NAME,LOCATION_FIELD_NAME,FILTERS);

    public static final List<String> EXTRACT_VET_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, SAMPLE_ID_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME, CALL_GT, CALL_GQ, AS_QUALapprox, QUALapprox, CALL_PL, CALL_AD);
    public static final List<String> EXTRACT_REF_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, SAMPLE_ID_FIELD_NAME, LENGTH_FIELD_NAME, STATE_FIELD_NAME);

    public static final List<String> SAMPLE_FIELDS = Arrays.asList(SchemaUtils.SAMPLE_NAME_FIELD_NAME, SchemaUtils.SAMPLE_ID_FIELD_NAME);
    public static final List<String> YNG_FIELDS = Arrays.asList(FILTER_SET_NAME, LOCATION_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME, VQSLOD, YNG_STATUS);
    public static final List<String> TRANCHE_FIELDS = Arrays.asList(TARGET_TRUTH_SENSITIVITY, MIN_VQSLOD, TRANCHE_FILTER_NAME, TRANCHE_MODEL);

    public static final List<String> ALT_ALLELE_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, SAMPLE_ID_FIELD_NAME, REF_ALLELE_FIELD_NAME, "allele", ALT_ALLELE_FIELD_NAME, "allele_pos", CALL_GT, AS_RAW_MQ, RAW_MQ, AS_RAW_MQRankSum, "raw_mqranksum_x_10", AS_QUALapprox, "qual", AS_RAW_ReadPosRankSum, "raw_readposranksum_x_10", AS_SB_TABLE, "SB_REF_PLUS","SB_REF_MINUS","SB_ALT_PLUS","SB_ALT_MINUS", CALL_AD, "ref_ad", "ad");
    public static final List<String> FEATURE_EXTRACT_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, REF_ALLELE_FIELD_NAME, "allele", RAW_QUAL, "ref_ad", AS_MQRankSum, "AS_MQRankSum_ft", AS_ReadPosRankSum, "AS_ReadPosRankSum_ft", RAW_MQ, SUM_AD, RAW_AD, "RAW_AD_GT_1", "SB_REF_PLUS","SB_REF_MINUS","SB_ALT_PLUS","SB_ALT_MINUS","num_het_samples","num_homvar_samples","distinct_alleles","hq_genotype_samples", "sum_qualapprox", "num_snp_alleles");

    public static final String LOAD_STATUS_FIELD_NAME = "status";
    public static final String LOAD_STATUS_EVENT_TIMESTAMP_NAME = "event_timestamp";

    public static final long chromAdjustment = 1000000000000L;

    public static long encodeLocation(String chrom, int position) {
        int chromosomeIndex = ChromosomeEnum.valueOfContig(chrom).index;
        return (long) chromosomeIndex * chromAdjustment + (long) position;
    }

    public static String decodeContig(long location) {
        return ChromosomeEnum.valueOfIndex((int)(location/chromAdjustment)).getContigName();
    }

    public static int decodePosition(long location) {
        return (int)(location % chromAdjustment);
    }

    public final static Comparator<GenericRecord> LOCATION_COMPARATOR = (o1, o2) -> {
        final long firstLocation = (Long) o1.get(SchemaUtils.LOCATION_FIELD_NAME);
        final long secondLocation = (Long) o2.get(SchemaUtils.LOCATION_FIELD_NAME);
        return Long.compare(firstLocation, secondLocation);
    };

}
