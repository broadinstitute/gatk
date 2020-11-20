package org.broadinstitute.hellbender.tools.variantdb;

import com.google.common.collect.ImmutableSet;
import org.apache.avro.generic.GenericRecord;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class SchemaUtils {

//    public static final String POSITION_FIELD_NAME = "position";
//    public static final String CHROM_FIELD_NAME = "chrom";
    public static final String LOCATION_FIELD_NAME = "location";

    public static final String SAMPLE_NAME_FIELD_NAME = "sample_name";
    public static final String SAMPLE_ID_FIELD_NAME = "sample_id";

    public static final String BASIC_ARRAY_DATA_FIELD_NAME = "basic_array_data";
    public static final String RAW_ARRAY_DATA_FIELD_NAME = "raw_array_data";

    // TODO remove this one - we should not have this ambiguous field
//    public static final String SAMPLE_FIELD_NAME = "sample";
    public static final String STATE_FIELD_NAME = "state";
    public static final String REF_ALLELE_FIELD_NAME = "ref";
    public static final String ALT_ALLELE_FIELD_NAME = "alt";

    public static final String GENOTYPE_FIELD_PREFIX = "call_";
    public static final String AS_FIELD_PREFIX = "AS_";
    public static final String MULTIVALUE_FIELD_DELIMITER = ",";

    public static final String CALL_GT = GENOTYPE_FIELD_PREFIX + "GT";
    public static final String CALL_GQ = GENOTYPE_FIELD_PREFIX + "GQ";
    public static final String CALL_RGQ = GENOTYPE_FIELD_PREFIX + "RGQ";

    public static final ImmutableSet<String> REQUIRED_FIELDS = ImmutableSet.of(
            SAMPLE_NAME_FIELD_NAME,
//            SAMPLE_FIELD_NAME,
            STATE_FIELD_NAME,
            REF_ALLELE_FIELD_NAME,
            ALT_ALLELE_FIELD_NAME
    );

    public static final List<String> COHORT_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, SAMPLE_NAME_FIELD_NAME, STATE_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME, CALL_GT, CALL_GQ, CALL_RGQ);
    public static final List<String> ARRAY_COHORT_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, SAMPLE_NAME_FIELD_NAME, STATE_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME, CALL_GT, CALL_GQ);

    public static final List<String> RAW_ARRAY_COHORT_FIELDS_COMPRESSED = 
        Arrays.asList(BASIC_ARRAY_DATA_FIELD_NAME, RAW_ARRAY_DATA_FIELD_NAME);
        
    public static final List<String> RAW_ARRAY_COHORT_FIELDS_UNCOMPRESSED = 
        Arrays.asList(SAMPLE_ID_FIELD_NAME, "probe_id", "GT_encoded","NORMX","NORMY","BAF","LRR");

    // exomes & genomes
    public static final String RAW_QUAL = "RAW_QUAL";
    public static final String RAW_MQ = "RAW_MQ";
    public static final String AS_RAW_MQ = AS_FIELD_PREFIX + RAW_MQ;
    public static final String AS_MQRankSum = AS_FIELD_PREFIX + "MQRankSum";
    public static final String AS_RAW_MQRankSum = AS_FIELD_PREFIX + "RAW_MQRankSum";
    public static final String AS_QUALapprox = AS_FIELD_PREFIX + "QUALapprox";
    public static final String AS_RAW_ReadPosRankSum = AS_FIELD_PREFIX + "RAW_ReadPosRankSum";
    public static final String AS_ReadPosRankSum = AS_FIELD_PREFIX + "ReadPosRankSum";
    public static final String AS_SB_TABLE = AS_FIELD_PREFIX + "SB_TABLE";
    public static final String AS_VarDP = AS_FIELD_PREFIX + "VarDP";
    public static final String CALL_AD = GENOTYPE_FIELD_PREFIX + "AD";
    public static final String RAW_AD = "RAW_AD";
    public static final String CALL_PGT = GENOTYPE_FIELD_PREFIX + "PGT";
    public static final String CALL_PID = GENOTYPE_FIELD_PREFIX + "PID";
    public static final String CALL_PL = GENOTYPE_FIELD_PREFIX + "PL";

    //Filtering table
    public static final String FILTER_SET_NAME = "filter_set_name";
    public static final String TYPE = "type";
    public static final String VQSLOD = "vqslod";
    public static final String CULPRIT = "culprit";
    public static final String TRAINING_LABEL = "training_label";
    public static final String YNG_STATUS = "yng_status";

    public static final List<String> SAMPLE_FIELDS = Arrays.asList(SchemaUtils.SAMPLE_NAME_FIELD_NAME, SchemaUtils.SAMPLE_ID_FIELD_NAME);
    public static final List<String> YNG_FIELDS = Arrays.asList(FILTER_SET_NAME, LOCATION_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME, VQSLOD, YNG_STATUS);


    public static final List<String> PET_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, SAMPLE_ID_FIELD_NAME, STATE_FIELD_NAME);
    public static final List<String> VET_FIELDS = Arrays.asList(SAMPLE_ID_FIELD_NAME, LOCATION_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME, AS_RAW_MQ,
            AS_RAW_MQRankSum, AS_QUALapprox, AS_RAW_ReadPosRankSum, AS_SB_TABLE, AS_VarDP, CALL_GT, CALL_AD, CALL_GQ, CALL_PGT, CALL_PID, CALL_PL);
    public static final List<String> ALT_ALLELE_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, SAMPLE_ID_FIELD_NAME, REF_ALLELE_FIELD_NAME, "allele", ALT_ALLELE_FIELD_NAME, "allele_pos", CALL_GT, AS_RAW_MQ, RAW_MQ, AS_RAW_MQRankSum, "raw_mqranksum_x_10", AS_QUALapprox, "qual", AS_RAW_ReadPosRankSum, "raw_readposranksum_x_10", AS_SB_TABLE, "SB_REF_PLUS","SB_REF_MINUS","SB_ALT_PLUS","SB_ALT_MINUS", CALL_AD, "ref_ad", "ad");
    public static final List<String> FEATURE_EXTRACT_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, REF_ALLELE_FIELD_NAME, "allele", RAW_QUAL, "ref_ad", AS_MQRankSum, "AS_MQRankSum_ft", AS_ReadPosRankSum, "AS_ReadPosRankSum_ft", RAW_MQ, RAW_AD, "RAW_AD_GT_1", "SB_REF_PLUS","SB_REF_MINUS","SB_ALT_PLUS","SB_ALT_MINUS");


//    private String getPositionTableForContig(final String contig ) {
//        return contigToPositionTableMap.get(contig);
//    }
//
//    private String getVariantTableForContig(final String contig ) {
//        return contigToVariantTableMap.get(contig);
//    }
//

//    private void validateSchema(final Set<String> columnNames) {
//        for ( final String requiredField : REQUIRED_FIELDS ) {
//            if ( ! columnNames.contains(requiredField) ) {
//                throw new UserException("Missing required column: " + requiredField +
//                        ". Actual columns encountered were: " + columnNames);
//            }
//        }
//    }


    public static final long chromAdjustment = 1000000000000L;

    public static long encodeLocation(String chrom, int position) {
        int chromosomeIndex = ChromosomeEnum.valueOfContig(chrom).index;
        long adjustedLocation = Long.valueOf(chromosomeIndex) * chromAdjustment + Long.valueOf(position);
        return adjustedLocation;
    }

    public static String decodeContig(long location) {
        return ChromosomeEnum.valueOfIndex((int)(location/chromAdjustment)).getContigName();
    }

    public static int decodePosition(long location) {
        return (int)(location % chromAdjustment);
    }

    public final static Comparator<GenericRecord> LOCATION_COMPARATOR = new Comparator<GenericRecord>() {
        @Override
        public int compare( GenericRecord o1, GenericRecord o2 ) {
            final long firstLocation = (Long) o1.get(SchemaUtils.LOCATION_FIELD_NAME);
            final long secondLocation = (Long) o2.get(SchemaUtils.LOCATION_FIELD_NAME);
            return Long.compare(firstLocation, secondLocation);
        }
    };

}
