package org.broadinstitute.hellbender.tools.variantdb;

import com.google.common.collect.ImmutableSet;

import java.util.Arrays;
import java.util.List;

public class SchemaUtils {

//    public static final String POSITION_FIELD_NAME = "position";
//    public static final String CHROM_FIELD_NAME = "chrom";
    public static final String LOCATION_FIELD_NAME = "location";

    public static final String SAMPLE_NAME_FIELD_NAME = "sample_name";
    public static final String SAMPLE_ID_FIELD_NAME = "sample_id";
    public static final String STATE_FIELD_NAME = "state";
    public static final String REF_ALLELE_FIELD_NAME = "ref";
    public static final String ALT_ALLELE_FIELD_NAME = "alt";

    public static final String GENOTYPE_FIELD_PREFIX = "call_";
    public static final String MULTIVALUE_FIELD_DELIMITER = ",";

    public static final ImmutableSet<String> REQUIRED_FIELDS = ImmutableSet.of(
            SAMPLE_NAME_FIELD_NAME,
            STATE_FIELD_NAME,
            REF_ALLELE_FIELD_NAME,
            ALT_ALLELE_FIELD_NAME
    );

    public static final List<String> COHORT_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, SAMPLE_NAME_FIELD_NAME, STATE_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME, "call_GT", "call_GQ", "call_RGQ");
    public static final List<String> ARRAY_COHORT_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, SAMPLE_NAME_FIELD_NAME, STATE_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME, "call_GT", "call_GQ");
    public static final List<String> SAMPLE_FIELDS = Arrays.asList(SAMPLE_NAME_FIELD_NAME);
    public static final List<String> YNG_FIELDS = Arrays.asList(LOCATION_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME);

//    private void validateSchema(final Set<String> columnNames) {
//        for ( final String requiredField : REQUIRED_FIELDS ) {
//            if ( ! columnNames.contains(requiredField) ) {
//                throw new UserException("Missing required column: " + requiredField +
//                        ". Actual columns encountered were: " + columnNames);
//            }
//        }
//    }

    public enum ChromosomeEnum {
        chr1(1),
        chr2(2),
        chr3(3),
        chr4(4),
        chr5(5),
        chr6(6),
        chr7(7),
        chr8(8),
        chr9(9),
        chr10(10),
        chr11(11),
        chr12(12),
        chr13(13),
        chr14(14),
        chr15(15),
        chr16(16),
        chr17(17),
        chr18(18),
        chr19(19),
        chr20(20),
        chr21(21),
        chr22(22),
        chrX(23),
        chrY(24),
        chrM(25);

        int index;

        ChromosomeEnum(int index) {
            this.index = index;
        }

        public static ChromosomeEnum valueOfIndex(int index) {
            switch(index) {
                case 1: return chr1;
                case 2: return chr2;
                case 3: return chr3;
                case 4: return chr4;
                case 5: return chr5;
                case 6: return chr6;
                case 7: return chr7;
                case 8: return chr8;
                case 9: return chr9;
                case 10: return chr10;
                case 11: return chr11;
                case 12: return chr12;
                case 13: return chr13;
                case 14: return chr14;
                case 15: return chr15;
                case 16: return chr16;
                case 17: return chr17;
                case 18: return chr18;
                case 19: return chr19;
                case 20: return chr20;
                case 21: return chr21;
                case 22: return chr22;
                case 23: return chrX;
                case 24: return chrY;
                case 25: return chrM;
                default: return null;
            }
        }
    }

    public static final long chromAdjustment = 1000000000000L;

    public static long encodeLocation(String chrom, int position) {
        int chromosomeIndex = ChromosomeEnum.valueOf(chrom.toUpperCase()).index;
        long adjustedLocation = Long.valueOf(chromosomeIndex) * chromAdjustment + Long.valueOf(position);
        return adjustedLocation;
    }

    public static ChromosomeEnum decodeContig(long location) {
        return ChromosomeEnum.valueOfIndex((int)(location/chromAdjustment));
    }

    public static int decodePosition(long location) {
        return (int)(location % chromAdjustment);
    }
}
