package org.broadinstitute.hellbender.tools.variantdb;

import com.google.common.collect.ImmutableSet;

import java.util.Arrays;
import java.util.List;

public class SchemaUtils {

    public static final String POSITION_FIELD_NAME = "position";
    public static final String CHROM_FIELD_NAME = "chrom";
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

    public static final List<String> COHORT_FIELDS = Arrays.asList(CHROM_FIELD_NAME, POSITION_FIELD_NAME, SAMPLE_NAME_FIELD_NAME, STATE_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME, "call_GT", "call_GQ", "call_RGQ");
    public static final List<String> SAMPLE_FIELDS = Arrays.asList(SAMPLE_NAME_FIELD_NAME);
    public static final List<String> YNG_FIELDS = Arrays.asList(CHROM_FIELD_NAME, POSITION_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME);

//    private void validateSchema(final Set<String> columnNames) {
//        for ( final String requiredField : REQUIRED_FIELDS ) {
//            if ( ! columnNames.contains(requiredField) ) {
//                throw new UserException("Missing required column: " + requiredField +
//                        ". Actual columns encountered were: " + columnNames);
//            }
//        }
//    }

    public enum ChromosomeEnum {
        Chr1(1),
        Chr2(2),
        Chr3(3),
        Chr4(4),
        Chr5(5),
        Chr6(6),
        Chr7(7),
        Chr8(8),
        Chr9(9),
        Chr10(10),
        Chr11(11),
        Chr12(12),
        Chr13(13),
        Chr14(14),
        Chr15(15),
        Chr16(16),
        Chr17(17),
        Chr18(18),
        Chr19(19),
        Chr20(20),
        Chr21(21),
        Chr22(22),
        ChrX(23),
        ChrY(24),
        ChrM(25);

        int index;

        ChromosomeEnum(int index) {
            this.index = index;
        }

        public static ChromosomeEnum valueOfIndex(int index) {
            switch(index) {
                case 1: return Chr1;
                case 2: return Chr2;
                case 3: return Chr3;
                case 4: return Chr4;
                case 5: return Chr5;
                case 6: return Chr6;
                case 7: return Chr7;
                case 8: return Chr8;
                case 9: return Chr9;
                case 10: return Chr10;
                case 11: return Chr11;
                case 12: return Chr12;
                case 13: return Chr13;
                case 14: return Chr14;
                case 15: return Chr15;
                case 16: return Chr16;
                case 17: return Chr17;
                case 18: return Chr18;
                case 19: return Chr19;
                case 20: return Chr20;
                case 21: return Chr21;
                case 22: return Chr22;
                case 23: return ChrX;
                case 24: return ChrY;
                case 25: return ChrM;
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
