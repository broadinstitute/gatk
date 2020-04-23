package org.broadinstitute.hellbender.tools.variantdb;

import com.google.common.collect.ImmutableSet;

import java.util.Arrays;
import java.util.List;

public class SchemaConstants {

    public static final String POSITION_FIELD_NAME = "position";
    public static final String VALUES_ARRAY_FIELD_NAME = "values";

    public static final String SAMPLE_FIELD_NAME = "sample";
    public static final String SAMPLE_ID_FIELD_NAME = "sample_id";
    public static final String STATE_FIELD_NAME = "state";
    public static final String REF_ALLELE_FIELD_NAME = "ref";
    public static final String ALT_ALLELE_FIELD_NAME = "alt";

    public static final String GENOTYPE_FIELD_PREFIX = "call_";
    public static final String MULTIVALUE_FIELD_DELIMITER = ",";

    public static final ImmutableSet<String> REQUIRED_FIELDS = ImmutableSet.of(
            SAMPLE_FIELD_NAME,
            STATE_FIELD_NAME,
            REF_ALLELE_FIELD_NAME,
            ALT_ALLELE_FIELD_NAME
    );

    public static final List<String> COHORT_FIELDS = Arrays.asList(POSITION_FIELD_NAME, SAMPLE_FIELD_NAME, STATE_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME, "call_GT", "call_GQ", "call_RGQ");
    public static final List<String> SAMPLE_FIELDS = Arrays.asList(SAMPLE_FIELD_NAME);
    public static final List<String> YNG_FIELDS = Arrays.asList(POSITION_FIELD_NAME, REF_ALLELE_FIELD_NAME, ALT_ALLELE_FIELD_NAME);

//    private void validateSchema(final Set<String> columnNames) {
//        for ( final String requiredField : REQUIRED_FIELDS ) {
//            if ( ! columnNames.contains(requiredField) ) {
//                throw new UserException("Missing required column: " + requiredField +
//                        ". Actual columns encountered were: " + columnNames);
//            }
//        }
//    }


}
