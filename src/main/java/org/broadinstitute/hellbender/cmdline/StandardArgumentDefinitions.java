package org.broadinstitute.hellbender.cmdline;

/**
 * A set of String constants in which the name of the constant (minus the _SHORT_NAME suffix)
 * is the standard long Option name, and the value of the constant is the standard shortName.
 */
public final class StandardArgumentDefinitions {

    public static final String INPUT_LONG_NAME = "input";
    public static final String OUTPUT_LONG_NAME = "output";
    public static final String REFERENCE_LONG_NAME = "reference";
    public static final String VARIANT_LONG_NAME = "variant";
    public static final String FEATURE_LONG_NAME = "feature";
    public static final String USE_ORIGINAL_QUALITIES_LONG_NAME = "useOriginalQualities";
    public static final String LENIENT_LONG_NAME = "lenient";

    public static final String INPUT_SHORT_NAME = "I";
    public static final String OUTPUT_SHORT_NAME = "O";
    public static final String REFERENCE_SHORT_NAME = "R";
    public static final String VARIANT_SHORT_NAME = "V";
    public static final String FEATURE_SHORT_NAME = "F";
    public static final String LENIENT_SHORT_NAME = "LE";
    public static final String SAMPLE_ALIAS_SHORT_NAME = "ALIAS";
    public static final String LIBRARY_NAME_SHORT_NAME = "LIB";
    public static final String EXPECTED_INSERT_SIZE_SHORT_NAME = "INSERT";
    public static final String SEQUENCE_DICTIONARY_SHORT_NAME = "SD";
    public static final String METRICS_FILE_SHORT_NAME = "M";
    public static final String ASSUME_SORTED_SHORT_NAME = "AS";
    public static final String PF_READS_ONLY_SHORT_NAME = "PF";
    public static final String MINIMUM_MAPPING_QUALITY_SHORT_NAME = "MQ";
    public static final String READ_GROUP_ID_SHORT_NAME = "RG";
    public static final String PROGRAM_RECORD_ID_SHORT_NAME = "PG";
    public static final String MINIMUM_LOD_SHORT_NAME = "LOD";
    public static final String SORT_ORDER_SHORT_NAME = "SO";
    public static final String USE_ORIGINAL_QUALITIES_SHORT_NAME = "OQ";
}
