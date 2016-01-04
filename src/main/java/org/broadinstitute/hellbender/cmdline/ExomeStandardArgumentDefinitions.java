package org.broadinstitute.hellbender.cmdline;

/**
 * Created by davidben on 11/30/15.
 */
public final class ExomeStandardArgumentDefinitions {
    public static final String SEGMENT_FILE_SHORT_NAME = "S";
    public static final String SEGMENT_FILE_LONG_NAME = "segments";

    //file of target intervals
    public static final String TARGET_FILE_SHORT_NAME = "T";
    public static final String TARGET_FILE_LONG_NAME = "targets";

    public static final String LEGACY_SEG_FILE_SHORT_NAME = "LEG";
    public static final String LEGACY_SEG_FILE_LONG_NAME = "legacy";

    public static final String SAMPLE_LONG_NAME = "sample";

    public static final String NORMAL_BAM_FILE_SHORT_NAME = "N";
    public static final String NORMAL_BAM_FILE_LONG_NAME = "normal";

    public static final String TUMOR_BAM_FILE_SHORT_NAME = "T";
    public static final String TUMOR_BAM_FILE_LONG_NAME = "tumor";

    public static final String SNP_FILE_SHORT_NAME = "SNP";
    public static final String SNP_FILE_LONG_NAME = "snpIntervals";

    public static final String NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME = "normalHets";
    public static final String NORMAL_ALLELIC_COUNTS_FILE_SHORT_NAME = "NHET";

    public static final String TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME = "tumorHets";
    public static final String TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME = "THET";

    public static final String PON_FILE_SHORT_NAME = "PON";
    public static final String PON_FILE_LONG_NAME = "panelOfNormals";

    //are inputs in log_2 space
    public static final String LOG2_LONG_NAME= "log2Input";
    public static final String LOG2_SHORT_NAME = "LOG";

    public static final String INCLUDE_SEX_CHROMOSOMES_LONG_NAME= "sexChromosomes";
    public static final String INCLUDE_SEX_CHROMOSOMES_SHORT_NAME = "schr";

    public static final String PRE_TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME = "preTangentNormalized";
    public static final String PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME = "PTN";

    public static final String TANGENT_NORMALIZED_COUNTS_LONG_NAME = "tangentNormalized";
    public static final String TANGENT_NORMALIZED_COUNTS_SHORT_NAME = "TN";
}
