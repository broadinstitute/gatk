package org.broadinstitute.hellbender.cmdline;

/**
 * Created by davidben on 11/30/15.
 */
public final class ExomeStandardArgumentDefinitions {
    public static final String SEGMENT_FILE_LONG_NAME = "segments";
    public static final String SEGMENT_FILE_SHORT_NAME = "S";

    //file of target intervals
    public static final String TARGET_FILE_LONG_NAME = "targets";
    public static final String TARGET_FILE_SHORT_NAME = "T";

    public static final String LEGACY_SEG_FILE_LONG_NAME = "legacy";
    public static final String LEGACY_SEG_FILE_SHORT_NAME = "LEG";

    public static final String NORMAL_BAM_FILE_LONG_NAME = "normal";
    public static final String NORMAL_BAM_FILE_SHORT_NAME = "N";

    public static final String TUMOR_BAM_FILE_LONG_NAME = "tumor";
    public static final String TUMOR_BAM_FILE_SHORT_NAME = "T";

    //IntervalList of common SNP sites
    public static final String SNP_FILE_LONG_NAME = "snpIntervals";
    public static final String SNP_FILE_SHORT_NAME = "SNP";

    public static final String SITES_FILE_LONG_NAME = "siteIntervals";
    public static final String SITES_FILE_SHORT_NAME = "SITES";

    public static final String ALLELIC_COUNTS_FILE_LONG_NAME = "hets";
    public static final String ALLELIC_COUNTS_FILE_SHORT_NAME = "HET";

    public static final String NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME = "normalHets";
    public static final String NORMAL_ALLELIC_COUNTS_FILE_SHORT_NAME = "NHET";

    public static final String TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME = "tumorHets";
    public static final String TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME = "THET";

    public static final String PON_FILE_LONG_NAME = "panelOfNormals";
    public static final String PON_FILE_SHORT_NAME = "PON";

    public static final String ALLELIC_PON_FILE_LONG_NAME = "allelicPanelOfNormals";
    public static final String ALLELIC_PON_FILE_SHORT_NAME = "aPON";

    public static final String AF_PARAMETER_FILE_LONG_NAME = "paramAF";
    public static final String AF_PARAMETER_FILE_SHORT_NAME = "AF";

    //are inputs in log_2 space
    public static final String LOG2_LONG_NAME= "log2Input";
    public static final String LOG2_SHORT_NAME = "LOG";

    public static final String PRE_TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME = "preTangentNormalized";
    public static final String PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME = "PTN";

    public static final String TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME = "tangentNormalized";
    public static final String TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME = "TN";
}
