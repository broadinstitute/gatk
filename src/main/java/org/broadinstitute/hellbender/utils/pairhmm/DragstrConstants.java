package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.barclay.argparser.Argument;

public final class DragstrConstants {

    public static final String DRAGSTRINFO_KEY = "DRAGstrInfo";
    public static final String DRAGSTRPARAMS_KEY = "DRAGstrParams";

    public static final String USE_ALL_EVIDENCE_FULL_NAME = "use-all-evidence";
    public static final String MAX_PERIOD_ARGUMENT_FULL_NAME = "max-period";
    public static final String MAX_REPEATS_ARGUMENT_FULL_NAME = "max-repeats";
    public static final String RANDOM_SEED_ARGUMENT_FULL_NAME = "random-seed";
    public static final String SAMPLING_LOCI_ARGUMENT_FULL_NAME = "sampling-loci-zip";
    public static final String SAMPLING_PADDING_ARGUMENT_FULL_NAME = "pileup-padding";
    public static final String SAMPLING_MIN_MQ_ARGUMENT_FULL_NAME = "sampling-min-mq";
    public static final String SAMPLING_MAX_MQ_ARGUMENT_FULL_NAME = "sampling-max-count";
    public static final String SAMPLING_MIN_BQ_THRESHOLD_ARGUMENT_FULL_NAME = "sampling-high-bq-threshold";
    public static final String SAMPLING_MAX_BQ_EXCEPTIONS_ALLOWED_ARGUMENT_FULL_NAME = "sampling-max-low-bq-exceptions-per-read";
    public static final String SAMPLING_MAX_CASE_COUNT_ARGUMENT_FULL_NAME = "sampling-max-case-count";
    public static final String SAMPLING_MIN_NON_REF_CONTAING_CASE_COUNT_ARGUMENT_FULL_NAME = "sampling-min-variant-case-count";
    public static final String SAMPLING_MIN_CASE_COUNT_ARGUMENT_FULL_NAME = "sampling-min-case-count";


    public static final int DEFAULT_MAX_PERIOD = 8;
    public static final int DEFAULT_MAX_REPEATS = 20;
    public static final int DEFAULT_RANDOM_SEED = 23;
    public static final int DEFAULT_SAMPLING_PADDING = 5;
    public static final int DEFAULT_SAMPLING_MIN_MQ = 20;
    public static final int DEFAULT_SAMPLING_MAX_COUNT = 2000;
    public static final int DEFAULT_SAMPLING_BASE_QUAL_THRESHOLD = 10;
    public static final int DEFAULT_SAMPLING_MAX_BQ_EXCEPTIONS_ALLOWED = 2;
    public static final int DEFAULT_SAMPLING_MIN_NON_REF_CONTAINING_CASE_COUNT = 200;
    public static final int DEFAULT_SAMPLING_MIN_CASE_COUNT = 400;
    public static final String MIN_DEPTH_ARGUMENT_FULL_NAME = "minimum-depth";
    public static final String MIN_DEPTH_ARGUMENT_SHORT_NAME = "md";

    public static final String GP_VALUES_ARGUMENT_FULL_NAME = "gp-values";
    public static final String API_VALUES_ARGUMENT_FULL_NAME = "api-values";
    public static final String HET_TO_HOM_RATIO_FULL_NAME = "het-to-hom-ratio";
    public static final String MIN_LOCI_COUNT_FULL_NAME = "min-loci-count";
    public static final String API_MONO_THRESHOLD_FULL_NAME = "api-mono-threshold";
    public static final String MIN_GOP_FULL_NAME = "min-gop";
    public static final String MAX_GOP_FULL_NAME = "max-gop";
    public static final DoubleSequence DEFAULT_PHRED_GP_VALUES = new DoubleSequence("10:1.0:50");
    public static final DoubleSequence DEFAULT_PHRED_API_VALUES = new DoubleSequence("0:1.0:40");
    public static final double DEFAULT_HET_TO_HOM_RATIO = 2.0;
    public static final int DEFAULT_MIN_LOCI_COUNT = 50;
    public static final int DEFAULT_API_MONO_THRESHOLD = 3;
    public static final double DEFAULT_MIN_GOP = 10;
    public static final double DEFAULT_MAX_GOP = 50;
}
