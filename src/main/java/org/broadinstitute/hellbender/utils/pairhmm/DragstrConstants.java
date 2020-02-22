package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.barclay.argparser.Argument;

public final class DragstrConstants {

    public static final String DRAGSTRINFO_KEY = "DRAGstrInfo";
    public static final String DRAGSTRPARAMS_KEY = "DRAGstrParams";

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


    @Argument(fullName = SAMPLING_MIN_BQ_THRESHOLD_ARGUMENT_FULL_NAME,
            doc = "Base quality thershold for base-call. Reads that contain a number bases that map on the STR will not qualify for sampling",
            optional = true)
    private int baseQualThreshold = DEFAULT_SAMPLING_BASE_QUAL_THRESHOLD; // 10;

    @Argument(fullName = SAMPLING_MAX_BQ_EXCEPTIONS_ALLOWED_ARGUMENT_FULL_NAME,
            doc = "Maximum number of STR overlapping base call with low quality allowed for any read to be considered for further analysis",
            optional = true)
    private int baseQualExceptionsAllowed = DEFAULT_SAMPLING_MAX_BQ_EXCEPTIONS_ALLOWED; // 2;

    @Argument(fullName = SAMPLING_MAX_CASE_COUNT_ARGUMENT_FULL_NAME,
            doc = "maximum number of sites sampled for each combination of period and repeat count",
            optional = true)
    private int maximumNumberOfCases = DEFAULT_SAMPLING_MAX_COUNT; // 2000.

    @Argument(fullName = SAMPLING_MIN_NON_REF_CONTAING_CASE_COUNT_ARGUMENT_FULL_NAME,
            doc = "targeted minimum number of sites sampled that contain non-ref reads for each combination of period and repeat count",
            optional = true)
    private int targetMinimumNonRefCases = DEFAULT_SAMPLING_MIN_NON_REF_CONTAINING_CASE_COUNT; //200.

    @Argument(fullName = SAMPLING_MIN_CASE_COUNT_ARGUMENT_FULL_NAME,
            doc = "minimum number of sites sampled for each combination of period and repeat count",
            optional = true)
    private int minimumNumberOfCases = DEFAULT_SAMPLING_MIN_CASE_COUNT;
}
