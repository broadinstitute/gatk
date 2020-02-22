package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;

public class DragstrCasesSamplerArgumentCollection {
    @Argument(fullName = DragstrConstants.RANDOM_SEED_ARGUMENT_FULL_NAME,
            doc = "random number generator base seed",
            optional = true)
    public int randomSeed;
    @Argument(fullName = DragstrConstants.SAMPLING_PADDING_ARGUMENT_FULL_NAME,
            doc = "bases on either side of the repeat that are included in the STR pileup",
            optional = true)
    public int pileupPadding;
    @Argument(fullName = DragstrConstants.SAMPLING_MIN_MQ_ARGUMENT_FULL_NAME,
            doc = "the minimum read mapping quality allowed in sampled loci. Any read with a lower MQ will result in discarding that locus",
            optional = true)
    public int samplingMinMQ;
    @Argument(fullName = DragstrConstants.SAMPLING_MAX_MQ_ARGUMENT_FULL_NAME,
            doc = "the maximum number of sites to sample per period and repeat count combination",
            optional = true)
    public int maxCount;
    @Argument(fullName = DragstrConstants.SAMPLING_MIN_BQ_THRESHOLD_ARGUMENT_FULL_NAME,
            doc = "Base quality thershold for base-call. Reads that contain a number bases that map on the STR will not qualify for sampling",
            optional = true)
    public int baseQualThreshold;

    @Argument(fullName = org.broadinstitute.hellbender.utils.pairhmm.DragstrConstants.SAMPLING_MAX_BQ_EXCEPTIONS_ALLOWED_ARGUMENT_FULL_NAME,
        doc ="Maximum number of STR overlapping base call with low quality allowed for any read to be considered for further analysis",
        optional =true)
    public int baseQualExceptionsAllowed;

    @Argument(fullName = org.broadinstitute.hellbender.utils.pairhmm.DragstrConstants.SAMPLING_MAX_CASE_COUNT_ARGUMENT_FULL_NAME,
        doc ="maximum number of sites sampled for each combination of period and repeat count",
        optional =true)
    public int maximumNumberOfCases;

    @Argument(fullName = org.broadinstitute.hellbender.utils.pairhmm.DragstrConstants.SAMPLING_MIN_NON_REF_CONTAING_CASE_COUNT_ARGUMENT_FULL_NAME,
        doc ="targeted minimum number of sites sampled that contain non-ref reads for each combination of period and repeat count",
        optional =true)
    public int targetMinimumNonRefCases;

    @Argument(fullName = org.broadinstitute.hellbender.utils.pairhmm.DragstrConstants.SAMPLING_MIN_CASE_COUNT_ARGUMENT_FULL_NAME,
        doc ="minimum number of sites sampled for each combination of period and repeat count",
        optional =true)
    public int minimumNumberOfCases;

    public DragstrCasesSamplerArgumentCollection() {
        this.randomSeed = DragstrConstants.DEFAULT_RANDOM_SEED;
        this.pileupPadding = DragstrConstants.DEFAULT_SAMPLING_PADDING;
        this.samplingMinMQ = DragstrConstants.DEFAULT_SAMPLING_MIN_MQ;
        this.maxCount = DragstrConstants.DEFAULT_SAMPLING_MAX_COUNT;
        this.baseQualThreshold = DragstrConstants.DEFAULT_SAMPLING_BASE_QUAL_THRESHOLD;
        this.baseQualExceptionsAllowed = DragstrConstants.DEFAULT_SAMPLING_MAX_BQ_EXCEPTIONS_ALLOWED;
        this.maximumNumberOfCases = DragstrConstants.DEFAULT_SAMPLING_MAX_COUNT;
        this.targetMinimumNonRefCases = DragstrConstants.DEFAULT_SAMPLING_MIN_NON_REF_CONTAINING_CASE_COUNT;
        this.minimumNumberOfCases = DragstrConstants.DEFAULT_SAMPLING_MIN_CASE_COUNT;
    }
}