package org.broadinstitute.hellbender.tools.walkers.vqsr;


import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;

/*
 * A bunch of arguments for VQSR.
 */
final class VariantRecalibratorArgumentCollection {

    public enum Mode {
        SNP,
        INDEL,
        BOTH
    }

    /**
     * Use either SNP for recalibrating only SNPs (emitting indels untouched in the output VCF) or INDEL for indels (emitting SNPs untouched in the output VCF). There is also a BOTH option for recalibrating both SNPs and indels simultaneously, but this is meant for testing purposes only and should not be used in actual analyses.
     */
    @Argument(fullName = "mode", shortName = "mode", doc = "Recalibration mode to employ", optional = false)
    public VariantRecalibratorArgumentCollection.Mode MODE = VariantRecalibratorArgumentCollection.Mode.SNP;

    /**
     * Generate a VQSR model using per-allele data instead of the default per-site data, assuming that the input VCF contains allele-specific annotations.
     * Annotations should be specified using their full names with AS_ prefix. Non-allele-specific (scalar) annotations will be applied to all alleles.
     */
    @Argument(fullName="use-allele-specific-annotations",
            shortName="AS",
            doc="If specified, the variant recalibrator will attempt to use the allele-specific versions of the specified annotations.", optional=true)
    public boolean useASannotations = false;

    /**
     * This parameter determines the maximum number of Gaussians that should be used when building a positive model
     * using the variational Bayes algorithm.
     */
    @Advanced
    @Argument(fullName = "max-gaussians", doc = "Max number of Gaussians for the positive model", optional = true)
    public int MAX_GAUSSIANS = 8;

    /**
     * This parameter determines the maximum number of Gaussians that should be used when building a negative model
     * using the variational Bayes algorithm. The actual maximum used is the smaller value between the mG and mNG
     * arguments, meaning that if -mG is smaller than -mNG, -mG will be used for both. Note that this number should
     * be small (e.g. 4) to achieve the best results.
     */
    @Advanced
    @Argument(fullName = "max-negative-gaussians", doc = "Max number of Gaussians for the negative model", optional = true)
    public int MAX_GAUSSIANS_FOR_NEGATIVE_MODEL = 2;

    /**
     * This parameter determines the maximum number of VBEM iterations to be performed in the variational Bayes algorithm.
     * The procedure will normally end when convergence is detected.
     */
    @Advanced
    @Argument(fullName = "max-iterations", doc = "Maximum number of VBEM iterations", optional = true)
    public int MAX_ITERATIONS = 150;

    /**
     * This parameter determines the number of k-means iterations to perform in order to initialize the means of
     * the Gaussians in the Gaussian mixture model.
     */
    @Advanced
    @Argument(fullName = "k-means-iterations", doc = "Number of k-means iterations", optional = true)
    public int NUM_KMEANS_ITERATIONS = 100;

    /**
     * If a variant has annotations more than -std standard deviations away from mean, it won't be used for building
     * the Gaussian mixture model.
     */
    @Advanced
    @Argument(fullName = "standard-deviation-threshold", shortName = "std", doc = "Annotation value divergence threshold (number of standard deviations from the means) ", optional = true)
    public double STD_THRESHOLD = 10.0;

    @Advanced
    @Argument(fullName = "shrinkage", doc = "The shrinkage parameter in the variational Bayes algorithm.", optional = true)
    public double SHRINKAGE = 1.0;

    @Advanced
    @Argument(fullName = "dirichlet",doc = "The dirichlet parameter in the variational Bayes algorithm.", optional = true)
    public double DIRICHLET_PARAMETER = 0.001;

    @Advanced
    @Argument(fullName = "prior-counts", doc = "The number of prior counts to use in the variational Bayes algorithm.", optional = true)
    public double PRIOR_COUNTS = 20.0;

    /**
     * The number of variants to use in building the Gaussian mixture model. Training sets larger than this will be randomly downsampled.
     */
    @Advanced
    @Argument(fullName = "maximum-training-variants", doc = "Maximum number of training data", optional = true)
    protected int MAX_NUM_TRAINING_DATA = 2500000;

    /**
     * This parameter determines the minimum number of variants that will be selected from the list of worst scoring
     * variants to use for building the Gaussian mixture model of bad variants.
     */
    @Advanced
    @Argument(fullName = "minimum-bad-variants", doc = "Minimum number of bad variants", optional = true)
    public int MIN_NUM_BAD_VARIANTS = 1000;

    /**
     * Variants scoring lower than this threshold will be used to build the Gaussian model of bad variants.
     */
    @Advanced
    @Argument(fullName = "bad-lod-score-cutoff", shortName = "bad-lod-cutoff", doc = "LOD score cutoff for selecting bad variants", optional = true)
    public double BAD_LOD_CUTOFF = -5.0;

    /**
     * MQ is capped at a "max" value (60 for bwa-mem) when the alignment is considered perfect. Typically, a huge
     * proportion of the reads in a dataset are perfectly mapped, which yields a distribution of MQ values with a
     * blob below the max value and a huge peak at the max value. This does not conform to the expectations of the
     * Gaussian mixture model of VQSR and has been observed to yield a ROC curve with a jump.
     *
     * This argument aims to mitigate this problem. Using MQCap = X has 2 effects:  (1) MQs are transformed by a scaled
     * logit on [0,X] (+ epsilon to avoid division by zero) to make the blob more Gaussian-like and (2) the transformed
     * MQ=X are jittered to break the peak into a narrow Gaussian.
     *
     * Beware that IndelRealigner, if used, adds 10 to MQ for successfully realigned indels. We recommend to either use
     * --read-filter ReassignOriginalMQAfterIndelRealignment with HaplotypeCaller or use a MQCap=max+10 to take that
     * into account.
     *
     * If this option is not used, or if MQCap is set to 0, MQ will not be transformed.
     */
    @Advanced
    @Argument(fullName="mq-cap-for-logit-jitter-transform", shortName = "mq-cap", doc="Apply logit transform and jitter to MQ values", optional=true)
    public int MQ_CAP = 0;

    /**
     * The following 2 arguments are hidden because they are only for testing different jitter amounts with and without logit transform.
     * Once this will have been tested, and the correct jitter amount chosen (perhaps as a function of the logit range [0,max]) they can be removed.
     */

    @Hidden
    @Advanced
    @Argument(fullName = "no-mq-logit", doc="MQ is by default transformed to log[(MQ_cap + epsilon - MQ)/(MQ + epsilon)] to make it more Gaussian-like.  Use this flag to not do that.", optional = true)
    public boolean NO_MQ_LOGIT = false;

    @Hidden
    @Advanced
    @Argument(fullName="mq-jitter", doc="Amount of jitter (as a factor to a Normal(0,1) noise) to add to the MQ capped values", optional = true)
    public double MQ_JITTER = 0.05;

}