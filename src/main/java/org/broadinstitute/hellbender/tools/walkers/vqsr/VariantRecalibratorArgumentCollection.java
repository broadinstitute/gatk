package org.broadinstitute.hellbender.tools.walkers.vqsr;


import org.broadinstitute.hellbender.cmdline.Argument;

/*
 * A bunch of annotations for VQSR.
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
     * This parameter determines the maximum number of Gaussians that should be used when building a positive model
     * using the variational Bayes algorithm.
     */
    @Argument(fullName = "maxGaussians", shortName = "mG", doc = "Max number of Gaussians for the positive model", optional = true)
    public int MAX_GAUSSIANS = 8;

    /**
     * This parameter determines the maximum number of Gaussians that should be used when building a negative model
     * using the variational Bayes algorithm. The actual maximum used is the smaller value between the mG and mNG
     * arguments, meaning that if -mG is smaller than -mNG, -mG will be used for both. Note that this number should
     * be small (e.g. 4) to achieve the best results.
     */
    @Argument(fullName = "maxNegativeGaussians", shortName = "mNG", doc = "Max number of Gaussians for the negative model", optional = true)
    public int MAX_GAUSSIANS_FOR_NEGATIVE_MODEL = 2;

    /**
     * This parameter determines the maximum number of VBEM iterations to be performed in the variational Bayes algorithm.
     * The procedure will normally end when convergence is detected.
     */
    @Argument(fullName = "maxIterations", shortName = "mI", doc = "Maximum number of VBEM iterations", optional = true)
    public int MAX_ITERATIONS = 150;

    /**
     * This parameter determines the number of k-means iterations to perform in order to initialize the means of
     * the Gaussians in the Gaussian mixture model.
     */
    @Argument(fullName = "numKMeans", shortName = "nKM", doc = "Number of k-means iterations", optional = true)
    public int NUM_KMEANS_ITERATIONS = 100;

    /**
     * If a variant has annotations more than -std standard deviations away from mean, it won't be used for building
     * the Gaussian mixture model.
     */
    @Argument(fullName = "stdThreshold", shortName = "std", doc = "Annotation value divergence threshold (number of standard deviations from the means) ", optional = true)
    public double STD_THRESHOLD = 10.0;

    @Argument(fullName = "shrinkage", shortName = "shrinkage", doc = "The shrinkage parameter in the variational Bayes algorithm.", optional = true)
    public double SHRINKAGE = 1.0;

    @Argument(fullName = "dirichlet", shortName = "dirichlet", doc = "The dirichlet parameter in the variational Bayes algorithm.", optional = true)
    public double DIRICHLET_PARAMETER = 0.001;

    @Argument(fullName = "priorCounts", shortName = "priorCounts", doc = "The number of prior counts to use in the variational Bayes algorithm.", optional = true)
    public double PRIOR_COUNTS = 20.0;

    /**
     * The number of variants to use in building the Gaussian mixture model. Training sets larger than this will be randomly downsampled.
     */
    @Argument(fullName = "maxNumTrainingData", shortName = "maxNumTrainingData", doc = "Maximum number of training data", optional = true)
    protected int MAX_NUM_TRAINING_DATA = 2500000;

    /**
     * This parameter determines the minimum number of variants that will be selected from the list of worst scoring
     * variants to use for building the Gaussian mixture model of bad variants.
     */
    @Argument(fullName = "minNumBadVariants", shortName = "minNumBad", doc = "Minimum number of bad variants", optional = true)
    public int MIN_NUM_BAD_VARIANTS = 1000;

    /**
     * Variants scoring lower than this threshold will be used to build the Gaussian model of bad variants.
     */
    @Argument(fullName = "badLodCutoff", shortName = "badLodCutoff", doc = "LOD score cutoff for selecting bad variants", optional = true)
    public double BAD_LOD_CUTOFF = -5.0;

}