package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;

/**
 * Represents the panel of normals used for allele-bias correction.  See docs/CNVs/CNV-methods.pdf.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPanelOfNormals {
    private static final Logger logger = LogManager.getLogger(AllelicPanelOfNormals.class);

    public static final AllelicPanelOfNormals EMPTY_PON = new AllelicPanelOfNormals();

    private final Map<SimpleInterval, HyperparameterValues> siteToHyperparameterPairMap = new HashMap<>();
    private final HyperparameterValues mleHyperparameterValues;
    private final double mleMeanBias;
    private final double mleBiasVariance;

    private AllelicPanelOfNormals() {
        mleHyperparameterValues = new HyperparameterValues(Double.NaN, Double.NaN);
        mleMeanBias = Double.NaN;
        mleBiasVariance = Double.NaN;
    }

    /**
     * Constructs an allelic panel of normals from an {@link AllelicCountCollection} that contains
     * total alt and ref counts observed across all normals at each site.
     * @param counts    total alt and ref counts observed across all normals at each site
     */
    public AllelicPanelOfNormals(final AllelicCountCollection counts) {
        mleHyperparameterValues = calculateMLEHyperparameterValues(counts);
        mleMeanBias = meanBias(mleHyperparameterValues.alpha, mleHyperparameterValues.beta);
        mleBiasVariance = biasVariance(mleHyperparameterValues.alpha, mleHyperparameterValues.beta);
        initializeSiteToHyperparameterPairMap(counts);
    }

    /**
     * Constructs an allelic panel of normals from a file that contains
     * total alt and ref counts observed across all normals at each site.
     * @param inputFile    contains total alt and ref counts observed across all normals at each site
     */
    public AllelicPanelOfNormals(final File inputFile) {
        this(new AllelicCountCollection(validateFile(inputFile)));
    }

    private static File validateFile(final File inputFile) {
        Utils.nonNull(inputFile);
        Utils.regularReadableUserFile(inputFile);
        return inputFile;
    }

    /**
     * Gets the reference-bias alpha hyperparameter at a given SNP site if it is in the panel of normals
     * and the MLE alpha hyperparameter across all sites if it is not.
     * @param site  SNP site
     * @return      reference-bias alpha hyperparameter if site is in panel of normals,
     *              MLE alpha hyperparameter across all sites if it is not
     */
    public double getAlpha(final SimpleInterval site) {
        throwExceptionIfPoNIsEmpty();
        Utils.nonNull(site);
        return siteToHyperparameterPairMap.getOrDefault(site, mleHyperparameterValues).alpha;
    }

    /**
     * Gets the reference-bias beta hyperparameter at a given SNP site if it is in the panel of normals
     * and the MLE beta hyperparameter across all sites if it is not.
     * @param site  SNP site
     * @return      reference-bias beta hyperparameter if site is in panel of normals,
     *              MLE beta hyperparameter across all sites if it is not
     */
    public double getBeta(final SimpleInterval site) {
        throwExceptionIfPoNIsEmpty();
        Utils.nonNull(site);
        return siteToHyperparameterPairMap.getOrDefault(site, mleHyperparameterValues).beta;
    }

    /**
     * Gets the MLE mean-bias hyperparameter across all sites.
     * @return  MLE mean-bias hyperparameter across all sites
     */
    public double getMLEMeanBias() {
        throwExceptionIfPoNIsEmpty();
        return mleMeanBias;
    }

    /**
     * Gets the MLE bias-variance hyperparameter across all sites.
     * @return  MLE bias-variance hyperparameter across all sites
     */
    public double getMLEBiasVariance() {
        throwExceptionIfPoNIsEmpty();
        return mleBiasVariance;
    }

    private class HyperparameterValues {
        private final double alpha;
        private final double beta;

        private HyperparameterValues(final double alpha, final double beta) {
            this.alpha = alpha;
            this.beta = beta;
        }

        /**
         * Initializes the hyperparameter values at a site given the observed counts in all normals.
         * @param a     total alt counts observed across all normals at site
         * @param r     total ref counts observed across all normals at site
         */
        private HyperparameterValues(final int a, final int r) {
            final double f = 0.5;
            final int n = a + r;
            final double alpha = mleHyperparameterValues.alpha;
            final double beta = mleHyperparameterValues.beta;
            final double lambda0 = AlleleFractionLikelihoods.biasPosteriorMode(alpha, beta, f, a, r);
            final double kappa = AlleleFractionLikelihoods.biasPosteriorCurvature(alpha, f, r, n, lambda0);
            this.alpha = AlleleFractionLikelihoods.biasPosteriorEffectiveAlpha(lambda0, kappa);
            this.beta = AlleleFractionLikelihoods.biasPosteriorEffectiveBeta(lambda0, kappa);
        }
    }

    /**
     * Find MLE hyperparameter values from counts in panel of normals.  See analogous code in {@link AlleleFractionInitializer}.
     */
    private HyperparameterValues calculateMLEHyperparameterValues(final AllelicCountCollection counts) {
        double meanBias = AlleleFractionInitializer.INITIAL_MEAN_BIAS;
        double biasVariance = AlleleFractionInitializer.INITIAL_BIAS_VARIANCE;
        double previousIterationLogLikelihood;
        double nextIterationLogLikelihood = Double.NEGATIVE_INFINITY;
        logger.info(String.format("Initializing MLE hyperparameter values for allelic panel of normals.  Iterating until log likelihood converges to within %.3f.",
                AlleleFractionInitializer.LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD));
        int iteration = 1;
        do {
            previousIterationLogLikelihood = nextIterationLogLikelihood;
            meanBias = estimateMeanBias(meanBias, biasVariance, counts);
            biasVariance = estimateBiasVariance(meanBias, biasVariance, counts);
            nextIterationLogLikelihood = AlleleFractionLikelihoods.logLikelihoodForAllelicPanelOfNormals(meanBias, biasVariance, counts);
            logger.info(String.format("Iteration %d, model log likelihood = %.3f.", iteration, nextIterationLogLikelihood));
            iteration++;
        } while (iteration < AlleleFractionInitializer.MAX_ITERATIONS &&
                nextIterationLogLikelihood - previousIterationLogLikelihood > AlleleFractionInitializer.LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD);

        final double alpha = alpha(meanBias, biasVariance);
        final double beta = beta(meanBias, biasVariance);
        logger.info("MLE hyperparameter values for allelic panel of normals found:");
        logger.info("alpha = " + alpha);
        logger.info("beta = " + beta);
        return new HyperparameterValues(alpha, beta);
    }

    private static double estimateMeanBias(final double meanBias, final double biasVariance, final AllelicCountCollection counts) {
        final Function<Double, Double> objective = proposedMeanBias ->
                AlleleFractionLikelihoods.logLikelihoodForAllelicPanelOfNormals(proposedMeanBias, biasVariance, counts);
        return OptimizationUtils.argmax(objective, 0.0, AlleleFractionInitializer.MAX_REASONABLE_MEAN_BIAS, meanBias);
    }

    private static double estimateBiasVariance(final double meanBias, final double biasVariance, final AllelicCountCollection counts) {
        final Function<Double, Double> objective = proposedBiasVariance ->
                AlleleFractionLikelihoods.logLikelihoodForAllelicPanelOfNormals(meanBias, proposedBiasVariance, counts);
        return OptimizationUtils.argmax(objective, 0.0, AlleleFractionInitializer.MAX_REASONABLE_BIAS_VARIANCE, biasVariance);
    }

    private void initializeSiteToHyperparameterPairMap(final AllelicCountCollection counts) {
        logger.info("Initializing allelic panel of normals...");
        for (final AllelicCount count : counts.getCounts()) {
            final SimpleInterval site = count.getInterval();
            final HyperparameterValues hyperparameterValues = new HyperparameterValues(count.getAltReadCount(), count.getRefReadCount());
            if (siteToHyperparameterPairMap.containsKey(site)) {
                throw new UserException.BadInput("Input file for allelic panel of normals contains duplicate sites.");
            } else {
                siteToHyperparameterPairMap.put(site, hyperparameterValues);
            }
        }
        logger.info("Allelic panel of normals initialized.");
    }

    private void throwExceptionIfPoNIsEmpty() {
        if (this.equals(EMPTY_PON)) {
            throw new UnsupportedOperationException("Cannot get MLE hyperparameters for empty panel of normals.");
        }
    }

    private double alpha(final double meanBias, final double biasVariance) {
        return meanBias * meanBias / biasVariance;
    }

    private static double beta(final double meanBias, final double biasVariance) {
        return meanBias / biasVariance;
    }

    private static double meanBias(final double alpha, final double beta) {
        return alpha / beta;
    }

    private static double biasVariance(final double alpha, final double beta) {
        return alpha / (beta * beta);
    }
}
