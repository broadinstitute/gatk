package org.broadinstitute.hellbender.tools.exome.allelefraction;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.special.Beta;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;

import java.util.EnumMap;
import java.util.List;
import java.util.stream.Collectors;

/**
 * The Allele Fraction Model (after marginalizing latent parameters as described in docs/AllelicCapSeg/ACS-methods.pdf)
 * contains the following parameters:
 *      1.  minor allele fractions for each segment
 *      2.  a global outlier probability
 *      3.  the mean allelic bias
 *      4.  the rate (mean / variance) of the allelic bias
 *
 * Note that 3 and 4 are hyperparameters specifying a gamma distribution prior on allelic bias -- the latent variables
 * for bias at each het site have been marginalized but the hyperparameter have not.
 *
 * The Allele Fraction Model samples the distribution of these parameters using Markov chain Monte Carlo and in principle
 * an initialization step is not necessary.  However, in practice this initialization find the mode of the posterior
 * distributions in only a few iterations, whereas sampling would require many more.  Thus we greatly reduce the
 * number of burn-in samples that we must discard.
 *
 * The initialization is similar to the EM algorithm, in that it iteratively optimizes the likelihood averaged over
 * posterior distributions of the latent variables (bias ratios and per-het indicators of alt minor / ref minor /outlier).
 * A likely source of confusion is that while these latent variables do not appear in the MCMC sampling (because we sample
 * a marginalized likelihood) they are necessary here to perform the E step.  Note also that some of the M step likelihood
 * maximizations are only approximate / heuristic -- this is fine because initialization does not need to be exact.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionInitializer {
    @VisibleForTesting
    protected static final double INITIAL_OUTLIER_PROBABILITY = 0.01;
    private static final double BIAS_MEAN_INITIAL = 1.0;
    private static final double BIAS_VARIANCE_INITIAL = 0.1;   // this is an overestimate, but starting small makes it slow for
                                                            // mean bias to escape a bad initial guess
    private static final int NUMBER_OF_EM_ITERATIONS = 10;

    private final AlleleFractionData data;
    private AlleleFractionState state;
    private final List<Responsibility> responsibilities;

    /**
     * Data structure containing the E step posterior distribution (a 3-D categorical distribution over
     * alt minor / ref minor / outlier) for the indicator variable of a single het site.
     */
    @SuppressWarnings("serial")
    @VisibleForTesting
    protected static final class Responsibility extends EnumMap<AlleleFractionIndicator, Double> {

        /**
         *  initialize naively based on assumption of no allelic bias and unknown minor fraction f
         *  (hence integrate over f to get the marginal posterior on the indicator)
         *
         * that is, likelihood of alt minor is int_{0 to 1/2} f^a (1-f)^r df
         *          likelihood of ref minor is int_{0 to 1/2} f^r (1-f)^a df
         * these are proportional to I(1/2, a + 1, r + 1) and I(1/2, r + 1, a + 1),
         * respectively, where I is the (incomplete) regularized Beta function
         * By definition these likelihoods sum to 1, ie they are already normalized.
         *
         * @param count an AllelicCount containing the number of alt and ref reads at this het site
         */
        public Responsibility(final AllelicCount count) {
            super(AlleleFractionIndicator.class);
            final int a = count.getAltReadCount();
            final int r = count.getRefReadCount();
            double altMinorLk;
            try {
                altMinorLk = Beta.regularizedBeta(0.5, a + 1, r + 1);
            } catch(final MaxCountExceededException e) {
                altMinorLk = a < r ? 1.0 : 0.0; //if the special function can't be computed, give an all-or-nothing responsibility
            }
            final double refMinorLk = 1 - altMinorLk;

            put(AlleleFractionIndicator.ALT_MINOR, altMinorLk * (1 - INITIAL_OUTLIER_PROBABILITY));
            put(AlleleFractionIndicator.REF_MINOR, refMinorLk * (1 - INITIAL_OUTLIER_PROBABILITY));
            put(AlleleFractionIndicator.OUTLIER, INITIAL_OUTLIER_PROBABILITY);
        }

        /**
         * given unnormalized log-space posteriors, set normalized real-space posteriors
         *
         * @param altMinorLogPosterior unnormalized log posterior of alt minor
         * @param refMinorLogPosterior unnormalized log posterior of ref minor
         * @param outlierLogPosterior  unnormalized log posterior of outlier
         */
        public void setFromLogPosteriors(final double altMinorLogPosterior, final double refMinorLogPosterior, final double outlierLogPosterior) {
            // in order to avoid overflow from exponentiation, subtract the maximum unnormalized log posterior
            // from each log posterior (subtracting in log space is division in real space)
            final double maxUnnormalizedLogPosterior = Doubles.max(altMinorLogPosterior, refMinorLogPosterior, outlierLogPosterior);

            final double altMinorUnnormalized = Math.exp(altMinorLogPosterior - maxUnnormalizedLogPosterior);
            final double refMinorUnnormalized = Math.exp(refMinorLogPosterior - maxUnnormalizedLogPosterior);
            final double outlierUnnormalized = Math.exp(outlierLogPosterior - maxUnnormalizedLogPosterior);

            final double normalization = altMinorUnnormalized + refMinorUnnormalized + outlierUnnormalized;

            put(AlleleFractionIndicator.ALT_MINOR, altMinorUnnormalized / normalization);
            put(AlleleFractionIndicator.REF_MINOR, refMinorUnnormalized / normalization);
            put(AlleleFractionIndicator.OUTLIER, outlierUnnormalized / normalization);
        }
    }

    /**
     * Do the initialization
     * @param data data
     */
    public AlleleFractionInitializer(final AlleleFractionData data) {
        this.data = data;
        responsibilities = data.getAllelicCounts().stream().map(c -> new Responsibility(c)).collect(Collectors.toList());
        state = new AlleleFractionState(BIAS_MEAN_INITIAL, BIAS_VARIANCE_INITIAL, INITIAL_OUTLIER_PROBABILITY, estimateMinorFractions());
        for (int iteration = 0; iteration < NUMBER_OF_EM_ITERATIONS; iteration++) {
            performSingleEMIteration();
        }
    }

    /**
     *
     * @return the initialized state of the Allele Fraction Model
     */
    public AlleleFractionState getInitializedState() { return state; }

    // END OF PUBLIC API

    private double pAltMinor(final int het) { return responsibilities.get(het).get(AlleleFractionIndicator.ALT_MINOR); }
    private double pRefMinor(final int het) { return responsibilities.get(het).get(AlleleFractionIndicator.REF_MINOR); }
    private double pOutlier(final int het) { return responsibilities.get(het).get(AlleleFractionIndicator.OUTLIER); }
    private double pNotOutlier(final int het) { return 1 - pOutlier(het); }

    //perform an E step to update responsibilities using high-level parameters from current AlleleFractionState
    private void updateResponsibilities() {
        for (int segment = 0; segment < data.numSegments(); segment++) {
            for (final int het : data.hetsInSegment(segment)) {
                final AllelicCount count= data.count(het);
                final double altMinorLogPosterior = AlleleFractionModeller.hetLogLikelihood(state, segment, count, AlleleFractionIndicator.ALT_MINOR);
                final double refMinorLogPosterior = AlleleFractionModeller.hetLogLikelihood(state, segment, count, AlleleFractionIndicator.REF_MINOR);
                final double outlierLogPosterior = AlleleFractionModeller.hetLogLikelihood(state, segment, count, AlleleFractionIndicator.OUTLIER);
                responsibilities.get(het).setFromLogPosteriors(altMinorLogPosterior, refMinorLogPosterior, outlierLogPosterior);
            }
        }
    }

    // M step for minor allele fractions -- set minor fraction to responsibility-weighted total count of
    // reads in minor allele divided by total reads, ignoring outliers
    private AlleleFractionState.MinorFractions estimateMinorFractions() {
        AlleleFractionState.MinorFractions result = new AlleleFractionState.MinorFractions();
        for (int segment = 0; segment < data.numSegments(); segment++) {
            final double minorReads = data.hetsInSegment(segment).stream()
                    .mapToDouble(het -> data.altCount(het) * pAltMinor(het) + data.refCount(het) * pRefMinor(het)).sum();
            final double totalReads = data.hetsInSegment(segment).stream().mapToDouble(het -> data.readCount(het) * pNotOutlier(het)).sum();

            // the "+ 1"s correspond to a flat prior on minor fractions
            result.add((minorReads + 1)/(totalReads + 1));
        }
        return result;
    }

    // M step for outlierProbability -- set it to the average of all OUTLIER responsibilities
    private double estimateOutlierProbability() {
        return responsibilities.stream().mapToDouble(r -> r.get(AlleleFractionIndicator.OUTLIER)).average().getAsDouble();
    }

    /**
     * Estimate the kth moment of the distribution of latent allelic bias parameters by averaging the posterior
     * moment over all het sites.  These are used in the M step for the bias hyperparameters, which
     * uses moment-matching as a heuristic substitute for likelihood maximization.
     *
     * @param order the order of the desired moment
     * @return the posterior expectation of bias^order, averaged over all het sites
     */
    private double estimateBiasMoment(final int order) {
        double numerator = 0.0;
        double denominator = 0.0;
        for (int segment = 0; segment < data.numSegments(); segment++) {
            for (int het : data.hetsInSegment(segment)) {
                final AllelicCount count= data.count(het);
                numerator += pAltMinor(het) * AlleleFractionModeller.biasPosteriorMoment(state, segment, count, AlleleFractionIndicator.ALT_MINOR, order)
                        + pRefMinor(het) * AlleleFractionModeller.biasPosteriorMoment(state, segment, count, AlleleFractionIndicator.REF_MINOR, order);
                denominator += pNotOutlier(het);
            }
        }
        return numerator / denominator;
    }

    private double estimateAverageBias() { return estimateBiasMoment(1); }
    private double estimateAverageSquaredBias() { return estimateBiasMoment(2); }

    private void performSingleEMIteration() {
        updateResponsibilities();

        //update minor allele fractions
        state = new AlleleFractionState(state.meanBias(), state.biasVariance(), state.outlierProbability(), estimateMinorFractions());

        //update minor allele fractions
        state = new AlleleFractionState(state.meanBias(), state.biasVariance(), estimateOutlierProbability(), state.minorFractions());

        //update gamma hyperparameters
        final double meanBias = estimateAverageBias();
        final double meanSquaredBias = estimateAverageSquaredBias();
        final double variance = meanSquaredBias - meanBias*meanBias;
        state = new AlleleFractionState(meanBias, variance, state.outlierProbability(), state.minorFractions());
    }
}
