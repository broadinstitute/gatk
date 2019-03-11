package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ThresholdCalculator {
    public enum Strategy {
        CONSTANT, FALSE_DISCOVERY_RATE, OPTIMAL_F_SCORE
    }

    private final Strategy strategy;
    private final double maxFalseDiscoveryRate;
    private final double fScoreBeta;

    private double threshold;

    final List<Double> artifactProbabilities = new ArrayList<>();

    public ThresholdCalculator(final Strategy strategy, final double initialThreshold, final double maxFalseDiscoveryRate, final double fScoreBeta) {
        this.strategy = strategy;
        this.threshold = initialThreshold;
        this.maxFalseDiscoveryRate = maxFalseDiscoveryRate;
        this.fScoreBeta = fScoreBeta;
    }

    public void addArtifactProbability(final double artifactProbability) {
        artifactProbabilities.add(artifactProbability);
    }

    public void relearnThresholdAndClearAcumulatedProbabilities() {
        switch (strategy) {
            case CONSTANT:  // don't adjust
                break;
            case FALSE_DISCOVERY_RATE:
                threshold = ThresholdCalculator.calculateThresholdBasedOnFalseDiscoveryRate(artifactProbabilities, maxFalseDiscoveryRate);
                break;
            case OPTIMAL_F_SCORE:
                threshold = ThresholdCalculator.calculateThresholdBasedOnOptimalFScore(artifactProbabilities, fScoreBeta);
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Invalid threshold strategy type: " + strategy + ".");
        }
        clear();
    }

    public void clear() {
        artifactProbabilities.clear();
    }

    public double getThreshold() {
        return threshold;
    }


    /**
     * Compute the filtering threshold that maximizes the F_beta score
     *
     * @param posteriors A list of posterior probabilities, which gets sorted
     * @param beta relative weight of recall to precision
     */
    @VisibleForTesting
    static double calculateThresholdBasedOnOptimalFScore(final List<Double> posteriors, final double beta){
        ParamUtils.isPositiveOrZero(beta, "requested F-score beta must be non-negative");

        Collections.sort(posteriors);

        final double expectedTruePositives = posteriors.stream()
                .mapToDouble(prob -> 1 - prob).sum();


        // starting from filtering everything (threshold = 0) increase the threshold to maximize the F score
        final MutableDouble truePositives = new MutableDouble(0);
        final MutableDouble falsePositives = new MutableDouble(0);
        final MutableDouble falseNegatives = new MutableDouble(expectedTruePositives);
        int optimalIndexInclusive = -1; // include all indices up to and including this. -1 mean filter all.
        double optimalFScore = 0;   // if you exclude everything, recall is zero

        final int N = posteriors.size();

        for (int n = 0; n < N; n++){
            truePositives.add(1 - posteriors.get(n));
            falsePositives.add(posteriors.get(n));
            falseNegatives.subtract(1 - posteriors.get(n));
            final double F = (1+beta*beta)*truePositives.getValue() /
                    ((1+beta*beta)*truePositives.getValue() + beta*beta*falseNegatives.getValue() + falsePositives.getValue());
            if (F >= optimalFScore) {
                optimalIndexInclusive = n;
                optimalFScore = F;
            }
        }

        return optimalIndexInclusive == -1 ? 0 : (optimalIndexInclusive == N - 1 ? 1 : posteriors.get(optimalIndexInclusive));
    }

    /**
     *
     * Compute the filtering threshold that ensures that the false positive rate among the resulting pass variants
     * will not exceed the requested false positive rate
     *
     * @param posteriors A list of posterior probabilities, which gets sorted
     * @param requestedFPR We set the filtering threshold such that the FPR doesn't exceed this value
     */
    @VisibleForTesting
    public static double calculateThresholdBasedOnFalseDiscoveryRate(final List<Double> posteriors, final double requestedFPR){
        ParamUtils.isPositiveOrZero(requestedFPR, "requested FPR must be non-negative");
        final double thresholdForFilteringNone = 1.0;
        final double thresholdForFilteringAll = 0.0;

        Collections.sort(posteriors);

        final int numPassingVariants = posteriors.size();
        double cumulativeExpectedFPs = 0.0;

        for (int i = 0; i < numPassingVariants; i++){
            final double posterior = posteriors.get(i);

            // One can show that the cumulative error rate is monotonically increasing in i
            final double expectedFPR = (cumulativeExpectedFPs + posterior) / (i + 1);
            if (expectedFPR > requestedFPR){
                return i > 0 ? posteriors.get(i-1) : thresholdForFilteringAll;
            }

            cumulativeExpectedFPs += posterior;
        }

        // If the expected FP rate never exceeded the max tolerable value, then we can let everything pass
        return thresholdForFilteringNone;
    }
}
