package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.util.Collection;
import java.util.Map;
import java.util.Objects;

public class EvidenceStatUtils {

    private static final Median MEDIAN = new Median();

    /**
     * Converts probability to Phred-scaled quality
     * @param pError  error probability
     * @param maxQual max quality score
     * @return        phred-scaled score, or null of pError is null
     */
    public static Double probToQual(final Double pError, final double maxQual) {
        return pError == null ? null : Math.min(-10.0*Math.log10(pError), maxQual);
    }

    /**
     * Gets carrier signal as percentage of total signal, i.e. 100 * carrierSignal / (carrierSignal + backgroundSignal).
     * Returns null if either input is null, otherwise 0 if either signal is 0.
     * @param carrierSignal     carrier evidence count
     * @param backgroundSignal  background evidence count
     * @return                  result on [0-100] or null
     */
    public static Double carrierSignalFraction(final Double carrierSignal, final Double backgroundSignal) {
        if (backgroundSignal == null || carrierSignal == null) {
            return null;
        }
        if (backgroundSignal == 0 && carrierSignal == 0) {
            return 0.;
        }
        return 100 * carrierSignal / (backgroundSignal + carrierSignal);
    }

    /**
     * Find the median of site counts in a subset of samples when normalized by sample coverage. Returns 0 if no samples.
     * @param samples               samples to restrict to
     * @param sampleCounts          sample counts
     * @param sampleCoverageMap     characteristic sample depth
     * @param targetCoverage        normalized characteristic coverage
     * @return                      result
     */
    @VisibleForTesting
    protected static double getMedianNormalizedCount(final Collection<String> samples,
                                                     final Map<String, Integer> sampleCounts,
                                                     final Map<String, Double> sampleCoverageMap,
                                                     final double targetCoverage) {
        if (samples.isEmpty() || sampleCounts.isEmpty()) {
            return 0;
        }
        final double[] normalizedCounts = new double[samples.size()];
        int i = 0;
        for (final String sample : samples) {
            normalizedCounts[i] = getNormalizedCount(sampleCounts.getOrDefault(sample, 0), sampleCoverageMap.get(sample), targetCoverage);
            i++;
        }
        return MEDIAN.evaluate(normalizedCounts);
    }

    /**
     * Normalizes counts by sample coverage to the target coverage.
     */
    public static double getNormalizedCount(final int count, final double sampleCoverage, final double targetCoverage) {
        return Math.round(targetCoverage * count / sampleCoverage);
    }

    /**
     * Poisson CDF
     */
    public static double cumulativePoissonProbability(final double mean, final int x) {
        if (x < 0) {
            return 0;
        }
        if (x == Integer.MAX_VALUE) {
            return 1;
        }
        return Gamma.regularizedGammaQ((double) x + 1, mean);
    }

    /**
     * Performs poisson test on a single site by computing the probability of observing the background counts
     * under a carrier count distribution
     *
     * @param sampleCounts      sample counts
     * @param carrierSamples    carrier sample IDs
     * @param backgroundSamples background sample IDs
     * @param meanCoverage      mean coverage of all samples
     * @return                  probability on [0,1]
     */
    public static PoissonTestResult calculateOneSamplePoissonTest(final Map<String, Integer> sampleCounts,
                                                                  final Collection<String> carrierSamples,
                                                                  final Collection<String> backgroundSamples,
                                                                  final Map<String, Double> sampleCoverageMap,
                                                                  final double meanCoverage) {
        final double medianNormalizedCarrierCount = EvidenceStatUtils.getMedianNormalizedCount(carrierSamples, sampleCounts, sampleCoverageMap, meanCoverage);
        final double medianBackgroundRate = EvidenceStatUtils.getMedianNormalizedCount(backgroundSamples, sampleCounts, sampleCoverageMap, meanCoverage);
        // If a common variant (AF > 0.5), clamp background to 0
        final int backgroundCount = carrierSamples.size() > backgroundSamples.size() ? 0 : (int) Math.round(medianBackgroundRate);
        final double p = EvidenceStatUtils.cumulativePoissonProbability(medianNormalizedCarrierCount, backgroundCount);
        return new PoissonTestResult(p, medianNormalizedCarrierCount, medianBackgroundRate);
    }

    public static final class PoissonTestResult {
        private final double p;
        private final double carrierSignal;
        private final double backgroundSignal;

        public PoissonTestResult(final double p, final double carrierSignal, final double backgroundSignal) {
            this.p = p;
            this.carrierSignal = carrierSignal;
            this.backgroundSignal = backgroundSignal;
        }

        /**
         * Probability
         */
        public double getP() {
            return p;
        }

        /**
         * Median normalized carrier signal
         */
        public double getCarrierSignal() {
            return carrierSignal;
        }

        /**
         * Median background carrier signal
         */
        public double getBackgroundSignal() {
            return backgroundSignal;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            PoissonTestResult that = (PoissonTestResult) o;
            return Double.compare(that.p, p) == 0 && Double.compare(that.carrierSignal, carrierSignal) == 0 && Double.compare(that.backgroundSignal, backgroundSignal) == 0;
        }

        @Override
        public int hashCode() {
            return Objects.hash(p, carrierSignal, backgroundSignal);
        }
    }
}
