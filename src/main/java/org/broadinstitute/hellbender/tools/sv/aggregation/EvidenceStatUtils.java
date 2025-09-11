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

    public static Double probToQual(final Double pError, final double maxQual) {
        return pError == null ? null : Math.min(-10.0*Math.log10(pError), maxQual);
    }

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
     * Find the median of site counts in a subset of samples when normalized by sample coverage
     *
     * @param samples sample ids to restrict to
     * @param sampleCounts
     * @return median
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
            normalizedCounts[i] = Math.round(targetCoverage * sampleCounts.getOrDefault(sample, 0) / sampleCoverageMap.get(sample));
            i++;
        }
        return MEDIAN.evaluate(normalizedCounts);
    }

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
     * @param sampleCounts
     * @param carrierSamples
     * @param backgroundSamples
     * @param meanCoverage mean coverage of all samples
     * @return probability on [0,1]
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

        public double getP() {
            return p;
        }

        public double getCarrierSignal() {
            return carrierSignal;
        }

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
