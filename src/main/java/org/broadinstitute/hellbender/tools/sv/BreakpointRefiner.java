package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Tuple;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Refines variant breakpoints using split read evidence.
 *
 * The start and end of the breakpoint are tested independently. At each end we perform a series of Poisson
 * tests using the following model:
 *
 *  Given:
 *      d_i : depth of sample i
 *      c_i : split read count of sample i
 *
 *   Define a Poisson model of split read counts:
 *      r_i = c_i / d_i : normalized split read count of sample i
 *      m_c : median carrier r_i
 *      m_b : median background (non-carrier) r_i
 *      mu : mean depth of all samples
 *
 *      lambda = mu * m_c : expected carrier count
 *      X ~ Poisson(lambda) : carrier count model
 *      x_b = round(mu * m_b) : adjusted median background count
 *
 *   Calculate probability of observing the background count:
 *      p = P(X < x_b)
 *
 *   We then select the site with the lowest score. Breakpoint end positions are restricted by a lowerbound that depends
 *   on the refined start position (see {@link BreakpointRefiner#getEndLowerBound(SVCallRecord, int)}.
 */
public class BreakpointRefiner {

    private final Map<String,Double> sampleCoverageMap;
    private final SAMSequenceDictionary dictionary;
    /**
     * Number bases that insertion split read positions can pass by the original breakpoint, when left-clipped
     * reads are to the left of the breakpoint and right-clipped reads to the right.
     */
    private int maxInsertionSplitReadCrossDistance;

    public static final int DEFAULT_MAX_INSERTION_CROSS_DISTANCE = 20;

    /**
     * @param sampleCoverageMap map with (sample id, per-base sample coverage) entries
     * @param dictionary reference dictionary
     */
    public BreakpointRefiner(final Map<String,Double> sampleCoverageMap, final SAMSequenceDictionary dictionary) {
        this.sampleCoverageMap = Utils.nonNull(sampleCoverageMap);
        this.dictionary = Utils.nonNull(dictionary);
        setMaxInsertionSplitReadCrossDistance(DEFAULT_MAX_INSERTION_CROSS_DISTANCE);
    }

    public void setMaxInsertionSplitReadCrossDistance(final int distance) {
        maxInsertionSplitReadCrossDistance = distance;
    }

    /**
     * Performs breakend refinement for a call
     *
     * @param call with split read evidence
     * @return record with new breakpoints
     */
    public SVCallRecordWithEvidence refineCall(final SVCallRecordWithEvidence call) {
        Utils.nonNull(call);
        SVCallRecordUtils.validateCoordinatesWithDictionary(call, dictionary);

        // Depth-only calls cannot be refined
        if (SVClusterEngine.isDepthOnlyCall(call)) {
            return call;
        }

        // Sample sets
        final Set<String> backgroundSamples = getBackgroundSamples(call);
        final Set<String> calledSamples = call.getCalledSamples();

        // Refine start
        final SplitReadSite refinedStartSite = getRefinedSite(call.getStartSplitReadSites(), calledSamples, backgroundSamples, call.getPositionA());

        // Refine end
        final int endLowerBound = getEndLowerBound(call, refinedStartSite.getPosition());
        final int defaultEndPosition = Math.max(endLowerBound, call.getPositionB());
        final List<SplitReadSite> validEndSites = getValidEndSplitReadSites(call, endLowerBound);
        final SplitReadSite refinedEndSite = getRefinedSite(validEndSites, calledSamples, backgroundSamples, defaultEndPosition);

        final int length;
        if (call.getType().equals(StructuralVariantType.DEL) || call.getType().equals(StructuralVariantType.DEL)
                || call.getType().equals(StructuralVariantType.CNV) || call.getType().equals(StructuralVariantType.INV)) {
            length = call.getPositionB() - call.getPositionA() + 1;
        } else {
            length = call.getLength();
        }

        // Create new record
        return new SVCallRecordWithEvidence(
                call.getId(), call.getContigA(), refinedStartSite.getPosition(), call.getStrandA(), call.getContigB(),
                refinedEndSite.getPosition(), call.getStrandB(), call.getType(), length, call.getAlgorithms(),
                call.getGenotypes(), call.getStartSplitReadSites(), call.getEndSplitReadSites(), call.getDiscordantPairs(),
                call.getCopyNumberDistribution());
    }

    /**
     * Filters end sites with position less than lower-bound
     *
     * @param call
     * @param lowerBound min position
     * @return filtered set of end sites
     */
    private List<SplitReadSite> getValidEndSplitReadSites(final SVCallRecordWithEvidence call, final int lowerBound) {
        return call.getEndSplitReadSites().stream().filter(s -> s.getPosition() >= lowerBound).collect(Collectors.toList());
    }

    /**
     * Gets non-carrier samples
     *
     * @param call
     * @return sample ids
     */
    private Set<String> getBackgroundSamples(final SVCallRecord call) {
        return sampleCoverageMap.keySet().stream().filter(s -> !call.getCarrierSamples().contains(s)).collect(Collectors.toSet());
    }

    /**
     * Determines lower-bound on end site position (inclusive). For inter-chromosomal variants, boundaries are at the
     * start of the chromsome (any position is valid). For insertions, {@link BreakpointRefiner#maxInsertionSplitReadCrossDistance}
     * is used to determine how far past the original breakpoint it can be. Otherwise, we just use the new start position.
     *
     * @param call
     * @param refinedStartPosition new start position of call
     * @return position
     */
    private int getEndLowerBound(final SVCallRecord call, final int refinedStartPosition) {
        if (!SVCallRecordUtils.isIntrachromosomal(call)) {
            return 1;
        }
        if (call.getType().equals(StructuralVariantType.INS)) {
            return refinedStartPosition - maxInsertionSplitReadCrossDistance;
        }
        return refinedStartPosition + 1;
    }

    /**
     * Performs refinement on one side of a breakpoint
     *
     * @param sites sites to test
     * @param carrierSamples carrier sample ids
     * @param backgroundSamples background sample ids
     * @param defaultPosition position to use if test cannot be performed (no sites or no carriers)
     * @return site with refined breakpoints
     */
    private SplitReadSite getRefinedSite(final List<SplitReadSite> sites,
                                         final Set<String> carrierSamples,
                                         final Set<String> backgroundSamples,
                                         final int defaultPosition) {
        Utils.validateArg(sampleCoverageMap.keySet().containsAll(carrierSamples),
                "One or more carrier samples not found in sample coverage map");
        Utils.validateArg(sampleCoverageMap.keySet().containsAll(backgroundSamples),
                "One or more non-carrier samples not found in sample coverage map");

        // Default case
        if (sites.isEmpty() || carrierSamples.isEmpty()) return new SplitReadSite(defaultPosition, Collections.emptyMap());

        final long meanCoverage = Math.round(sampleCoverageMap.values().stream().mapToDouble(c -> c).sum()) / sampleCoverageMap.size();
        final List<Tuple<SplitReadSite,Double>> siteScorePairs = sites.stream()
                .map(s -> new Tuple<>(s, Double.valueOf(calculateOneSamplePoissonTest(s, carrierSamples, backgroundSamples, meanCoverage))))
                .collect(Collectors.toList());
        final double minLogP = siteScorePairs.stream().mapToDouble(p -> p.b).min().getAsDouble();
        return siteScorePairs.stream()
                .filter(p -> p.b == minLogP)
                .min(Comparator.comparingInt(s -> Math.abs(s.a.getPosition() - defaultPosition)))
                .get().a;
    }

    /**
     * Performs poisson test on a single site by computing the probability of observing the background counts
     * under a carrier count distribution
     *
     * @param site
     * @param carrierSamples
     * @param backgroundSamples
     * @param meanCoverage mean coverage of all samples
     * @return probability
     */
    private double calculateOneSamplePoissonTest(final SplitReadSite site,
                                                 final Set<String> carrierSamples,
                                                 final Set<String> backgroundSamples,
                                                 final double meanCoverage) {
        final double medianNormalizedCarrierCount = getMedianNormalizedCount(carrierSamples, site);
        if (medianNormalizedCarrierCount == 0) {
            return 1;  // Degenerate case in which the Poisson distribution is undefined
        }
        final double medianBackgroundRate = getMedianNormalizedCount(backgroundSamples, site);
        final int backgroundCount = (int) Math.round(medianBackgroundRate * meanCoverage);
        return cumulativePoissonProbability(meanCoverage * medianNormalizedCarrierCount, backgroundCount);
    }

    /**
     * Find the median of site counts in a subset of samples when normalized by sample coverage
     *
     * @param samples sample ids to restrict to
     * @param site
     * @return median
     */
    private double getMedianNormalizedCount(final Set<String> samples, final SplitReadSite site) {
        final Collection<Double> normalizedCounts = samples.stream()
                .map(s -> site.getCount(s) / sampleCoverageMap.get(s))
                .collect(Collectors.toList());
        return median(normalizedCounts);
    }

    private static double median(final Collection<Double> values) {
        if (values.isEmpty()) return 0.;
        final List<Double> sortedValues = values.stream().sorted().collect(Collectors.toList());
        if (sortedValues.size() % 2 == 1) {
            return sortedValues.get(sortedValues.size() /2 );
        }
        return 0.5 * (sortedValues.get(sortedValues.size() / 2) + sortedValues.get((sortedValues.size() / 2) - 1));
    }

    private static double cumulativePoissonProbability(final double mean, final int x) {
        if (x < 0) {
            return 0;
        }
        if (x == Integer.MAX_VALUE) {
            return 1;
        }
        return Gamma.regularizedGammaQ((double) x + 1, mean);
    }
}
