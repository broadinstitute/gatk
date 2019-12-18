package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.StructuralVariantType;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

public class BreakpointRefiner {

    private final Map<String,Double> sampleCoverageMap;
    private int maxInsertionSplitReadCrossDistance;

    public static final int DEFAULT_MAX_INSERTION_CROSS_DISTANCE = 20;

    public BreakpointRefiner(final Map<String,Double> sampleCoverageMap) {
        this.sampleCoverageMap = sampleCoverageMap;
        this.maxInsertionSplitReadCrossDistance = DEFAULT_MAX_INSERTION_CROSS_DISTANCE;
    }

    public void setMaxInsertionSplitReadCrossDistance(final int distance) {
        maxInsertionSplitReadCrossDistance = distance;
    }

    public SVCallRecordWithEvidence refineCall(final SVCallRecordWithEvidence call) {
        Utils.nonNull(call);
        Utils.nonNull(call.getStartSplitReadSites());
        final SVCallRecordWithEvidence refinedCall;
        if (SVClusterEngine.isDepthOnlyCall(call)) {
            refinedCall = call;
        } else {
            final Set<String> backgroundSamples = getBackgroundSamples(call);
            final SplitReadSite refinedStartSite = getRefinedSite(call.getStartSplitReadSites(), call.getSamples(), backgroundSamples, call.getStart());
            final int endLowerBound = getEndLowerBound(call.getType(), call.getContig(), refinedStartSite.getPosition(), call.getEndContig());
            final int defaultEndPosition = Math.max(endLowerBound, call.getEnd());
            final List<SplitReadSite> validSites = getValidEndSplitReadSites(call, endLowerBound);
            final SplitReadSite splitReadEndSite = getRefinedSite(validSites, call.getSamples(), backgroundSamples, defaultEndPosition);
            refinedCall = new SVCallRecordWithEvidence(
                    call.getContig(), refinedStartSite.getPosition(), call.getStartStrand(),
                    call.getEndContig(), splitReadEndSite.getPosition(), call.getEndStrand(),
                    call.getType(), call.getLength(), call.getAlgorithms(), call.getSamples(),
                    call.getStartSplitReadSites(), call.getEndSplitReadSites(), call.getDiscordantPairs());
        }
        return refinedCall;
    }

    private List<SplitReadSite> getValidEndSplitReadSites(final SVCallRecordWithEvidence call, final int endLowerBound) {
        return call.getEndSplitReadSites().stream()
                .filter(s -> s.getPosition() >= endLowerBound)
                .collect(Collectors.toList());
    }

    private Set<String> getBackgroundSamples(final SVCallRecord call) {
        return sampleCoverageMap.keySet().stream().filter(s -> !call.getSamples().contains(s)).collect(Collectors.toSet());
    }

    private int getEndLowerBound(final StructuralVariantType type, final String startContig, final int startPosition, final String endContig) {
        if (!startContig.equals(endContig)) return 0;
        return type.equals(StructuralVariantType.INS) ?
                startPosition - maxInsertionSplitReadCrossDistance :
                startPosition + 1;
    }

    private SplitReadSite getRefinedSite(final List<SplitReadSite> sites,
                                         final Set<String> carrierSamples,
                                         final Set<String> backgroundSamples,
                                         final int defaultPosition) {
        if (!sampleCoverageMap.keySet().containsAll(carrierSamples)) {
            throw new IllegalArgumentException("One or more carrier samples not found in sample coverage map");
        }
        if (!sampleCoverageMap.keySet().containsAll(backgroundSamples)) {
            throw new IllegalArgumentException("One or more non-carrier samples not found in sample coverage map");
        }
        return testSitesByOneSamplePoissonTest(sites, defaultPosition, carrierSamples, backgroundSamples);
    }

    private SplitReadSite testSitesByOneSamplePoissonTest(final List<SplitReadSite> sites,
                                                          final int defaultPosition,
                                                          final Set<String> carrierSamples,
                                                          final Set<String> backgroundSamples) {
        if (sites.isEmpty() || carrierSamples.isEmpty()) return new SplitReadSite(defaultPosition, Collections.emptyMap());
        final long meanCoverage = Math.round(sampleCoverageMap.values().stream().mapToDouble(c -> c).sum()) / sampleCoverageMap.size();
        final List<Tuple2<SplitReadSite,Double>> siteScorePairs = sites.stream()
                .map(s -> new Tuple2<>(s, calculateOneSamplePoissonTest(s, carrierSamples, backgroundSamples, meanCoverage)))
                .collect(Collectors.toList());
        final double maxLogP = siteScorePairs.stream().mapToDouble(Tuple2::_2).max().getAsDouble();
        return siteScorePairs.stream()
                .filter(p -> p._2 == maxLogP)
                .min(Comparator.comparingInt(s -> Math.abs(s._1.getPosition() - defaultPosition)))
                .get()._1;
    }

    private double calculateOneSamplePoissonTest(final SplitReadSite site,
                                                 final Set<String> carrierSamples,
                                                 final Set<String> backgroundSamples,
                                                 final double meanCoverage) {
        final double medianNormalizedCarrierCount = getMedianNormalizedCount(carrierSamples, site);
        if (medianNormalizedCarrierCount == 0) {
            return 0;
        }
        final double medianBackgroundRate = getMedianNormalizedCount(backgroundSamples, site);
        final int backgroundCount = (int) Math.round(medianBackgroundRate * meanCoverage);
        return 1.0 - cumulativePoissonProbability(meanCoverage * medianNormalizedCarrierCount, backgroundCount);
    }

    private double getMedianNormalizedCount(final Set<String> samples, final SplitReadSite site) {
        final Collection<Double> normalizedCounts = samples.stream()
                .map(s -> site.getSampleCountsMap().containsKey(s) ? site.getSampleCountsMap().get(s) / sampleCoverageMap.get(s) : 0.)
                .collect(Collectors.toList());
        return median(normalizedCounts);
    }

    private static double median(final Collection<Double> values) {
        if (values.isEmpty()) return 0.;
        final List<Double> sortedValues = values.stream().sorted().collect(Collectors.toList());
        if (sortedValues.size() % 2 == 1) {
            return sortedValues.get(sortedValues.size() /2 );
        }
        return (sortedValues.get(sortedValues.size() / 2) + sortedValues.get((sortedValues.size() / 2) - 1)) / 2.;
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
}
