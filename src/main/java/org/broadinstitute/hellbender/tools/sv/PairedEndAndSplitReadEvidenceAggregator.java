package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.utils.*;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class PairedEndAndSplitReadEvidenceAggregator {

    private final FeatureDataSource<SplitReadEvidence> splitReadSource;
    private final FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private final SAMSequenceDictionary dictionary;

    private SimpleInterval splitReadCacheInterval;
    private SimpleInterval discordantPairCacheInterval;
    private int splitReadPadding;
    private int discordantPairInnerPadding;
    private int discordantPairOuterPadding;
    private ProgressMeter progressMeter;

    public static final int DEFAULT_SPLIT_READ_PADDING = 50;
    public static final int DEFAULT_DISCORDANT_PAIR_INNER_PADDING = 50;
    public static final int DEFAULT_DISCORDANT_PAIR_OUTER_PADDING = 500;

    public PairedEndAndSplitReadEvidenceAggregator(final FeatureDataSource<SplitReadEvidence> splitReadSource,
                                                   final FeatureDataSource<DiscordantPairEvidence> discordantPairSource,
                                                   final SAMSequenceDictionary dictionary,
                                                   final ProgressMeter progressMeter) {
        this.splitReadSource = splitReadSource;
        this.discordantPairSource = discordantPairSource;
        this.dictionary = dictionary;
        this.splitReadCacheInterval = null;
        this.discordantPairCacheInterval = null;
        this.splitReadPadding = DEFAULT_SPLIT_READ_PADDING;
        this.discordantPairInnerPadding = DEFAULT_DISCORDANT_PAIR_INNER_PADDING;
        this.discordantPairOuterPadding = DEFAULT_DISCORDANT_PAIR_OUTER_PADDING;
        this.progressMeter = progressMeter;
    }

    public PairedEndAndSplitReadEvidenceAggregator(final FeatureDataSource<SplitReadEvidence> splitReadSource,
                                                   final FeatureDataSource<DiscordantPairEvidence> discordantPairSource,
                                                   final SAMSequenceDictionary dictionary) {
        this(splitReadSource, discordantPairSource, dictionary, null);
    }

    public void setSplitReadPadding(final int padding) {
        splitReadPadding = padding;
    }

    public void setDiscordantPairInnerPadding(final int padding) {
        discordantPairInnerPadding = padding;
    }

    public void setDiscordantPairOuterPadding(final int padding) {
        discordantPairOuterPadding = padding;
    }

    public int getSplitReadPadding() {
        return splitReadPadding;
    }

    public int getDiscordantPairInnerPadding() {
        return discordantPairInnerPadding;
    }

    public int getDiscordantPairOuterPadding() {
        return discordantPairOuterPadding;
    }

    public List<SVCallRecordWithEvidence> collectEvidence(final List<SVCallRecord> calls) {
        Utils.nonNull(calls);
        return processEndPositions(processStartPositions(calls));
    }

    private List<SVCallRecordWithEvidence> processStartPositions(final List<SVCallRecord> calls) {
        final OverlapDetector<SimpleInterval> splitReadStartIntervalOverlapDetector = getEvidenceOverlapDetector(calls, this::getStartSplitReadInterval);
        final OverlapDetector<SimpleInterval> discordantPairIntervalOverlapDetector = getEvidenceOverlapDetector(calls, this::getDiscordantPairStartInterval);
        return calls.stream()
                .map(c -> processDiscordantPairs(c, discordantPairIntervalOverlapDetector))
                .collect(Collectors.toList())
                .stream()
                .map(c -> processStartSplitReads(c, splitReadStartIntervalOverlapDetector))
                .collect(Collectors.toList());
    }

    private List<SVCallRecordWithEvidence> processEndPositions(final List<SVCallRecordWithEvidence> calls) {
        final OverlapDetector<SimpleInterval> splitReadEndIntervalOverlapDetector = getEvidenceOverlapDetector(calls, this::getEndSplitReadInterval);
        return calls.stream().sorted(Comparator.comparing(c -> c.getPositionBInterval(), IntervalUtils.getDictionaryOrderComparator(dictionary)))
                .map(c -> processEndSplitReads(c, splitReadEndIntervalOverlapDetector))
                .collect(Collectors.toList());
    }

    private <T extends SVCallRecord> OverlapDetector<SimpleInterval> getEvidenceOverlapDetector(final List<T> calls,
                                                                                final Function<T, SimpleInterval> evidenceFunction) {
        final List<SimpleInterval> rawIntervals = calls.stream()
                .map(c -> evidenceFunction.apply(c))
                .sorted(IntervalUtils.getDictionaryOrderComparator(dictionary))
                .collect(Collectors.toList());
        final GenomeLocParser parser = new GenomeLocParser(dictionary);
        final List<GenomeLoc> rawLocs = IntervalUtils.genomeLocsFromLocatables(parser, rawIntervals);
        final List<GenomeLoc> mergedLocs = IntervalUtils.mergeIntervalLocations(rawLocs, IntervalMergingRule.ALL);
        final List<SimpleInterval> mergedIntervals = IntervalUtils.convertGenomeLocsToSimpleIntervals(mergedLocs);
        return OverlapDetector.create(mergedIntervals);
    }

    private SVCallRecordWithEvidence processDiscordantPairs(final SVCallRecord call,
                                                            final OverlapDetector<SimpleInterval> discordantPairIntervalOverlapDetector) {
        Utils.nonNull(call);
        final SVCallRecordWithEvidence callWithEvidence;
        if (SVClusterEngine.isDepthOnlyCall(call)) {
            callWithEvidence = new SVCallRecordWithEvidence(call, Collections.emptyList(), Collections.emptyList(), Collections.emptyList(), null);
        } else {
            final List<DiscordantPairEvidence> discordantPairs = getDiscordantPairs(call, discordantPairIntervalOverlapDetector);
            callWithEvidence = new SVCallRecordWithEvidence(call, Collections.emptyList(), Collections.emptyList(), discordantPairs, null);
        }
        if (progressMeter != null) {
            progressMeter.update(call.getPositionAInterval());
        }
        return callWithEvidence;
    }

    private SVCallRecordWithEvidence processStartSplitReads(final SVCallRecordWithEvidence call,
                                                            final OverlapDetector<SimpleInterval> splitReadStartIntervalOverlapDetector) {
        Utils.nonNull(call);
        final SVCallRecordWithEvidence callWithEvidence;
        if (SVClusterEngine.isDepthOnlyCall(call)) {
            callWithEvidence = call;
        } else {
            final List<SplitReadEvidence> startSplitReads = getStartSplitReads(call, splitReadStartIntervalOverlapDetector);
            final List<SplitReadSite> startSitesList = computeSites(startSplitReads, call.getStrandA());
            callWithEvidence = new SVCallRecordWithEvidence(call, startSitesList, call.getEndSplitReadSites(), call.getDiscordantPairs(), call.getCopyNumberDistribution());
        }
        if (progressMeter != null) {
            progressMeter.update(call.getPositionAInterval());
        }
        return callWithEvidence;
    }

    private SVCallRecordWithEvidence processEndSplitReads(final SVCallRecordWithEvidence call,
                                                          final OverlapDetector<SimpleInterval> splitReadEndIntervalOverlapDetector) {
        Utils.nonNull(call);
        final SVCallRecordWithEvidence refinedCall;
        if (SVClusterEngine.isDepthOnlyCall(call)) {
            refinedCall = call;
        } else {
            final List<SplitReadEvidence> endSplitReads = getEndSplitReads(call, splitReadEndIntervalOverlapDetector);
            final List<SplitReadSite> endSitesList = computeSites(endSplitReads, call.getStrandB());
            refinedCall = new SVCallRecordWithEvidence(call, call.getStartSplitReadSites(), endSitesList, call.getDiscordantPairs(), call.getCopyNumberDistribution());
        }
        if (progressMeter != null) {
            progressMeter.update(call.getPositionBInterval());
        }
        return refinedCall;
    }

    private List<SplitReadEvidence> getSplitReads(final SimpleInterval interval,
                                                  final boolean strand,
                                                  final OverlapDetector<SimpleInterval> splitReadOverlapDetector) {
        if (invalidCacheInterval(splitReadCacheInterval, interval)) {
            Utils.nonNull(splitReadOverlapDetector, "Split read cache missed but overlap detector is null");
            final Set<SimpleInterval> queryIntervalSet = splitReadOverlapDetector.getOverlaps(interval);
            if (queryIntervalSet.size() != 1) {
                throw new IllegalArgumentException("Call split read interval " + interval + " overlapped " + queryIntervalSet.size() + " query intervals");
            }
            splitReadCacheInterval = queryIntervalSet.iterator().next();
            splitReadSource.queryAndPrefetch(splitReadCacheInterval);
        }
        return splitReadSource.queryAndPrefetch(interval).stream()
                .filter(e -> e.getStrand() == strand)
                .collect(Collectors.toList());
    }

    private List<SplitReadEvidence> getStartSplitReads(final SVCallRecordWithEvidence call,
                                                       final OverlapDetector<SimpleInterval> splitReadStartOverlapDetector) {
        return getSplitReads(getStartSplitReadInterval(call), call.getStrandA(), splitReadStartOverlapDetector);
    }

    private List<SplitReadEvidence> getEndSplitReads(final SVCallRecordWithEvidence call,
                                                     final OverlapDetector<SimpleInterval> splitReadEndOverlapDetector) {
        return getSplitReads(getEndSplitReadInterval(call), call.getStrandB(), splitReadEndOverlapDetector);
    }

    private SimpleInterval getStartSplitReadInterval(final SVCallRecord call) {
        return call.getPositionAInterval().expandWithinContig(splitReadPadding, dictionary);
    }

    private SimpleInterval getEndSplitReadInterval(final SVCallRecord call) {
        return call.getPositionBInterval().expandWithinContig(splitReadPadding, dictionary);
    }

    private boolean invalidCacheInterval(final SimpleInterval cacheInterval, final SimpleInterval queryInterval) {
        return cacheInterval == null
                || !queryInterval.getContig().equals(cacheInterval.getContig())
                || !queryInterval.spanWith(cacheInterval).equals(cacheInterval);
    }

    private List<DiscordantPairEvidence> getDiscordantPairs(final SVCallRecord call,
                                                            final OverlapDetector<SimpleInterval> discordantPairStartOverlapDetector) {
        final SimpleInterval startInterval = getDiscordantPairStartInterval(call);
        if (invalidCacheInterval(discordantPairCacheInterval, startInterval)) {
            final Set<SimpleInterval> queryIntervalSet = discordantPairStartOverlapDetector.getOverlaps(startInterval);
            if (queryIntervalSet.size() != 1) {
                throw new IllegalArgumentException("Call discordant pair interval " + startInterval + " overlapped " + queryIntervalSet.size() + " query intervals");
            }
            discordantPairCacheInterval = queryIntervalSet.iterator().next();
            discordantPairSource.queryAndPrefetch(discordantPairCacheInterval);
        }
        final SimpleInterval endInterval = getDiscordantPairEndInterval(call);
        return discordantPairSource.queryAndPrefetch(startInterval).stream()
                .filter(e -> discordantPairOverlapsInterval(e, startInterval, endInterval))
                .filter(e -> e.getStartStrand() == call.getStrandA() && e.getEndStrand() == call.getStrandB())
                .collect(Collectors.toList());
    }

    private SimpleInterval getDiscordantPairStartInterval(final SVCallRecord call) {
        final String contig = call.getContigA();
        if (call.getStrandA()) {
            return IntervalUtils.trimIntervalToContig(contig, call.getPositionA() - discordantPairOuterPadding, call.getPositionA() + discordantPairInnerPadding, dictionary.getSequence(contig).getSequenceLength());
        } else {
            return IntervalUtils.trimIntervalToContig(contig, call.getPositionA() - discordantPairInnerPadding, call.getPositionA() + discordantPairOuterPadding, dictionary.getSequence(contig).getSequenceLength());
        }
    }

    private SimpleInterval getDiscordantPairEndInterval(final SVCallRecord call) {
        final String contig = call.getContigB();
        if (call.getStrandB()) {
            return IntervalUtils.trimIntervalToContig(contig, call.getPositionB() - discordantPairOuterPadding, call.getPositionB() + discordantPairInnerPadding, dictionary.getSequence(contig).getSequenceLength());
        } else {
            return IntervalUtils.trimIntervalToContig(contig, call.getPositionB() - discordantPairInnerPadding, call.getPositionB() + discordantPairOuterPadding, dictionary.getSequence(contig).getSequenceLength());
        }
    }

    private boolean discordantPairOverlapsInterval(final DiscordantPairEvidence evidence,
                                                   final SimpleInterval startInterval,
                                                   final SimpleInterval endInterval) {
        return evidence.getContig().equals(startInterval.getContig())
                && evidence.getStart() >= startInterval.getStart()
                && evidence.getStart() < startInterval.getEnd()
                && evidence.getEndContig().equals(endInterval.getContig())
                && evidence.getEnd() >= endInterval.getStart()
                && evidence.getEnd() < endInterval.getEnd();
    }

    private List<SplitReadSite> computeSites(final List<SplitReadEvidence> evidenceList, final boolean strand) {
        if (!Ordering.from(IntervalUtils.getDictionaryOrderComparator(dictionary)).isOrdered(evidenceList)) {
            throw new IllegalArgumentException("Evidence list is not dictionary sorted");
        }
        final ArrayList<SplitReadSite> sites = new ArrayList<>();
        int position = 0;
        Map<String,Integer> sampleCounts = new HashMap<>();
        for (final SplitReadEvidence e : evidenceList) {
            if (e.getStart() != position) {
                if (!sampleCounts.isEmpty()) {
                    sites.add(new SplitReadSite(position, sampleCounts));
                    sampleCounts = new HashMap<>();
                }
                position = e.getStart();
            }
            if (e.getStrand() == strand && e.getCount() > 0) {
                final String sample = e.getSample();
                sampleCounts.put(sample, e.getCount());
            }
        }
        if (!sampleCounts.isEmpty()) {
            sites.add(new SplitReadSite(position, sampleCounts));
        }
        sites.trimToSize();
        return sites;
    }
}
