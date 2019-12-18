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

public class SVEvidenceCollector {

    private final FeatureDataSource<SplitReadEvidence> splitReadSource;
    private final FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private final SAMSequenceDictionary dictionary;

    private SimpleInterval splitReadCacheInterval;
    private SimpleInterval discordantPairCacheInterval;
    private int splitReadPadding;
    private int discordantPairPadding;
    private ProgressMeter progressMeter;

    public static final int DEFAULT_SPLIT_READ_PADDING = 50;
    public static final int DEFAULT_DISCORDANT_PAIR_PADDING = 500;

    public SVEvidenceCollector(final FeatureDataSource<SplitReadEvidence> splitReadSource,
                               final FeatureDataSource<DiscordantPairEvidence> discordantPairSource,
                               final SAMSequenceDictionary dictionary,
                               final ProgressMeter progressMeter) {
        this.splitReadSource = splitReadSource;
        this.discordantPairSource = discordantPairSource;
        this.dictionary = dictionary;
        this.splitReadCacheInterval = null;
        this.discordantPairCacheInterval = null;
        this.splitReadPadding = DEFAULT_SPLIT_READ_PADDING;
        this.discordantPairPadding = DEFAULT_DISCORDANT_PAIR_PADDING;
        this.progressMeter = progressMeter;
    }

    public SVEvidenceCollector(final FeatureDataSource<SplitReadEvidence> splitReadSource,
                               final FeatureDataSource<DiscordantPairEvidence> discordantPairSource,
                               final SAMSequenceDictionary dictionary) {
        this(splitReadSource, discordantPairSource, dictionary, null);
    }

    public void setSplitReadPadding(final int padding) {
        splitReadPadding = padding;
    }

    public void setDiscordantPairPadding(final int padding) {
        discordantPairPadding = padding;
    }

    public int getSplitReadPadding() {
        return splitReadPadding;
    }

    public int getDiscordantPairPadding() {
        return discordantPairPadding;
    }

    public List<SVCallRecordWithEvidence> collectEvidence(final List<SVCallRecordWithEvidence> calls) {
        Utils.nonNull(calls);
        return processEndPositions(processStartPositions(calls));
    }

    private List<SVCallRecordWithEvidence> processStartPositions(final List<SVCallRecordWithEvidence> calls) {
        final OverlapDetector splitReadStartIntervalOverlapDetector = getEvidenceOverlapDetector(calls, this::getStartSplitReadInterval);
        final OverlapDetector discordantPairIntervalOverlapDetector = getEvidenceOverlapDetector(calls, this::getDiscordantPairStartInterval);
        return calls.stream()
                .map(c -> processDiscordantPairs(c, discordantPairIntervalOverlapDetector))
                .collect(Collectors.toList())
                .stream()
                .map(c -> processStartSplitReads(c, splitReadStartIntervalOverlapDetector))
                .collect(Collectors.toList());
    }

    private List<SVCallRecordWithEvidence> processEndPositions(final List<SVCallRecordWithEvidence> calls) {
        final OverlapDetector splitReadEndIntervalOverlapDetector = getEvidenceOverlapDetector(calls, this::getEndSplitReadInterval);
        return calls.stream().sorted(Comparator.comparing(c -> c.getEndAsInterval(), IntervalUtils.getDictionaryOrderComparator(dictionary)))
                .map(c -> processEndSplitReads(c, splitReadEndIntervalOverlapDetector))
                .collect(Collectors.toList());
    }

    private <T extends SVCallRecord> OverlapDetector getEvidenceOverlapDetector(final List<T> calls,
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

    private SVCallRecordWithEvidence processDiscordantPairs(final SVCallRecordWithEvidence call,
                                                            final OverlapDetector discordantPairIntervalOverlapDetector) {
        Utils.nonNull(call);
        final SVCallRecordWithEvidence callWithEvidence;
        if (SVClusterEngine.isDepthOnlyCall(call)) {
            callWithEvidence = call;
        } else {
            final List<DiscordantPairEvidence> discordantPairs = getDiscordantPairs(call, discordantPairIntervalOverlapDetector);
            callWithEvidence = new SVCallRecordWithEvidence(
                    call.getContig(), call.getStart(), call.getStartStrand(), call.getEndContig(), call.getEnd(), call.getEndStrand(),
                    call.getType(), call.getLength(), call.getAlgorithms(), call.getSamples(), call.getStartSplitReadSites(), call.getEndSplitReadSites(), discordantPairs);
        }
        if (progressMeter != null) {
            progressMeter.update(call.getStartAsInterval());
        }
        return callWithEvidence;
    }

    private SVCallRecordWithEvidence processStartSplitReads(final SVCallRecordWithEvidence call,
                                                             final OverlapDetector splitReadStartIntervalOverlapDetector) {
        Utils.nonNull(call);
        final SVCallRecordWithEvidence callWithEvidence;
        if (SVClusterEngine.isDepthOnlyCall(call)) {
            callWithEvidence = call;
        } else {
            final List<SplitReadEvidence> startSplitReads = getStartSplitReads(call, splitReadStartIntervalOverlapDetector);
            final List<SplitReadSite> startSitesList = computeSites(startSplitReads, call.getStartStrand());
            callWithEvidence = new SVCallRecordWithEvidence(
                    call.getContig(), call.getStart(), call.getStartStrand(), call.getEndContig(), call.getEnd(), call.getEndStrand(),
                    call.getType(), call.getLength(), call.getAlgorithms(), call.getSamples(), startSitesList, call.getEndSplitReadSites(), call.getDiscordantPairs());
        }
        if (progressMeter != null) {
            progressMeter.update(call.getStartAsInterval());
        }
        return callWithEvidence;
    }

    private SVCallRecordWithEvidence processEndSplitReads(final SVCallRecordWithEvidence call,
                                                           final OverlapDetector splitReadEndIntervalOverlapDetector) {
        Utils.nonNull(call);
        final SVCallRecordWithEvidence refinedCall;
        if (SVClusterEngine.isDepthOnlyCall(call)) {
            refinedCall = call;
        } else {
            final List<SplitReadEvidence> endSplitReads = getEndSplitReads(call, splitReadEndIntervalOverlapDetector);
            final List<SplitReadSite> endSitesList = computeSites(endSplitReads, call.getEndStrand());
            refinedCall = new SVCallRecordWithEvidence(
                    call.getContig(), call.getStart(), call.getStartStrand(), call.getEndContig(), call.getEnd(), call.getEndStrand(),
                    call.getType(), call.getLength(), call.getAlgorithms(), call.getSamples(), call.getStartSplitReadSites(), endSitesList, call.getDiscordantPairs());
        }
        if (progressMeter != null) {
            progressMeter.update(call.getEndAsInterval());
        }
        return refinedCall;
    }

    private List<SplitReadEvidence> getSplitReads(final SimpleInterval interval,
                                                  final boolean strand,
                                                  final OverlapDetector splitReadOverlapDetector) {
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
                                                       final OverlapDetector splitReadStartOverlapDetector) {
        return getSplitReads(getStartSplitReadInterval(call), call.getStartStrand(), splitReadStartOverlapDetector);
    }

    private List<SplitReadEvidence> getEndSplitReads(final SVCallRecordWithEvidence call,
                                                     final OverlapDetector splitReadEndOverlapDetector) {
        return getSplitReads(getEndSplitReadInterval(call), call.getEndStrand(), splitReadEndOverlapDetector);
    }

    private SimpleInterval getStartSplitReadInterval(final SVCallRecord call) {
        return call.getStartAsInterval().expandWithinContig(splitReadPadding, dictionary);
    }

    private SimpleInterval getEndSplitReadInterval(final SVCallRecord call) {
        return call.getEndAsInterval().expandWithinContig(splitReadPadding, dictionary);
    }

    private boolean invalidCacheInterval(final SimpleInterval cacheInterval, final SimpleInterval queryInterval) {
        return cacheInterval == null
                || !queryInterval.getContig().equals(cacheInterval.getContig())
                || !queryInterval.spanWith(cacheInterval).equals(cacheInterval);
    }

    private List<DiscordantPairEvidence> getDiscordantPairs(final SVCallRecord call,
                                                            final OverlapDetector discordantPairStartOverlapDetector) {
        final SimpleInterval startInterval = getDiscordantPairStartInterval(call);
        if (invalidCacheInterval(discordantPairCacheInterval, startInterval)) {
            final Set<SimpleInterval> queryIntervalSet = discordantPairStartOverlapDetector.getOverlaps(startInterval);
            if (queryIntervalSet.size() != 1) {
                throw new IllegalArgumentException("Call end split read interval " + startInterval + " overlapped " + queryIntervalSet.size() + " query intervals");
            }
            discordantPairCacheInterval = queryIntervalSet.iterator().next();
            discordantPairSource.queryAndPrefetch(discordantPairCacheInterval);
        }
        final SimpleInterval endInterval = getDiscordantPairEndInterval(call);
        return discordantPairSource.queryAndPrefetch(startInterval).stream()
                .filter(e -> discordantPairOverlapsInterval(e, startInterval, endInterval))
                .filter(e -> e.getStartStrand() == call.getStartStrand() && e.getEndStrand() == call.getEndStrand())
                .collect(Collectors.toList());
    }

    private SimpleInterval getDiscordantPairStartInterval(final SVCallRecord call) {
        return call.getStartAsInterval().expandWithinContig(discordantPairPadding, dictionary);
    }

    private SimpleInterval getDiscordantPairEndInterval(final SVCallRecord call) {
        return call.getEndAsInterval().expandWithinContig(discordantPairPadding, dictionary);
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
            if (e.getStrand() == strand) {
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
