package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Streams;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.utils.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public abstract class CachingSVEvidenceAggregator<T extends Feature> {

    private final FeatureDataSource<T> source;
    private SimpleInterval cacheInterval;
    private IntervalTree<List<T>> cacheEvidence;
    private final ProgressMeter progressMeter;
    protected final SAMSequenceDictionary dictionary;

    private static final long RECORDS_BETWEEN_TIME_CHECKS = 100L;

    public CachingSVEvidenceAggregator(final FeatureDataSource<T> source,
                                       final SAMSequenceDictionary dictionary,
                                       final String progressLabel) {
        this.source = source;
        this.dictionary = dictionary;
        this.cacheInterval = null;
        this.cacheEvidence = null;
        this.progressMeter = new ProgressMeter();
        progressMeter.setRecordLabel(progressLabel);
        progressMeter.setRecordsBetweenTimeChecks(RECORDS_BETWEEN_TIME_CHECKS);
    }

    abstract protected SimpleInterval getEvidenceQueryInterval(final SVCallRecordWithEvidence record);
    abstract protected SVCallRecordWithEvidence assignEvidence(final SVCallRecordWithEvidence call, final List<T> evidence);
    protected boolean evidenceFilter(final SVCallRecord record, final T evidence) { return true; }

    public Stream<SVCallRecordWithEvidence> collectEvidence(final List<SVCallRecordWithEvidence> calls) {
        Utils.nonNull(calls);
        final OverlapDetector<SimpleInterval> overlapDetector = getEvidenceOverlapDetector(calls);
        return calls.stream().map(r -> collectRecordEvidence(r, overlapDetector));
    }

    public void startProgressMeter() {
        progressMeter.start();
    }

    public void stopProgressMeter() {
        progressMeter.stop();
    }

    private final SVCallRecordWithEvidence collectRecordEvidence(final SVCallRecordWithEvidence record,
                                                                 final OverlapDetector<SimpleInterval> overlapDetector) {
        final SimpleInterval evidenceInterval = getEvidenceQueryInterval(record);
        final List<T> evidence = getEvidenceOnInterval(evidenceInterval, overlapDetector).stream()
                .filter(e -> evidenceFilter(record, e))
                .collect(Collectors.toList());
        if (progressMeter.started()) {
            progressMeter.update(evidenceInterval);
        }
        return assignEvidence(record, evidence);
    }

    private final List<T> getEvidenceOnInterval(final SimpleInterval interval,
                                                final OverlapDetector<SimpleInterval> overlapDetector) {
        if (invalidCacheInterval(cacheInterval, interval)) {
            Utils.nonNull(overlapDetector, "Evidence cache missed but overlap detector is null");
            final Set<SimpleInterval> queryIntervalSet = overlapDetector.getOverlaps(interval);
            if (queryIntervalSet.size() != 1) {
                throw new IllegalArgumentException("Evidence interval " + interval + " overlapped " + queryIntervalSet.size() + " query intervals");
            }
            cacheInterval = queryIntervalSet.iterator().next();
            cacheEvidence = new IntervalTree<>();
            source.queryAndPrefetch(cacheInterval).stream().forEachOrdered(this::addItemToCacheTree);
        } else {
            cacheEvidence.remove(0, interval.getStart() - 1);
        }
        return Streams.stream(cacheEvidence.overlappers(interval.getStart(), interval.getEnd()))
                .map(IntervalTree.Node::getValue).flatMap(List::stream).collect(Collectors.toList()); //source.queryAndPrefetch(interval);
    }

    private void addItemToCacheTree(final T item) {
        final IntervalTree.Node<List<T>> overlapper = cacheEvidence.find(item.getStart(), item.getEnd());
        if (overlapper == null) {
            cacheEvidence.put(item.getStart(), item.getEnd(), Collections.singletonList(item));
        } else {
            final List<T> currentList = overlapper.getValue();
            final List<T> newList = new ArrayList<>(currentList.size());
            newList.addAll(currentList);
            newList.add(item);
            overlapper.setValue(newList);
        }
    }

    private OverlapDetector<SimpleInterval> getEvidenceOverlapDetector(final List<SVCallRecordWithEvidence> calls) {
        final List<SimpleInterval> rawIntervals = calls.stream()
                .map(this::getEvidenceQueryInterval)
                .sorted(IntervalUtils.getDictionaryOrderComparator(dictionary))
                .collect(Collectors.toList());
        final GenomeLocParser parser = new GenomeLocParser(dictionary);
        final List<GenomeLoc> rawLocs = IntervalUtils.genomeLocsFromLocatables(parser, rawIntervals);
        final List<GenomeLoc> mergedLocs = IntervalUtils.mergeIntervalLocations(rawLocs, IntervalMergingRule.ALL);
        final List<SimpleInterval> mergedIntervals = IntervalUtils.convertGenomeLocsToSimpleIntervals(mergedLocs);
        return OverlapDetector.create(mergedIntervals);
    }

    private final boolean invalidCacheInterval(final SimpleInterval cacheInterval, final SimpleInterval queryInterval) {
        return cacheInterval == null || !cacheInterval.contains(queryInterval);
    }
}
