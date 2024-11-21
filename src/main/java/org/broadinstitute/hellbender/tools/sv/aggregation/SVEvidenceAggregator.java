package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVFeature;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public abstract class SVEvidenceAggregator<T extends SVFeature> {

    private final FeatureDataSource<T> source;
    private SimpleInterval cacheInterval;
    private Deque<T> cacheEvidence;
    private OverlapDetector<SimpleInterval> cacheIntervalTree;
    private final int maxWindowSize;  // Evidence types like PE have different window sizes, should be largest
    protected final SAMSequenceDictionary dictionary;

    public SVEvidenceAggregator(final FeatureDataSource<T> source,
                                final int maxWindowSize,
                                final SAMSequenceDictionary dictionary) {
        Utils.nonNull(source);
        Utils.nonNull(dictionary);
        this.source = source;
        this.maxWindowSize = maxWindowSize;
        this.dictionary = dictionary;
    }

    public void setCacheIntervals(final Collection<SimpleInterval> intervals) {
        cacheIntervalTree = OverlapDetector.create(
                IntervalUtils.sortAndMergeIntervalsToStream(intervals, dictionary, IntervalMergingRule.ALL)
                        .collect(Collectors.toList())
        );
    }

    abstract public SimpleInterval getEvidenceQueryInterval(final SVCallRecord record);
    abstract public boolean evidenceFilter(final SVCallRecord record, final T evidence);

    private SimpleInterval getRegionInterval(final SimpleInterval interval) {
        final Set<SimpleInterval> evidenceIntervals = cacheIntervalTree.getOverlaps(interval);
        Utils.validate(evidenceIntervals.size() == 1, "Expected exactly 1 evidence interval but " +
                "found " + evidenceIntervals.size());
        return evidenceIntervals.iterator().next();
    }

    public List<T> collectEvidence(final SVCallRecord call) {
        Utils.nonNull(call);
        final SimpleInterval callInterval = getEvidenceQueryInterval(call);
        if (callInterval == null) {
            return Collections.emptyList();
        }
        final Collection<T> rawEvidence;
        if (cacheIntervalTree == null) {
            rawEvidence = source.queryAndPrefetch(callInterval);
        } else {
            final SimpleInterval regionInterval = getRegionInterval(callInterval);
            if (!regionInterval.equals(cacheInterval)) {
                cacheEvidence = new ArrayDeque<>(source.queryAndPrefetch(regionInterval));
                cacheInterval = regionInterval;
            }
            // Expect to encounter variants in sorted order, but window size may vary and cause the start positions
            // to be unsorted. For example, PE has inner/outer window that depends on orientation.
            while (!cacheEvidence.isEmpty() && callInterval.getStart() - maxWindowSize > cacheEvidence.peek().getStart()) {
                cacheEvidence.pop();
            }
            rawEvidence = cacheEvidence;
        }
        final List<T> callEvidence = new ArrayList<>();
        boolean foundOverlap = false;
        for (final T evidence : rawEvidence) {
            if (callInterval.overlaps(evidence)) {
                foundOverlap = true;
                if (evidenceFilter(call, evidence)) {
                    callEvidence.add(evidence);
                }
            } else if (foundOverlap) {
                break;
            }
        }
        return callEvidence;
    }
}
