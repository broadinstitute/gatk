package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Common utility functions for SVInterval and SVIntervalTree objects
 */
public class SVIntervalUtils {

    /**
     * Given a link with target intervals on the same contig, returns the interval from the left target interval's
     * start to the right target interval's end.
     */
    public static SVInterval getOuterIntrachromosomalLinkInterval(final EvidenceTargetLink link) {
        Utils.nonNull(link, "EvidenceTargetLink cannot be null");
        final PairedStrandedIntervals intervals = link.getPairedStrandedIntervals();
        final SVInterval leftInterval = intervals.getLeft().getInterval();
        final SVInterval rightInterval = intervals.getRight().getInterval();
        if (leftInterval.getContig() != rightInterval.getContig()) {
            throw new IllegalArgumentException("Link target intervals must be on the same contig");
        }
        return new SVInterval(leftInterval.getContig(), leftInterval.getStart(), rightInterval.getEnd());
    }

    /**
     * Returns new interval with the specified amount of padding
     */
    public static SVInterval getPaddedInterval(final SVInterval interval, final int padding, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(interval, "SVInterval cannot be null");
        if (padding < 0) {
            throw new IllegalArgumentException("Padding must be non-negative");
        }
        Utils.nonNull(dictionary, "SAMSequenceDictionary cannot be null");
        final SAMSequenceRecord sequence = dictionary.getSequence(interval.getContig());
        if (sequence == null) {
            throw new IllegalArgumentException("Could not find interval contig with index " + interval.getContig() + " in the sequence dictionary");
        }
        final int leftBoundary = Math.max(0, interval.getStart() - padding);
        final int contigLength = dictionary.getSequence(interval.getContig()).getSequenceLength();
        final int rightBoundary = Math.min(contigLength - 1, interval.getEnd() + padding);
        return new SVInterval(interval.getContig(), leftBoundary, rightBoundary);
    }

    /**
     * Gets all EvidenceTargetLink objects in the tree that overlap the given interval
     */
    public static Collection<EvidenceTargetLink> getOverlappingLinksOnInterval(final SVInterval interval, final SVIntervalTree<EvidenceTargetLink> tree) {
        Utils.nonNull(interval, "SVInterval cannot be null");
        return getTreeOverlapperStream(interval, tree).filter(node -> node.getInterval().getContig() == interval.getContig()).map(SVIntervalTree.Entry::getValue).collect(Collectors.toList());
    }

    /**
     * Returns true if all of interval b lies inside of interval a (inclusive)
     */
    public static boolean containsInterval(final SVInterval a, final SVInterval b) {
        Utils.nonNull(a, "SVInterval A cannot be null");
        Utils.nonNull(b, "SVInterval B cannot be null");
        return a.getContig() == b.getContig() && a.getStart() <= b.getStart() && a.getEnd() >= b.getEnd();
    }

    /**
     * Returns stream of tree Entry objects that overlap the given interval
     */
    private static <T> Stream<SVIntervalTree.Entry<T>> getTreeOverlapperStream(final SVInterval interval, final SVIntervalTree<T> tree) {
        Utils.nonNull(interval, "SVInterval cannot be null");
        Utils.nonNull(tree, "SVIntervalTree cannot be null");
        return Utils.stream(tree.overlappers(interval));
    }

    /**
     * Returns true if the tree contains an interval that has minimum recriprocal overlap fraction with the given interval
     */
    public static <T> boolean hasReciprocalOverlapInTree(final SVInterval interval, final SVIntervalTree<T> tree, final double minFractionOverlap) {
        return getTreeOverlapperStream(interval, tree).anyMatch(node -> hasReciprocalOverlap(interval, node.getInterval(), minFractionOverlap));
    }

    /**
     * Returns true if the two intervals have at least the given amount of reciprocal overlap
     */
    public static boolean hasReciprocalOverlap(final SVInterval a, final SVInterval b, final double minFractionOverlap) {
        Utils.validateArg(minFractionOverlap >= 0, "Overlap fraction must be non-negative");
        Utils.validateArg(minFractionOverlap <= 1, "Overlap fraction must less than or equal to 1");
        return reciprocalOverlap(a, b) >= minFractionOverlap;
    }

    /**
     * Calculates fraction reciprocal overlap of the given intervals
     */
    public static double reciprocalOverlap(final SVInterval a, final SVInterval b) {
        Utils.nonNull(a, "SVInterval A cannot be null");
        Utils.nonNull(b, "SVInterval B cannot be null");
        final int overlapLen = a.overlapLen(b);
        final double fractionOverlapA = overlapLen / (double) a.getLength();
        final double fractionOverlapB = overlapLen / (double) b.getLength();
        return Math.min(fractionOverlapA, fractionOverlapB);
    }

    /**
     * Converts SimpleInterval to SVInterval
     */
    public static SVInterval convertInterval(final SimpleInterval interval, final SAMSequenceDictionary dictionary) {
        final int sequenceIndex = dictionary.getSequenceIndex(interval.getContig());
        if (sequenceIndex == -1) {
            throw new IllegalArgumentException("Could not find contig " + interval.getContig() + " in sequence ditionary");
        }
        return new SVInterval(sequenceIndex, interval.getStart(), interval.getEnd());
    }

    /**
     * Converts SVInterval to SimpleInterval
     */
    public static SimpleInterval convertInterval(final SVInterval interval, final SAMSequenceDictionary dictionary) {
        final SAMSequenceRecord sequenceRecord = dictionary.getSequence(interval.getContig());
        if (sequenceRecord == null) {
            throw new IllegalArgumentException("Could not find contig " + interval.getContig() + " in sequence ditionary");
        }
        return new SimpleInterval(sequenceRecord.getSequenceName(), interval.getStart(), interval.getEnd());
    }

}
