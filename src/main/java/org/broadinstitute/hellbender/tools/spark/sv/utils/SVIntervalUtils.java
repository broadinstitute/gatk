package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth.SVCopyNumberInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.Comparator;

/**
 * Common utility functions for {@link SVInterval} and {@link SVIntervalTree}
 */
public class SVIntervalUtils {

    /**
     * Builds an interval tree of {@link SVInterval}s with null entry values
     */
    public static SVIntervalTree<Object> buildIntervalTreeWithNullValues(final Collection<SVInterval> intervals) {
        final SVIntervalTree<Object> tree = new SVIntervalTree<>();
        if (intervals != null) {
            for (final SVInterval interval : intervals) {
                tree.put(new SVInterval(interval.getContig(), interval.getStart(), interval.getEnd()), null);
            }
        }
        return tree;
    }

    /**
     * Gets comparator for {@link SVCopyNumberInterval} dictionary order
     */
    public static Comparator<SVCopyNumberInterval> getCopyNumberIntervalDictionaryOrderComparator() {
        return (o1, o2) -> compareIntervals(o1.getInterval(), o2.getInterval());
    }

    /**
     * Comparator function for {@link SVInterval} dictionary order
     */
    public static int compareIntervals(final SVInterval first, final SVInterval second) {
        Utils.nonNull(first);
        Utils.nonNull(second);

        int result = 0;
        if (first != second) {
            // compare the contigs
            result = Integer.compare(first.getContig(), second.getContig());
            if (result == 0) {
                // compare start position
                result = Integer.compare(first.getStart(), second.getStart());
                if (result == 0) {
                    // compare end position
                    result = Integer.compare(first.getEnd(), second.getEnd());
                }
            }
        }
        return result;
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
     * Returns true if the two intervals have at least the given amount of reciprocal overlap
     */
    public static boolean hasReciprocalOverlap(final SVInterval a, final SVInterval b, final double minOverlapFraction) {
        Utils.validateArg(minOverlapFraction >= 0, "Overlap fraction must be non-negative");
        Utils.validateArg(minOverlapFraction <= 1, "Overlap fraction must less than or equal to 1");
        return reciprocalOverlap(a, b) >= minOverlapFraction;
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
     * Creates new {@link SimpleInterval} from contig index, assuming 0-based half-open input coordinates
     */
    public static SimpleInterval createSimpleInterval(final int contig, final int start, final int end, final SAMSequenceDictionary dictionary) {
        final SAMSequenceRecord sequenceRecord = dictionary.getSequence(contig);
        if (sequenceRecord == null) {
            throw new IllegalArgumentException("Could not find contig " + contig + " in sequence ditionary");
        }
        return new SimpleInterval(sequenceRecord.getSequenceName(), start + 1, end);
    }

    /**
     * Converts {@link SVInterval} to {@link SimpleInterval}, assuming 0-based half-open input coordinates
     */
    public static SimpleInterval convertToSimpleInterval(final SVInterval interval, final SAMSequenceDictionary dictionary) {
        return createSimpleInterval(interval.getContig(), interval.getStart(), interval.getEnd(), dictionary);
    }

    /**
     * Generates an interval tree from {@link SVCopyNumberInterval} collection
     */
    public static SVIntervalTree<SVCopyNumberInterval> buildCopyNumberIntervalTree(final Collection<SVCopyNumberInterval> copyNumberIntervals) {
        final SVIntervalTree<SVCopyNumberInterval> tree = new SVIntervalTree<>();
        for (final SVCopyNumberInterval interval : copyNumberIntervals) {
            tree.put(interval.getInterval(), interval);
        }
        return tree;
    }
}
