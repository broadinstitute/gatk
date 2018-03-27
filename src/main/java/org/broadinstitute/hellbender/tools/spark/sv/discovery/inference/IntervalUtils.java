package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervals;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class IntervalUtils {


    public static int getLeftPaddedIntervalBoundary(final int position, final int padding) {
        if (padding < 0) {
            throw new IllegalArgumentException("Padding must be non-negative");
        }
        return Math.max(0, position - padding);
    }

    public static int getRightPaddedIntervalBoundary(final int position, final int padding, final int contigLength) {
        if (padding < 0) {
            throw new IllegalArgumentException("Padding must be non-negative");
        }
        return Math.min(contigLength - 1, position + padding);
    }

    public static SVInterval getOuterLinkInterval(final EvidenceTargetLink link) {
        final PairedStrandedIntervals intervals = link.getPairedStrandedIntervals();
        final SVInterval leftInterval = intervals.getLeft().getInterval();
        final SVInterval rightInterval = intervals.getRight().getInterval();
        return new SVInterval(leftInterval.getContig(), leftInterval.getStart(), rightInterval.getEnd());
    }

    public static SVInterval getPaddedInterval(final SVInterval interval, final int padding, final SAMSequenceDictionary dictionary) {
        final int contigLength = dictionary.getSequence(interval.getContig()).getSequenceLength();
        final int leftBoundary = getLeftPaddedIntervalBoundary(interval.getStart(), padding);
        final int rightBoundary = getRightPaddedIntervalBoundary(interval.getEnd(), padding, contigLength);
        return new SVInterval(interval.getContig(), leftBoundary, rightBoundary);
    }

    public static boolean linkEndsOverlapIntervals(final SVInterval leftInterval, final SVInterval rightInterval, final EvidenceTargetLink link) {
        return leftInterval.overlaps(link.getPairedStrandedIntervals().getLeft().getInterval()) && rightInterval.overlaps(link.getPairedStrandedIntervals().getRight().getInterval());
    }

    public static List<EvidenceTargetLink> getOverlappingLinksOnInterval(final SVInterval interval, final SVIntervalTree<EvidenceTargetLink> tree) {
        return getTreeOverlapperStream(interval, tree).filter(node -> node.getInterval().getContig() == interval.getContig()).map(SVIntervalTree.Entry::getValue).collect(Collectors.toList());
    }

    public static boolean containsInterval(final SVInterval a, final SVInterval b) {
        return a.getContig() == b.getContig() && a.getStart() <= b.getStart() && a.getEnd() >= b.getEnd();
    }

    public static <T> Stream<SVIntervalTree.Entry<T>> getTreeOverlapperStream(final SVInterval interval, final SVIntervalTree<T> tree) {
        return Utils.stream(tree.overlappers(interval));
    }

    public static <T> boolean hasReciprocalOverlapInTree(final SVInterval interval, final SVIntervalTree<T> tree, final double minFractionOverlap) {
        return getTreeOverlapperStream(interval, tree).anyMatch(node -> hasReciprocalOverlap(interval, node.getInterval(), minFractionOverlap));
    }

    public static boolean hasReciprocalOverlap(final SVInterval a, final SVInterval b, final double minFractionOverlap) {
        return reciprocalOverlap(a, b) >= minFractionOverlap;
    }

    public static double reciprocalOverlap(SVInterval a, SVInterval b) {
        if (a.getContig() != b.getContig()) return 0;
        if (a.getStart() > b.getStart()) {
            final SVInterval swap = a;
            a = b;
            b = swap;
        }
        if (a.getEnd() < b.getStart()) return 0;
        if (a.getEnd() > b.getEnd()) return b.getLength() / (double) a.getLength();
        final double overlap = a.getEnd() - b.getStart();
        return Math.min(overlap / a.getLength(), overlap / b.getLength());
    }

    public static int overlappingBases(final SimpleInterval a, final SVInterval b, final SAMSequenceDictionary dictionary) {
        return overlappingBases(convertInterval(a, dictionary), b);
    }

    public static double percentOverlap(final SimpleInterval a, final SVInterval b, final SAMSequenceDictionary dictionary) {
        return overlappingBases(a, b, dictionary) / (double) a.getLengthOnReference();
    }

    public static SVInterval convertInterval(final SimpleInterval interval, final SAMSequenceDictionary dictionary) {
        return new SVInterval(dictionary.getSequenceIndex(interval.getContig()), interval.getStart(), interval.getEnd());
    }

    public static int overlappingBases(final SVInterval a, final SVInterval b) {
        if (a.getContig() != b.getContig()) return 0;
        final int overlapStart = Math.max(a.getStart(), b.getStart());
        final int overlapEnd = Math.min(a.getEnd(), b.getEnd());
        return Math.max(0, overlapEnd - overlapStart);
    }

    public static List<CopyRatio> getCopyRatiosOnInterval(final SVInterval interval, final OverlapDetector<CopyRatio> overlapDetector, final int binsToTrim, final int contig, final SAMSequenceDictionary dictionary) {
        final int start = interval.getStart();
        final int end = interval.getEnd();
        final String contigName = dictionary.getSequence(contig).getSequenceName();
        final SimpleInterval simpleInterval = new SimpleInterval(contigName, start + 1, end + 1);
        if (simpleInterval.size() == 0) {
            return Collections.emptyList();
        }
        final List<CopyRatio> copyRatios = overlapDetector.getOverlaps(simpleInterval).stream().collect(Collectors.toList());

        if (copyRatios.size() <= 2 * binsToTrim) return Collections.emptyList();
        Collections.sort(copyRatios, Comparator.comparing(CopyRatio::getStart));
        return copyRatios.subList(binsToTrim, copyRatios.size() - binsToTrim);
    }

}
