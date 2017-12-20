package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
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
        return Math.max(0, position - padding);
    }

    public static int getRightPaddedIntervalBoundary(final int position, final int padding, final int contigLength) {
        return Math.min(contigLength - 1, position + padding);
    }

    public static SVInterval getLargeLinkEventInterval(final EvidenceTargetLink link) {
        final PairedStrandedIntervals intervals = link.getPairedStrandedIntervals();
        final SVInterval leftInterval = intervals.getLeft().getInterval();
        final SVInterval rightInterval = intervals.getRight().getInterval();
        return new SVInterval(leftInterval.getContig(), leftInterval.getStart(), rightInterval.getEnd());
    }

    public static SVInterval getSmallLinkEventInterval(final EvidenceTargetLink link) {
        final PairedStrandedIntervals intervals = link.getPairedStrandedIntervals();
        final SVInterval leftInterval = intervals.getLeft().getInterval();
        final SVInterval rightInterval = intervals.getRight().getInterval();
        return new SVInterval(leftInterval.getContig(), leftInterval.getEnd(), rightInterval.getStart());
    }

    public static SVInterval getPaddedInterval(final SVInterval interval, final int padding, final SAMSequenceDictionary dict) {
        final int contigLength = dict.getSequence(interval.getContig()).getSequenceLength();
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

    public static <T> Stream<T> getTreeNodesWithReciprocalOverlap(final SVInterval interval, final SVIntervalTree<T> tree, final double minFractionOverlap) {
        return getReciprocalOverlapsInTree(interval, tree, minFractionOverlap).map(SVIntervalTree.Entry::getValue);
    }

    private static <T> Stream<SVIntervalTree.Entry<T>> getReciprocalOverlapsInTree(final SVInterval interval, final SVIntervalTree<T> tree, final double minFractionOverlap) {
        return getTreeOverlapperStream(interval, tree).filter(event -> hasReciprocalOverlap(interval, event.getInterval(), minFractionOverlap));
    }

    public static <T> boolean hasReciprocalOverlapInTree(final SVInterval interval, final SVIntervalTree<T> tree, final double minFractionOverlap) {
        return getTreeOverlapperStream(interval, tree).anyMatch(node -> hasReciprocalOverlap(interval, node.getInterval(), minFractionOverlap));
    }

    public static boolean hasReciprocalOverlap(final SVInterval a, final SVInterval b, final double minFractionOverlap) {
        return reciprocalOverlap(a, b) >= minFractionOverlap;
    }

    private static double reciprocalOverlap(SVInterval a, SVInterval b) {
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

    public static boolean segmentCallMatchesInterval(final CalledCopyRatioSegment segment, final SVInterval interval, final CalledCopyRatioSegment.Call callType, final ReadMetadata readMetadata, final double minSegmentReciprocalOverlap) {
        final SVInterval segmentInterval = new SVInterval(readMetadata.getContigID(segment.getContig()), segment.getStart(), segment.getEnd());
        return segment.getCall() == callType && IntervalUtils.hasReciprocalOverlap(segmentInterval, interval, minSegmentReciprocalOverlap);
    }

    public static List<CopyRatio> getCopyRatiosOnInterval(final SVInterval interval, final OverlapDetector<CopyRatio> overlapDetector, final int binsToTrim, final int contig, final SAMSequenceDictionary dict) {
        final int start = interval.getStart();
        final int end = interval.getEnd();
        final String contigName = dict.getSequence(contig).getSequenceName();
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
