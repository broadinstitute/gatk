package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Helper/utility class for merging segments in a {@link SegmentedModel}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SegmentMergeUtils {
    private enum MergeDirection {
        LEFT, RIGHT, NONE
    }

    /**
     * Returns a new segment specified by the outermost breakpoints of the two segments to be joined.
     * The sample name of the current segment is used.  Segments must be on the same chromosome.
     * @param segment1  first segment to be joined
     * @param segment2  second segment to be joined
     * @return          a new segment constructed by joining the two segments
     */
    public static SimpleInterval mergeSegments(final SimpleInterval segment1, final SimpleInterval segment2) {
        if (!segment1.getContig().equals(segment2.getContig())) {
            throw new RuntimeException(String.format("Cannot join segments " +
                            segment1.getContig() + ":%d-%d and " +
                            segment2.getContig() + ":%d-%d on different chromosomes.",
                    segment1.getStart(), segment1.getEnd(), segment2.getStart(), segment2.getEnd()));
        }
        final int start = Math.min(segment1.getStart(), segment2.getStart());
        final int end = Math.max(segment1.getEnd(), segment2.getEnd());

        return new SimpleInterval(segment1.getContig(), start, end);
    }

    /**
     * Returns a distance between two segments based on copy-ratio distance.
     * This distance is infinite if the segments are on different chromosomes.
     * @param segment1  first segment
     * @param segment2  second segment
     * @param targets   targets (with log_2 copy ratios) used to calculate median copy ratios
     * @return          median-copy-ratio distance
     *                  (Double.POSITIVE_INFINITY if segments are on different chromosomes)
     */
    //TODO change to copy-ratio--allele-fraction Manhattan distance once allele-fraction initialization implemented
    private static double distance(final SimpleInterval segment1, final SimpleInterval segment2,
                                   final TargetCollection<TargetCoverage> targets) {
        if (segment1.getContig().equals(segment2.getContig())) {
            final double [] copyRatiosOfSegment1 =
                    targets.targets(segment1).stream().mapToDouble(t -> Math.pow(2., t.getCoverage())).toArray();
            final double [] copyRatiosOfSegment2 =
                    targets.targets(segment2).stream().mapToDouble(t -> Math.pow(2., t.getCoverage())).toArray();
            final double medianCopyRatioOfSegment1 = new Median().evaluate(copyRatiosOfSegment1);
            final double medianCopyRatioOfSegment2 = new Median().evaluate(copyRatiosOfSegment2);
            final double copyRatioDistance =
                    Math.abs(medianCopyRatioOfSegment1 - medianCopyRatioOfSegment2);
            return copyRatioDistance;
        }
        return Double.POSITIVE_INFINITY;
    }

    /**
     * Returns a distance between two segments (specified by their indices in a list) based on copy-ratio distance.
     * This distance is infinite if the segments are on different chromosomes.
     * @param segments  list of segments
     * @param targets   targets (with log_2 copy ratios) used to calculate median copy ratios
     * @param index1    index of first segment
     * @param index2    index of second segment
     * @return          median-copy-ratio distance
     *                  (Double.POSITIVE_INFINITY if segments are on different chromosomes or
     *                  when checking distance to left (right) of first (last) segment)
     *
     */
    private static double distance(final List<SimpleInterval> segments, final TargetCollection<TargetCoverage> targets,
                                   final int index1, final int index2) {
        if (index1 < 0 || index1 >= segments.size() || index2 < 0 || index2 >= segments.size()) {
            return Double.POSITIVE_INFINITY;
        }
        final SimpleInterval segment1 = segments.get(index1);
        final SimpleInterval segment2 = segments.get(index2);
        return distance(segment1, segment2, targets);
    }

    /**
     * Given a segment specified by an index, returns the direction of the adjacent segment that is closer in
     * distance given by {@link #distance}.  If distance to both adjacent segments is infinite, returns
     * MergeDirection.NONE.
     * @param segments  list of segments
     * @param targets   targets (with log_2 copy ratios) used to calculate median copy ratios
     * @param index     index of the segment to consider
     * @return          direction of the adjacent segment that is closer in median copy ratio
     *                  (MergeDirection.NONE if distance between segments is infinite)
     */
    private static MergeDirection calculateMergeDirection(final List<SimpleInterval> segments,
                                                          final TargetCollection<TargetCoverage> targets,
                                                          final int index) {
        Utils.validIndex(index, segments.size());
        final double leftDistance = distance(segments, targets, index, index - 1);
        final double rightDistance = distance(segments, targets, index, index + 1);
        if (rightDistance <= leftDistance && rightDistance != Double.POSITIVE_INFINITY) {
            return MergeDirection.RIGHT;
        } else if (leftDistance != Double.POSITIVE_INFINITY) {
            return MergeDirection.LEFT;
        }
        return MergeDirection.NONE;
    }

    /**
     * Returns the number of segments that contain a number of targets below a given threshold.
     * @param segments              list of segments
     * @param targets               targets to be segmented
     * @param targetNumberThreshold number of targets below which a segment is considered small
     * @return                      number of segments containing a number of targets strictly less than the threshold
     */
    public static int countSmallSegments(final List<SimpleInterval> segments,
                                         final TargetCollection <TargetCoverage> targets,
                                         final int targetNumberThreshold) {
        Utils.nonNull(segments, "The list of segments cannot be null.");
        Utils.nonNull(targets, "The list of targets cannot be null.");
        return (int) segments.stream()
                .filter(s -> targets.targetCount(s) < targetNumberThreshold)
                .count();
    }

    /**
     * Given a list of segments, returns a new, modifiable list of segments with the small segments (i.e., those
     * containing less than a specified number of targets) dropped.  Does not modify the input collections.
     * @param segments              list of segments (will not be modified)
     * @param targets               targets to be segmented
     * @param targetNumberThreshold number of targets below which a segment is considered small
     * @return                      new list of segments with small segments dropped, never {@code null}
     */
    public static List<SimpleInterval> dropSmallSegments(final List<SimpleInterval> segments,
                                                         final TargetCollection <TargetCoverage> targets,
                                                         final int targetNumberThreshold) {
        Utils.nonNull(segments, "The list of segments cannot be null.");
        Utils.nonNull(targets, "The list of targets cannot be null.");
        return segments.stream().filter(s -> targets.targetCount(s) >= targetNumberThreshold)
                .collect(Collectors.toList());
    }

    /**
     * Returns a new, modifiable  list of segments with small segments (i.e., those containing less than a specified
     * number of targets) merged.  The list of segments is traversed from beginning to end.  Upon arrival at a small
     * segment, the segment is repeatedly merged with its closest neighboring segment until it is above threshold,
     * then the traversal resumes.  Does not modify the input collections.
     * @param segments              original list of segments (will not be modified)
     * @param targets               targets (with log_2 copy ratios) to be segmented
     * @param targetNumberThreshold number of targets below which a segment is considered small
     * @return                      new list of segments with small segments merged, never {@code null}
     */
    public static List<SimpleInterval> mergeSmallSegments(final List<SimpleInterval> segments,
                                                          final TargetCollection <TargetCoverage> targets,
                                                          final int targetNumberThreshold) {
        Utils.nonNull(segments, "The list of segments cannot be null.");
        Utils.nonNull(targets, "The list of targets cannot be null.");
        final List<SimpleInterval> mergedSegments = new ArrayList<>(segments);
        int index = 0;
        while (index < mergedSegments.size()) {
            //if current segment is small, merge it with an adjacent segment
            if (targets.targetCount(mergedSegments.get(index)) < targetNumberThreshold) {
                final MergeDirection direction = calculateMergeDirection(mergedSegments, targets, index);
                if (direction == MergeDirection.LEFT) {
                    //current = merge(left, current), remove left, stay on current during next iteration
                    mergedSegments.set(index, mergeSegments(mergedSegments.get(index - 1), mergedSegments.get(index)));
                    mergedSegments.remove(index - 1);
                    index -= 2;
                } else if (direction == MergeDirection.RIGHT) {
                    //current = merge(current, right), remove right, stay on current during next iteration
                    mergedSegments.set(index, mergeSegments(mergedSegments.get(index), mergedSegments.get(index + 1)));
                    mergedSegments.remove(index + 1);
                    index--;
                }
            }
            index++; //if no merge performed, go to next segment during next iteration
        }
        //contigs containing only a single small segment do not get merged; drop these segments
        return dropSmallSegments(mergedSegments, targets, targetNumberThreshold);
    }
}
