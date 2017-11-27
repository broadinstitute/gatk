package org.broadinstitute.hellbender.tools.copynumber.utils.germlinetagging;

import com.google.common.collect.ImmutableSortedMap;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.PeekableIterator;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.tools.copynumber.coverage.caller.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * This utility class performs a simple tagging of germline segments in a tumor segments file.
 */
public class SimpleGermlineTagger {

    private SimpleGermlineTagger(){}

    /**
     * @param initialTumorSegments segment file from the tumor sample.  Cannot contain intervals that overlap.
     *                              Never {@code null}.
     * @param initialNormalSegments segment file from the normal sample.  Must have calls that are using the annotation
     *                              from the parameter {@code callAnnotation}.  No calls can be the empty string.
     *                              Cannot contain intervals that overlap.  Never {@code null}.
     * @param callAnnotation the annotation name (e.g. "CALL") used in both the tumor and normal segment files.
     * @param dictionary sequence dictionary used to generate the segment files.  Used for sorting.  Never {@code null}.
     * @param outputAnnotationName annotation to append to the tumor segments.  If you specify one that exists, it will overwrite.
     *                             Cannot be {@code null} nor empty string.
     * @return the tumor segments (sorted) with the additional germline tag annotation (the name of which is outputAnnotationName).
     */
    public static List<SimpleAnnotatedGenomicRegion> tagTumorSegmentsWithGermlineActivity(final List<SimpleAnnotatedGenomicRegion> initialTumorSegments,
                                                                                          final List<SimpleAnnotatedGenomicRegion> initialNormalSegments,
                                                                                          final String callAnnotation,
                                                                                          final SAMSequenceDictionary dictionary,
                                                                                          final String outputAnnotationName, final int paddingInBp) {

        Utils.nonNull(dictionary);
        Utils.nonNull(initialTumorSegments);
        Utils.nonNull(initialNormalSegments);
        IntervalUtils.validateNoOverlappingIntervals(initialTumorSegments);
        IntervalUtils.validateNoOverlappingIntervals(initialNormalSegments);

        final List<SimpleAnnotatedGenomicRegion> tumorSegments = IntervalUtils.sortLocatablesBySequenceDictionary(initialTumorSegments, dictionary);
        final List<SimpleAnnotatedGenomicRegion> normalSegments = IntervalUtils.sortLocatablesBySequenceDictionary(initialNormalSegments, dictionary);

        Utils.validateArg(!StringUtils.isEmpty(outputAnnotationName), "Output annotation name cannot be empty.");
        Utils.validateArg(normalSegments.stream().noneMatch(s -> StringUtils.isEmpty(s.getAnnotationValue(callAnnotation))),
                "All normal segments must have a call.  Column name must be: " + callAnnotation);

        final List<SimpleAnnotatedGenomicRegion> mergedNormalSegments = determineMergedSegmentsByAnnotation(callAnnotation, normalSegments);

        // Grab the merged normal segments that do not have a neutral call.
        final List<SimpleAnnotatedGenomicRegion> nonZeroMergedNormalSegments = mergedNormalSegments.stream()
                .filter(s -> !StringUtils.isEmpty(s.getAnnotations().get(callAnnotation)))
                .filter(s -> !s.getAnnotations().get(callAnnotation).equals(CalledCopyRatioSegment.Call.NEUTRAL.getOutputString()))
                .collect(Collectors.toList());

        final OverlapDetector<SimpleAnnotatedGenomicRegion> overlapDetector = OverlapDetector.create(tumorSegments);

        // First initialize an annotation for all of the tumor segments that state that there is no germline influence
        tumorSegments.forEach(s -> s.setAnnotation(outputAnnotationName, CalledCopyRatioSegment.Call.NEUTRAL.getOutputString()));

        for (final SimpleAnnotatedGenomicRegion nonZeroMergedNormalSegment : nonZeroMergedNormalSegments) {

            // This code assumes that the overlap detector will provide references to the tumor segments (as opposed to copies)
            final Set<SimpleAnnotatedGenomicRegion> overlappingTumorSegments = overlapDetector.getOverlaps(nonZeroMergedNormalSegment);

            // We need to see that normal segment start and end position (with padding) represented in start and end positions of the tumor segments.
            final boolean isStartPositionSeen = overlappingTumorSegments.stream()
                    .anyMatch(s -> Math.abs(s.getStart() - nonZeroMergedNormalSegment.getStart()) < paddingInBp);

            final boolean isEndPositionSeen = overlappingTumorSegments.stream()
                    .anyMatch(s -> Math.abs(s.getEnd() - nonZeroMergedNormalSegment.getEnd()) < paddingInBp);

            // TODO: There are still minor bugs here.  Mostly if a segment is smaller than the padding.
            if (isEndPositionSeen && isStartPositionSeen) {
                overlappingTumorSegments.forEach(s -> s.setAnnotation(outputAnnotationName,
                        nonZeroMergedNormalSegment.getAnnotationValue(callAnnotation)));
            }
        }
        return tumorSegments;
    }

    private static List<SimpleAnnotatedGenomicRegion> determineMergedSegmentsByAnnotation(String annotationToMerge, List<SimpleAnnotatedGenomicRegion> normalSegments) {
        final List<SimpleAnnotatedGenomicRegion> mergedSegments = new ArrayList<>();
        final PeekableIterator<SimpleAnnotatedGenomicRegion> normalSegmentsIterator = new PeekableIterator<>(normalSegments.iterator());
        while (normalSegmentsIterator.hasNext()) {
            SimpleAnnotatedGenomicRegion normalSegment = normalSegmentsIterator.next();
            SimpleAnnotatedGenomicRegion nextNormalSegment = normalSegmentsIterator.peek();
            final SimpleAnnotatedGenomicRegion segmentToAddToResult = new SimpleAnnotatedGenomicRegion(normalSegment.getInterval(),
                    ImmutableSortedMap.of(annotationToMerge, normalSegment.getAnnotationValue(annotationToMerge)));

            // Merge (if any to merge)
            while (normalSegmentsIterator.hasNext() && isMergeableByAnnotation(annotationToMerge, normalSegment, nextNormalSegment)) {
                segmentToAddToResult.setEnd(nextNormalSegment.getEnd());
                normalSegmentsIterator.next();
                nextNormalSegment = normalSegmentsIterator.peek();
            }
            mergedSegments.add(segmentToAddToResult);
        }
        return mergedSegments;
    }

    private static boolean isMergeableByAnnotation(String annotationToMerge, SimpleAnnotatedGenomicRegion normalSegment, SimpleAnnotatedGenomicRegion nextNormalSegment) {
        return normalSegment.getAnnotationValue(annotationToMerge).equals(nextNormalSegment.getAnnotationValue(annotationToMerge)) &&
                normalSegment.getContig().equals(nextNormalSegment.getContig()) &&
                normalSegment.getAnnotationValue(annotationToMerge) != null;
    }

}
