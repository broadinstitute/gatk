package org.broadinstitute.hellbender.tools.copynumber.utils.germlinetagging;

import com.google.common.collect.ImmutableSortedMap;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.PeekableIterator;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * This utility class performs a simple tagging of germline segments in a tumor segments file.
 */
public class SimpleGermlineTagger {

    private SimpleGermlineTagger(){}

    /** Look for concordant endpoints for determining whether an event in the tumor regions are in the germline regions.
     *
     * This method will not modify the input tumor segments.
     *
     * @param initialTumorSegments segment file from the tumor sample.  Cannot contain intervals that overlap.
     *                              Never {@code null}.
     * @param initialNormalSegments segment file from the normal sample.  Must have calls that are using the annotation
     *                              from the parameter {@code callAnnotation}.  No calls can be the empty string.
     *                              Cannot contain intervals that overlap.  Never {@code null}.
     * @param callAnnotation the annotation name (e.g. "CALL") used in both the tumor and normal segment files.
     * @param dictionary sequence dictionary used to generate the segment files.  Used for sorting.  Never {@code null}.
     * @param outputAnnotationName annotation to append to the tumor segments.  If you specify one that exists, it will overwrite.
     *                             Cannot be {@code null} nor empty string.
     * @return copy of the tumor segments (sorted) with the additional germline tag annotation (the name of which is outputAnnotationName).
     */
    public static List<AnnotatedInterval> tagTumorSegmentsWithGermlineActivity(final List<AnnotatedInterval> initialTumorSegments,
                                                                               final List<AnnotatedInterval> initialNormalSegments,
                                                                               final String callAnnotation,
                                                                               final SAMSequenceDictionary dictionary,
                                                                               final String outputAnnotationName, final int paddingInBp) {
        Utils.nonNull(dictionary);
        Utils.nonNull(initialTumorSegments);
        Utils.nonNull(initialNormalSegments);
        IntervalUtils.validateNoOverlappingIntervals(initialTumorSegments);
        IntervalUtils.validateNoOverlappingIntervals(initialNormalSegments);

        final List<AnnotatedInterval> tumorSegments = IntervalUtils.sortLocatablesBySequenceDictionary(initialTumorSegments, dictionary);
        final List<AnnotatedInterval> normalSegments = IntervalUtils.sortLocatablesBySequenceDictionary(initialNormalSegments, dictionary);

        Utils.validateArg(!StringUtils.isEmpty(outputAnnotationName), "Output annotation name cannot be empty.");
        Utils.validateArg(normalSegments.stream().noneMatch(s -> StringUtils.isEmpty(s.getAnnotationValue(callAnnotation))),
                "All normal segments must have a call.  Column name must be: " + callAnnotation);

        final List<AnnotatedInterval> mergedNormalSegments = mergedRegionsByAnnotation(callAnnotation, normalSegments);

        // Grab the merged normal segments that do not have a neutral call.
        final List<AnnotatedInterval> nonZeroMergedNormalSegments = mergedNormalSegments.stream()
                .filter(s -> !StringUtils.isEmpty(s.getAnnotationValue(callAnnotation)))
                .filter(s -> !s.getAnnotationValue(callAnnotation).equals(CalledCopyRatioSegment.Call.NEUTRAL.getOutputString()))
                .collect(Collectors.toList());

        final Map<AnnotatedInterval, List<AnnotatedInterval>> nonZeroMergedNormalSegmentsToTumorSegments = IntervalUtils.createOverlapMap(nonZeroMergedNormalSegments, tumorSegments, dictionary);
        final Map<AnnotatedInterval, CalledCopyRatioSegment.Call> tumorSegsToGermlineTag = createTumorSegmentsToGermlineTagMap(nonZeroMergedNormalSegmentsToTumorSegments, paddingInBp, callAnnotation);

        return tumorSegments.stream()
                .map(s -> createTumorTaggedSegment(s, outputAnnotationName, tumorSegsToGermlineTag.getOrDefault(s, CalledCopyRatioSegment.Call.NEUTRAL).getOutputString()))
                .collect(Collectors.toList());
    }

    private static Map<AnnotatedInterval, CalledCopyRatioSegment.Call> createTumorSegmentsToGermlineTagMap(final Map<AnnotatedInterval, List<AnnotatedInterval>> nonZeroMergedNormalSegmentsToTumorSegments, int paddingInBp, final String callAnnotation) {
        final Map<AnnotatedInterval,CalledCopyRatioSegment.Call> result = new HashMap<>();
        for (final AnnotatedInterval normalSeg : nonZeroMergedNormalSegmentsToTumorSegments.keySet()) {
            final List<AnnotatedInterval> overlappingTumorSegments = nonZeroMergedNormalSegmentsToTumorSegments.get(normalSeg);
            final boolean isStartPositionSeen = overlappingTumorSegments.stream()
                    .anyMatch(s -> Math.abs(s.getStart() - normalSeg.getStart()) <= paddingInBp);

            final boolean isEndPositionSeen = overlappingTumorSegments.stream()
                    .anyMatch(s -> Math.abs(s.getEnd() - normalSeg.getEnd()) <= paddingInBp);
            if (isStartPositionSeen && isEndPositionSeen) {
                final CalledCopyRatioSegment.Call normalCall = Arrays.stream(CalledCopyRatioSegment.Call.values())
                        .filter(c -> c.getOutputString().equals(normalSeg.getAnnotationValue(callAnnotation))).findFirst().orElse(null);
                if (normalCall == null) {
                    throw new UserException.BadInput("No call exists in normal segment.  Does normal input have a call field?");
                }
                result.putAll(overlappingTumorSegments.stream()
                        .filter(s -> (Math.abs(s.getStart() - normalSeg.getStart()) <= paddingInBp) || (Math.abs(normalSeg.getEnd() - s.getEnd()) <= paddingInBp)
                                || ((normalSeg.getStart() < s.getStart()) && (normalSeg.getEnd() > s.getEnd())) )
                        .collect(Collectors.toMap(Function.identity(), s -> normalCall)));
            }
        }

        return result;
    }

    private static AnnotatedInterval createTumorTaggedSegment(final AnnotatedInterval tumorSegment, final String updateAnnotationName, final String updateAnnotationValue) {
        final Pair<String, String> updatePair = Pair.of(updateAnnotationName, updateAnnotationValue);
        return createTumorTaggedSegment(tumorSegment, Arrays.asList(updatePair));
    }

    private static AnnotatedInterval createTumorTaggedSegment(final AnnotatedInterval tumorSegment, final List<Pair<String,String>> updateAnnotationNameValuePairs) {

        // Copy the input segment
        final SortedMap<String, String> newAnnotations = new TreeMap<>();
        newAnnotations.putAll(tumorSegment.getAnnotations());
        updateAnnotationNameValuePairs.forEach(p -> newAnnotations.put(p.getKey(), p.getValue()));

        return new AnnotatedInterval(tumorSegment.getInterval(), newAnnotations);
    }

    private static List<AnnotatedInterval> mergedRegionsByAnnotation(final String annotationToMerge, final List<AnnotatedInterval> regions) {
        final List<AnnotatedInterval> mergedSegments = new ArrayList<>();
        final PeekableIterator<AnnotatedInterval> segmentsIterator = new PeekableIterator<>(regions.iterator());
        while (segmentsIterator.hasNext()) {
            AnnotatedInterval normalSegment = segmentsIterator.next();
            AnnotatedInterval nextSegment = segmentsIterator.peek();

            int updatedEndPoint =  normalSegment.getEnd();

            // Merge (if any to merge)
            while (segmentsIterator.hasNext() && isMergeableByAnnotation(annotationToMerge, normalSegment, nextSegment)) {
                updatedEndPoint = nextSegment.getEnd();
                segmentsIterator.next();
                nextSegment = segmentsIterator.peek();
            }
            final AnnotatedInterval segmentToAddToResult = new AnnotatedInterval(
                    new SimpleInterval(normalSegment.getContig(), normalSegment.getStart(), updatedEndPoint),
                    ImmutableSortedMap.of(annotationToMerge, normalSegment.getAnnotationValue(annotationToMerge)));
            mergedSegments.add(segmentToAddToResult);
        }
        return mergedSegments;
    }

    private static boolean isMergeableByAnnotation(final String annotationToMerge, final AnnotatedInterval segment, final AnnotatedInterval nextSegment) {
        return segment.getAnnotationValue(annotationToMerge).equals(nextSegment.getAnnotationValue(annotationToMerge)) &&
                segment.getContig().equals(nextSegment.getContig()) &&
                segment.getAnnotationValue(annotationToMerge) != null;
    }
}
