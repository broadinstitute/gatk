package org.broadinstitute.hellbender.tools.copynumber.utils.germlinetagging;

import com.google.common.collect.ImmutableSortedMap;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.PeekableIterator;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.geometry.euclidean.oned.Interval;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * This utility class performs a simple tagging of germline segments in a tumor segments file.
 */
public class SimpleGermlineTagger {

    public static final String CNLOH_IN_GERMLINE = "C";

    private SimpleGermlineTagger(){}

    /** Look for concordant endpoints for determining whether an event in the tumor regions are in the germline regions.
     *
     * This method will not modify the input segments.
     *
     * @param initialTumorSegments segment file from the tumor sample.  Cannot contain intervals that overlap.  Must have
     *                             calls that are using the annotation from the parameter {@code callAnnotation}.
     *                             No calls can be the empty string.
     *                              Never {@code null}.
     * @param initialNormalSegments segment file from the normal sample.  Must have calls that are using the annotation
     *                              from the parameter {@code callAnnotation}.  No calls can be the empty string.
     *                              Cannot contain intervals that overlap.  Never {@code null}.
     * @param callAnnotation the annotation name (e.g. "CALL") used in both the tumor and normal segment files.  Never {@code null}.
     * @param dictionary sequence dictionary used to generate the segment files.  Used for sorting.  Never {@code null}.
     * @param outputAnnotationName annotation to append to the tumor segments.  If you specify one that exists, it will overwrite.
     *                             Cannot be {@code null} nor empty string.
     * @param paddingInBp The amount of slop to give (in bp) for matching segment endpoints.  Must be >= 0
     * @param reciprocalThreshold  The reciprocal threshold between the normal and tumor segment that must match in order
     *                             to match (whether or not matching breakpoints are found).  Must be 0.0 to 1.0.
     *
     * @return copy of the tumor segments (sorted) with the additional germline tag annotation (the name of which is outputAnnotationName).  Never {@code null}
     */
    public static List<AnnotatedInterval> tagTumorSegmentsWithGermlineActivity(final List<AnnotatedInterval> initialTumorSegments,
                                                                               final List<AnnotatedInterval> initialNormalSegments,
                                                                               final String callAnnotation,
                                                                               final SAMSequenceDictionary dictionary,
                                                                               final String outputAnnotationName, final int paddingInBp,
                                                                               final double reciprocalThreshold) {
        validateBasicInputParameters(initialTumorSegments, initialNormalSegments, callAnnotation, dictionary, outputAnnotationName, paddingInBp, reciprocalThreshold);

        final List<AnnotatedInterval> tumorSegments = IntervalUtils.sortLocatablesBySequenceDictionary(initialTumorSegments, dictionary);
        final List<AnnotatedInterval> normalSegments = IntervalUtils.sortLocatablesBySequenceDictionary(initialNormalSegments, dictionary);

        final Map<AnnotatedInterval, String> tumorSegsToGermlineCallMap = createAnnotatedIntervalToGermlineCallMap(tumorSegments, normalSegments, callAnnotation, paddingInBp, reciprocalThreshold, dictionary);

        return createdUpdatedAnnotatedIntervals(tumorSegments, Collections.singletonList(tumorSegsToGermlineCallMap), outputAnnotationName);
    }

    /**
     *  See {@link SimpleGermlineTagger#tagTumorSegmentsWithGermlineActivity(List, List, String, SAMSequenceDictionary, String, int, double)}
     *
     *  This will also attempt to look for areas that appear to be Copy-Neutral Loss-of-Heterozygosity regions in the normal segments.
     *  In this determinination, only the MAF in the normal sample are relevant.
     *
     * @param initialTumorSegments See {@link SimpleGermlineTagger#tagTumorSegmentsWithGermlineActivity(List, List, String, SAMSequenceDictionary, String, int, double)}
     * @param initialNormalSegments See {@link SimpleGermlineTagger#tagTumorSegmentsWithGermlineActivity(List, List, String, SAMSequenceDictionary, String, int, double)}
     * @param callAnnotation See {@link SimpleGermlineTagger#tagTumorSegmentsWithGermlineActivity(List, List, String, SAMSequenceDictionary, String, int, double)}
     * @param dictionary See {@link SimpleGermlineTagger#tagTumorSegmentsWithGermlineActivity(List, List, String, SAMSequenceDictionary, String, int, double)}
     * @param outputAnnotationName See {@link SimpleGermlineTagger#tagTumorSegmentsWithGermlineActivity(List, List, String, SAMSequenceDictionary, String, int, double)}
     * @param paddingInBp See {@link SimpleGermlineTagger#tagTumorSegmentsWithGermlineActivity(List, List, String, SAMSequenceDictionary, String, int, double)}
     * @param reciprocalThreshold See {@link SimpleGermlineTagger#tagTumorSegmentsWithGermlineActivity(List, List, String, SAMSequenceDictionary, String, int, double)}
     * @param mafMaxThreshold Maximum minor allelic fraction to accept as un-balanced segment.  Must be within 0.0 to 0.5.  Recommended value of 0.47
     * @param mafHiAnnotation Name of the annotation in the given annotated interval that has the upper bound of the Minor allelic fraction estimate.  Cannot be {@code null} nor "".
     * @param mafLoAnnotation Name of the annotation in the given annotated interval that has the lower bound of the Minor allelic fraction estimate.  Cannot be {@code null} nor "".
     * @return See {@link SimpleGermlineTagger#tagTumorSegmentsWithGermlineActivity(List, List, String, SAMSequenceDictionary, String, int, double)}
     */
    public static List<AnnotatedInterval> tagTumorSegmentsWithGermlineActivity(final List<AnnotatedInterval> initialTumorSegments,
                                                                               final List<AnnotatedInterval> initialNormalSegments,
                                                                               final String callAnnotation,
                                                                               final SAMSequenceDictionary dictionary,
                                                                               final String outputAnnotationName, final int paddingInBp,
                                                                               final double reciprocalThreshold, final double mafMaxThreshold,
                                                                               final int cnLoHCheckMaxSize, final String mafHiAnnotation, final String mafLoAnnotation) {
        validateBasicInputParameters(initialTumorSegments, initialNormalSegments, callAnnotation, dictionary, outputAnnotationName, paddingInBp, reciprocalThreshold);

        Utils.validateArg(initialNormalSegments.stream().noneMatch(s -> StringUtils.isEmpty(s.getAnnotationValue(mafHiAnnotation))),
                "All normal segments must have a maf high annotation.  Annotation (column header) name must be: " + mafHiAnnotation);
        Utils.validateArg(initialTumorSegments.stream().noneMatch(s -> StringUtils.isEmpty(s.getAnnotationValue(mafHiAnnotation))),
                "All tumor segments must have a maf high annotation.  Annotation (column header) name must be: " + mafHiAnnotation);
        Utils.validateArg(initialNormalSegments.stream().noneMatch(s -> StringUtils.isEmpty(s.getAnnotationValue(mafLoAnnotation))),
                "All normal segments must have a maf low annotation.  Annotation (column header) name must be: " + mafLoAnnotation);
        Utils.validateArg(initialTumorSegments.stream().noneMatch(s -> StringUtils.isEmpty(s.getAnnotationValue(mafLoAnnotation))),
                "All tumor segments must have a maf low annotation.  Annotation (column header) name must be: " + mafLoAnnotation);
        ParamUtils.inRange(mafMaxThreshold, 0.0, 0.5, "Max MAF threshold must be between 0.0 and 0.5");

        final List<AnnotatedInterval> tumorSegments = IntervalUtils.sortLocatablesBySequenceDictionary(initialTumorSegments, dictionary);
        final List<AnnotatedInterval> normalSegments = IntervalUtils.sortLocatablesBySequenceDictionary(initialNormalSegments, dictionary);

        final Map<AnnotatedInterval, String> tumorSegsToGermlineCallMap = createAnnotatedIntervalToGermlineCallMap(tumorSegments, normalSegments, callAnnotation, paddingInBp, reciprocalThreshold, dictionary);

        final Map<AnnotatedInterval, String> tumorSegsToGermlineCNLoHTag = createTumorNormalCNLoHCallMap(tumorSegments,
                normalSegments, reciprocalThreshold, callAnnotation, dictionary, mafMaxThreshold, cnLoHCheckMaxSize,
                mafHiAnnotation, mafLoAnnotation);
        final List<Map<AnnotatedInterval, String>> tagMaps = Arrays.asList(tumorSegsToGermlineCallMap, tumorSegsToGermlineCNLoHTag);

        return createdUpdatedAnnotatedIntervals(tumorSegments, tagMaps, outputAnnotationName);
    }

    private static List<AnnotatedInterval> createdUpdatedAnnotatedIntervals(final List<AnnotatedInterval> tumorSegments, final List<Map<AnnotatedInterval, String>> tagMaps, final String outputAnnotationName) {
        final Map<AnnotatedInterval, String> mergedMap = tagMaps.stream()
                .map(Map::entrySet)
                .flatMap(Set::stream)
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (valOriginal, valNew) -> valOriginal));

        return tumorSegments.stream()
                .map(s -> createTumorTaggedSegment(s, outputAnnotationName, mergedMap.getOrDefault(s, CalledCopyRatioSegment.Call.NEUTRAL.getOutputString())))
                .collect(Collectors.toList());
    }

    private static void validateBasicInputParameters(final List<AnnotatedInterval> initialTumorSegments, final List<AnnotatedInterval> initialNormalSegments, final String callAnnotation, final SAMSequenceDictionary dictionary, final String outputAnnotationName, final int paddingInBp, final double reciprocalThreshold) {
        Utils.nonNull(dictionary);
        Utils.nonNull(initialTumorSegments);
        Utils.nonNull(initialNormalSegments);
        Utils.nonNull(callAnnotation);
        IntervalUtils.validateNoOverlappingIntervals(initialTumorSegments);
        IntervalUtils.validateNoOverlappingIntervals(initialNormalSegments);
        ParamUtils.isPositiveOrZero(paddingInBp, "padding must be greater than or equal to zero.");
        ParamUtils.inRange(reciprocalThreshold, 0.0, 1.0, "Reciprocal threshold must be between 0.0 and 1.0");
        Utils.validateArg(!StringUtils.isEmpty(outputAnnotationName), "Output annotation name cannot be empty.");

        Utils.validateArg(initialNormalSegments.stream().noneMatch(s -> StringUtils.isEmpty(s.getAnnotationValue(callAnnotation))),
                "All normal segments must have a call.  Call annotation (column header) name must be: " + callAnnotation);
        Utils.validateArg(initialTumorSegments.stream().noneMatch(s -> StringUtils.isEmpty(s.getAnnotationValue(callAnnotation))),
                "All tumor segments must have a call.  Call annotation (column header) must be: " + callAnnotation);
    }

    private static Map<AnnotatedInterval, String> createAnnotatedIntervalToGermlineCallMap(final List<AnnotatedInterval> tumorSegments, final List<AnnotatedInterval> normalSegments, final String callAnnotation, final int paddingInBp, final double reciprocalThreshold, final SAMSequenceDictionary dictionary) {

        final List<AnnotatedInterval> mergedNormalSegments = mergedRegionsByAnnotation(callAnnotation, normalSegments);
        return createTumorNormalCopyRatioCallMap(tumorSegments, mergedNormalSegments, paddingInBp, reciprocalThreshold, callAnnotation, dictionary);
    }

    private static Map<AnnotatedInterval, String> createTumorNormalCopyRatioCallMap(final List<AnnotatedInterval> tumorSegments, final List<AnnotatedInterval> mergedNormalSegments, final int paddingInBp, final double reciprocalThreshold, final String callAnnotation, final SAMSequenceDictionary dictionary) {
        // Grab the merged normal segments that do not have a neutral call.
        final List<AnnotatedInterval> nonZeroMergedNormalSegments = mergedNormalSegments.stream()
                .filter(s -> !StringUtils.isEmpty(s.getAnnotationValue(callAnnotation)))
                .filter(s -> !s.getAnnotationValue(callAnnotation).equals(CalledCopyRatioSegment.Call.NEUTRAL.getOutputString()))
                .collect(Collectors.toList());

        final Map<AnnotatedInterval, List<AnnotatedInterval>> nonZeroMergedNormalSegmentsToTumorSegments = IntervalUtils.createOverlapMap(nonZeroMergedNormalSegments, tumorSegments, dictionary);
        return createTumorSegmentsToGermlineTagMap(nonZeroMergedNormalSegmentsToTumorSegments, paddingInBp, callAnnotation, reciprocalThreshold);
    }

    //TODO: Can't be string
    private static Map<AnnotatedInterval, String> createTumorNormalCNLoHCallMap(final List<AnnotatedInterval> tumorSegments, final List<AnnotatedInterval> normalSegments, final double reciprocalThreshold, final String callAnnotation, final SAMSequenceDictionary dictionary,
                                                                                final double mafMaxThreshold, final int cnLoHCheckMaxSize, final String mafHiAnnotation, final String mafLowAnnotation) {
        // Grab the normal segments that have a neutral call.
        final List<AnnotatedInterval> copyNeutralNormalSegments = normalSegments.stream()
                .filter(s -> !StringUtils.isEmpty(s.getAnnotationValue(callAnnotation)))
                .filter(s -> s.getAnnotationValue(callAnnotation).equals(CalledCopyRatioSegment.Call.NEUTRAL.getOutputString()))
                .collect(Collectors.toList());

        final Map<AnnotatedInterval, List<AnnotatedInterval>> copyNeutralNormalSegmentsToTumorSegments = IntervalUtils.createOverlapMap(copyNeutralNormalSegments, tumorSegments, dictionary);
        return createTumorSegmentsToGermlineCnLohMap(copyNeutralNormalSegmentsToTumorSegments, reciprocalThreshold, mafMaxThreshold, cnLoHCheckMaxSize, mafHiAnnotation, mafLowAnnotation);
    }


    private static Map<AnnotatedInterval, String> createTumorSegmentsToGermlineTagMap(final Map<AnnotatedInterval, List<AnnotatedInterval>> nonZeroMergedNormalSegmentsToTumorSegments, int paddingInBp, final String callAnnotation, final double reciprocalThreshold) {
        final Map<AnnotatedInterval, String> result = new HashMap<>();
        for (final AnnotatedInterval normalSeg : nonZeroMergedNormalSegmentsToTumorSegments.keySet()) {
            final List<AnnotatedInterval> overlappingTumorSegments = nonZeroMergedNormalSegmentsToTumorSegments.get(normalSeg);

            final boolean isSegmentPositionMatch = isSegmentPositionMatch(normalSeg, overlappingTumorSegments, paddingInBp, reciprocalThreshold, callAnnotation);

            if (isSegmentPositionMatch) {
                final CalledCopyRatioSegment.Call normalCall = Arrays.stream(CalledCopyRatioSegment.Call.values())
                        .filter(c -> c.getOutputString().equals(normalSeg.getAnnotationValue(callAnnotation))).findFirst().orElse(null);
                if (normalCall == null) {
                    throw new UserException.BadInput("No call exists in normal segment.  Does normal input have a call field \"" + callAnnotation + "\"?");
                }
                result.putAll(overlappingTumorSegments.stream()
                        .filter(s -> ((Math.abs(s.getStart() - normalSeg.getStart()) <= paddingInBp) || (Math.abs(normalSeg.getEnd() - s.getEnd()) <= paddingInBp)
                                || ((normalSeg.getStart() < s.getStart()) && (normalSeg.getEnd() > s.getEnd())))
                                && (normalSeg.getInterval().intersect(s).size() > (s.getInterval().size() * reciprocalThreshold))
                        )
                        .collect(Collectors.toMap(Function.identity(), s -> normalCall.getOutputString())));
            }
        }

        return result;
    }

    private static boolean isSegmentPositionMatch(final AnnotatedInterval normalSeg, final List<AnnotatedInterval> overlappingTumorSegments, final int paddingInBp, final double reciprocalThreshold, final String annotationToMergeTumorSegments) {
        final List<AnnotatedInterval> mergedTumorSegments = mergedRegionsByAnnotation(annotationToMergeTumorSegments, overlappingTumorSegments);
        final boolean isReciprocalOverlapSeen = mergedTumorSegments.stream()
            .anyMatch(s -> IntervalUtils.isReciprocalOverlap(s.getInterval(), normalSeg.getInterval(), reciprocalThreshold));

        final boolean isStartPositionSeen = overlappingTumorSegments.stream()
                .anyMatch(s -> Math.abs(s.getStart() - normalSeg.getStart()) <= paddingInBp);

        final boolean isEndPositionSeen = overlappingTumorSegments.stream()
                .anyMatch(s -> Math.abs(s.getEnd() - normalSeg.getEnd()) <= paddingInBp);

        return (isStartPositionSeen && isEndPositionSeen) || isReciprocalOverlapSeen;
    }

    //TODO: This cannot return a map to a string.  Must be an enum
    private static Map<AnnotatedInterval, String> createTumorSegmentsToGermlineCnLohMap(final Map<AnnotatedInterval, List<AnnotatedInterval>> copyNeutralNormalSegmentsToTumorSegments, final double reciprocalThreshold, final double mafMaxThreshold, final int cnLoHCheckMaxSize, final String mafHiAnnotation, final String mafLowAnnotation) {
        final Map<AnnotatedInterval, String> result = new HashMap<>();
        for (final AnnotatedInterval normalSeg : copyNeutralNormalSegmentsToTumorSegments.keySet()) {
            final boolean isNormalSegmentShort = normalSeg.getInterval().size() <= cnLoHCheckMaxSize;
            if (isNormalSegmentShort) {
                final List<AnnotatedInterval> overlappingTumorSegments = copyNeutralNormalSegmentsToTumorSegments.get(normalSeg);
                for (final AnnotatedInterval overlappingTumorSegment : overlappingTumorSegments) {

                    final boolean isReciprocalOverlapRegion = IntervalUtils.isReciprocalOverlap(normalSeg.getInterval(), overlappingTumorSegment.getInterval(), reciprocalThreshold);
                    final boolean isMafUnlikelyBalanced = Double.parseDouble(normalSeg.getAnnotationValue(mafHiAnnotation)) < mafMaxThreshold;

                    if (isReciprocalOverlapRegion && isMafUnlikelyBalanced) {
                        result.put(overlappingTumorSegment, CNLOH_IN_GERMLINE);
                    }
                }
            }
        }
        return result;
    }

    private static AnnotatedInterval createTumorTaggedSegment(final AnnotatedInterval tumorSegment, final String updateAnnotationName, final String updateAnnotationValue) {
        final Pair<String, String> updatePair = Pair.of(updateAnnotationName, updateAnnotationValue);
        return createTumorTaggedSegment(tumorSegment, Collections.singletonList(updatePair));
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

    private static Interval intersect(final Interval interval1, final Interval interval2) {
        // Inf is lower bound, sup is upper bound
        return new Interval(Math.max(interval1.getInf(), interval2.getInf()), Math.min(interval1.getSup(), interval2.getSup()));
    }

    private static boolean overlap(final Interval interval1, final Interval interval2) {
        // Inf is lower bound, sup is upper bound
        return (interval1.getInf() <= interval2.getSup()) && (interval2.getInf() <= interval1.getSup());
    }

    private static boolean isReciprocalOverlap(final Interval interval1, final Interval interval2, final double reciprocalOverlapThreshold) {
        return overlap(interval1, interval2) &&
                (intersect(interval1, interval2).getSize() >= (interval2.getSize() * reciprocalOverlapThreshold)) &&
                (intersect(interval1, interval2).getSize() >= (interval1.getSize() * reciprocalOverlapThreshold));
    }

    private static boolean isReciprocalOverlap(final AnnotatedInterval annotatedInterval1, final AnnotatedInterval annotatedInterval2, final double reciprocalOverlapThreshold,
                                               final String mafAnnotationLo, final String mafAnnotationHi) {
        final Interval interval1 = new Interval(Double.parseDouble(annotatedInterval1.getAnnotationValue(mafAnnotationLo)),
            Double.parseDouble(annotatedInterval1.getAnnotationValue(mafAnnotationHi)));
        final Interval interval2 = new Interval(Double.parseDouble(annotatedInterval2.getAnnotationValue(mafAnnotationLo)),
            Double.parseDouble(annotatedInterval2.getAnnotationValue(mafAnnotationHi)));

        return isReciprocalOverlap(interval1, interval2, reciprocalOverlapThreshold);
    }
}
