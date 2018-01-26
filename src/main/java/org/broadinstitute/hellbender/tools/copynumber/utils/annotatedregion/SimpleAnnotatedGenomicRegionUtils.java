package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Consumer;

public class SimpleAnnotatedGenomicRegionUtils {

    // Can't instantiate this class.
    private SimpleAnnotatedGenomicRegionUtils() {}

    /** Merge SimpleAnnotatedGenomicRegions by whether the regions overlap.  Sorting is performed as well, so input
     * ordering is lost.
     *
     * When two overlapping regions are merged, annotations are merged as follows:
     *  - The annotation appears in one region:  Resulting region will have the annotation with the value of that one region
     *  - The annotation appears in both regions:  Resulting region will have both values separated by the {@code annotationSeparator}
     *  parameter.
     *
     * Only merges overlaps, not abutters.
     *
     * @param initialRegions regions to merge.  Never {@code null}
     * @param dictionary sequence dictionary to use for sorting.  Never {@code null}
     * @param annotationSeparator separator to use in the case of annotation conflicts.  Never {@code null}
     * @param progressUpdater Consuming function that can callback with the latest region that has been processed.
     *                        Never {@code null}.  Use a no-op function if you do not wish to provide a progress update.
     * @return Segments will be sorted by the sequence dictionary
     */
    public static List<SimpleAnnotatedGenomicRegion> mergeRegions(final List<SimpleAnnotatedGenomicRegion> initialRegions,
                                                                  final SAMSequenceDictionary dictionary, final String annotationSeparator,
                                                                  final Consumer<Locatable> progressUpdater) {

        Utils.nonNull(initialRegions);
        Utils.nonNull(dictionary);
        Utils.nonNull(annotationSeparator);
        Utils.nonNull(progressUpdater);

        final List<SimpleAnnotatedGenomicRegion> segments = IntervalUtils.sortLocatablesBySequenceDictionary(initialRegions,
                dictionary);

        final List<SimpleAnnotatedGenomicRegion> finalSegments = new ArrayList<>();
        final PeekableIterator<SimpleAnnotatedGenomicRegion> segmentsIterator = new PeekableIterator<>(segments.iterator());
        while (segmentsIterator.hasNext()) {
            SimpleAnnotatedGenomicRegion currentRegion = segmentsIterator.next();
            while (segmentsIterator.peek() != null && IntervalUtils.overlaps(currentRegion, segmentsIterator.peek())) {
                final SimpleAnnotatedGenomicRegion toBeMerged = segmentsIterator.next();
                currentRegion = merge(currentRegion, toBeMerged, annotationSeparator);
            }
            progressUpdater.accept(currentRegion);
            finalSegments.add(currentRegion);
        }
        return finalSegments;
    }

    /** Merges two simple annotated genomic regions.
     * Throws exception if the two regions cannot be merged.  This is usually due to being on different contigs.
     * When annotations conflict, use the separator to put separate values.
     */
    private static SimpleAnnotatedGenomicRegion merge(final SimpleAnnotatedGenomicRegion region1, final SimpleAnnotatedGenomicRegion region2,
                                                      final String separator) {
        final SimpleInterval interval = IntervalUtils.mergeSegments(region1.getInterval(), region2.getInterval());

        final Set<Map.Entry<String, String>> allEntries = Sets.union(region1.getAnnotations().entrySet(),
                region2.getAnnotations().entrySet());

        // For each remaining entry, if the annotation name only exists in one region, then just pass it through.
        //     if it exists in both entries, then merge it using the separator.
        final BiFunction<String, String,String> conflictFunction = (s1, s2) -> renderConflict(s1, s2, separator);
        final SortedMap<String, String> annotations = new TreeMap<>();
        allEntries.forEach(e -> annotations.put(e.getKey(), mergeAnnotationValue(e.getKey(), region1, region2, conflictFunction)));

        return new SimpleAnnotatedGenomicRegion(interval, annotations);
    }

    /**
     *  Return a merged annotation value for the two regions and given annotation name.  Automatically solves conflicts.
     *
     * @param annotationName the annotation to determine
     * @param region1 first region to merge
     * @param region2 second region to merge
     * @param conflictFunction the function to run to solve conflicts
     * @return string with the new, merged value of the annotation.
     */
    private static String mergeAnnotationValue(final String annotationName, final SimpleAnnotatedGenomicRegion region1,
                                               final SimpleAnnotatedGenomicRegion region2, final BiFunction<String, String, String> conflictFunction) {
        final boolean doesRegion1ContainAnnotation = region1.hasAnnotation(annotationName);
        final boolean doesRegion2ContainAnnotation = region2.hasAnnotation(annotationName);

        if (doesRegion1ContainAnnotation && doesRegion2ContainAnnotation) {

            // Both regions contain an annotation and presumably these are of different values.
            return conflictFunction.apply(region1.getAnnotationValue(annotationName),
                    region2.getAnnotationValue(annotationName));
        } else if (doesRegion1ContainAnnotation) {
            return region1.getAnnotationValue(annotationName);
        } else if (doesRegion2ContainAnnotation) {
            return region2.getAnnotationValue(annotationName);
        }

        return null;
    }

    /**
     * Renders a single string from two strings that (presumably) represent the same annotation in different regions of the genome.
     * This will simply separate the two given values with the separator.
     *
     * @param s1 the value from region 1
     * @param s2 the value from region 2
     * @param separator string to use for value separation
     * @return new string
     */
    private static String renderConflict(final String s1, final String s2, final String separator) {
        final String[] s1Vals = StringUtils.splitByWholeSeparator(s1, separator);
        final String[] s2Vals = StringUtils.splitByWholeSeparator(s2, separator);

        final Set<String> allValsSet = new HashSet<>(Arrays.asList(s1Vals));
        allValsSet.addAll(Arrays.asList(s2Vals));

        final List<String> allVals = new ArrayList<>(allValsSet);
        allVals.sort(String::compareTo);

        return Utils.join(separator, allVals);
    }
}
