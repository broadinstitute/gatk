package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.collections.Lists;

import java.io.File;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        oneLineSummary = "(EXPERIMENTAL) Combine the breakpoints of two segment files and annotate the resulting intervals with chosen columns from each file.",
        summary = "Combine the breakpoints of two segment files while preserving annotations.\n" +
                "This tool will load all segments into RAM.\n"+
        "Expected interval columns are: " + SimpleAnnotatedGenomicRegion.CONTIG_HEADER + ", " +
        SimpleAnnotatedGenomicRegion.START_HEADER + ", " + SimpleAnnotatedGenomicRegion.END_HEADER,
        programGroup = CopyNumberProgramGroup.class)
@BetaFeature
public class CombineSegmentBreakpoints extends GATKTool {

    public static final String COLUMNS_OF_INTEREST_LONG_NAME = "columnsOfInterest";
    public static final String COLUMNS_OF_INTEREST_SHORT_NAME = "cols";

    public static final String LABELS_LONG_NAME = "labels";
    public static final String LABELS_SHORT_NAME = "lbls";

    @Argument(
            doc = "Input segment files -- must be specified twice, but order does not matter.",
            fullName = CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.SEGMENTS_FILE_SHORT_NAME,
            maxElements = 2,
            minElements = 2
    )
    private List<File> segmentFiles = new ArrayList<>();

    @Argument(
            doc = "Input segment file labels -- these will appear as suffixes in case of collisions.  The specification order must correspond to the input segment files.",
            fullName = LABELS_LONG_NAME,
            shortName = LABELS_SHORT_NAME,
            maxElements = 2,
            minElements = 2,
            optional = true
    )
    private List<String> segmentFileLabels = Lists.newArrayList("1", "2");

    @Argument(
            doc = "List of columns in either segment file that should be reported in the output file.  If the column header exists in both, it will have the appropriate label appended as a suffix.",
            fullName = COLUMNS_OF_INTEREST_LONG_NAME,
            shortName = COLUMNS_OF_INTEREST_SHORT_NAME,
            minElements = 1
    )
    private Set<String> columnsOfInterest = new HashSet<>();

    @Argument(
            doc = "Output TSV file with combined segment breakpoints",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private File outputFile;

    @Override
    public boolean requiresReference() { return true;}

    @Override
    public void traverse() {

        final List<SimpleAnnotatedGenomicRegion> segments1 = SimpleAnnotatedGenomicRegion.readAnnotatedRegions(segmentFiles.get(0), columnsOfInterest);
        final List<SimpleAnnotatedGenomicRegion> segments2 = SimpleAnnotatedGenomicRegion.readAnnotatedRegions(segmentFiles.get(1), columnsOfInterest);

        // Check to see if we should warn the user that one or more columns of interest were not seen in any input file.
        final Set<String> allSeenAnnotations = Sets.union(segments1.get(0).getAnnotations().keySet(), segments2.get(0).getAnnotations().keySet());
        final Set<String> unusedColumnsOfInterest = Sets.difference(columnsOfInterest, allSeenAnnotations);
        if (unusedColumnsOfInterest.size() > 0) {
            final List<String> missingColumns = new ArrayList<>(unusedColumnsOfInterest);
            throw new UserException.BadInput("Some columns of interest specified by the user were not seen in any input files: " + StringUtils.join(missingColumns, ", "));
        }

        // Create a map of input to output headers.  I.e. the annotation in the segment1 to the output that should be written in the final file.
        //  This assumes that the keys in each entry of the list is the same.
        // TODO: If we want to support more than two segment files, this is the only bit that requires thinking.  Once this is solved, then this logic can go into a utility class.
        final Set<String> intersectingAnnotations = Sets.intersection(segments1.get(0).getAnnotations().keySet(), segments2.get(0).getAnnotations().keySet());

        // Create the obvious mappings that are identity then tack on new annotations for conflicts.
        // These are mappings that take the header name from the segment files and map to an output header to avoid conflicts.
        final Map<String, String> input1ToOutputHeaderMap = segments1.get(0).getAnnotations().keySet().stream().filter(a -> !intersectingAnnotations.contains(a))
                .collect(Collectors.toMap(Function.identity(), Function.identity()));
        intersectingAnnotations.forEach(a -> input1ToOutputHeaderMap.put(a, a + "_" + segmentFileLabels.get(0)));

        final Map<String, String> input2ToOutputHeaderMap = segments2.get(0).getAnnotations().keySet().stream().filter(a -> !intersectingAnnotations.contains(a))
                .collect(Collectors.toMap(Function.identity(), Function.identity()));
        intersectingAnnotations.forEach(a -> input2ToOutputHeaderMap.put(a, a + "_" + segmentFileLabels.get(1)));

        final List<SimpleAnnotatedGenomicRegion> finalList = annotateCombinedIntervals(segments1, segments2,
                Arrays.asList(input1ToOutputHeaderMap, input2ToOutputHeaderMap), getBestAvailableSequenceDictionary());

        SimpleAnnotatedGenomicRegion.writeAnnotatedRegionsAsTsv(finalList, outputFile);
    }

    /**
     *  Create intervals with breakpoints of segments1 and segments2 as described in {@link IntervalUtils::combineAndSortBreakpoints} and annotate
     *   the new intervals with annotations of interest in segments1 and segments2.
     *
     * @param segments1 a list of simple annotated regions
     * @param segments2 a list of simple annotated regions
     * @param inputToOutputHeaderMaps a list of maps (of length 2) for segments1 and segments2 respectively.
     *                                Each maps the annotation name as it appears in the segment list to the output annotation name it should get.
     *                                This is required typically to avoid conflicts in the output annotation names.
     *
     * @return a list of simple annotated regions that is a combination of the breakpoints in segments2 and contains all the annotations specified
     *  in inputToOutputHeaderMaps.
     */
    private List<SimpleAnnotatedGenomicRegion> annotateCombinedIntervals(final List<SimpleAnnotatedGenomicRegion> segments1, final List<SimpleAnnotatedGenomicRegion> segments2,
                                                                         final List<Map<String, String>> inputToOutputHeaderMaps,
                                                                         final SAMSequenceDictionary dictionary) {

        final List<List<SimpleAnnotatedGenomicRegion>> segmentLists = Arrays.asList(segments1, segments2);

        final List<Locatable> combinedIntervals = IntervalUtils.combineAndSortBreakpoints(
                segmentLists.get(0).stream().map(SimpleAnnotatedGenomicRegion::getInterval)
                        .collect(Collectors.toList()),
                segmentLists.get(1).stream().map(SimpleAnnotatedGenomicRegion::getInterval)
                        .collect(Collectors.toList()), dictionary);

        // Create a list of maps where each entry corresponds to overlapping intervals of the combined intervals.
        final List<Map<Locatable, List<SimpleAnnotatedGenomicRegion>>> combinedIntervalsToSegmentsMaps = segmentLists.stream()
                .map(segs -> IntervalUtils.createOverlapMap(combinedIntervals, segs, dictionary)).collect(Collectors.toList());

        final List<SimpleAnnotatedGenomicRegion> result = new ArrayList<>();

        for (final Locatable interval: combinedIntervals) {
            final SortedMap<String, String> intervalAnnotationMap = new TreeMap<>();

            for (int i = 0; i < combinedIntervalsToSegmentsMaps.size(); i ++) {
                final Map<Locatable, List<SimpleAnnotatedGenomicRegion>> combinedIntervalsToSegmentsMap = combinedIntervalsToSegmentsMaps.get(i);
                final Map<String, String> inputToOutputHeaderMap =  inputToOutputHeaderMaps.get(i);
                final List<SimpleAnnotatedGenomicRegion> matchingSegments = combinedIntervalsToSegmentsMap.get(interval);
                inputToOutputHeaderMap.forEach((k,v) -> {
                                if (matchingSegments.size() > 1) {
                                    logger.warn(interval + " had more than one segment: " + matchingSegments + " only annotating with the first match.");
                                }
                                if (matchingSegments.size() >= 1) {
                                    intervalAnnotationMap.put(v, matchingSegments.get(0).getAnnotations().getOrDefault(k, ""));
                                } else {
                                    intervalAnnotationMap.put(v, "");
                                }
                    });
            }

            result.add(new SimpleAnnotatedGenomicRegion(new SimpleInterval(interval), intervalAnnotationMap));
        }
        return result;
    }
}
