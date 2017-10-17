package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.VersatileAnnotatedRegionParser;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

@CommandLineProgramProperties(
        oneLineSummary = "(EXPERIMENTAL) Do a breakpoint union of two segment files and annotate with chosen columns from each file.",
        summary = "Breakpoint union of two segment files while preserving annotations.\n" +
                "This tool will load all segments into RAM.",
        programGroup = CopyNumberProgramGroup.class)
@BetaFeature
public class UnionSegments extends GATKTool {

    static final String COLUMNS_OF_INTEREST_LONG_NAME = "columnsOfInterest";
    static final String COLUMNS_OF_INTEREST_SHORT_NAME = "cols";
    static final String LABELS_LONG_NAME = "labels";
    static final String LABELS_SHORT_NAME = "lbls";
    @Argument(
            doc = "Input segment files -- must be specified twice, but order does not matter.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            maxElements = 2, minElements = 2
    )
    protected List<File> segmentsFile;
    @Argument(
            doc = "Input segment file labels -- these will appear as suffixes in case of collisions.  The specification order must correspond to the input segment files.",
            fullName = LABELS_LONG_NAME,
            shortName = LABELS_SHORT_NAME,
            maxElements = 2, minElements = 2, optional = true
    )
    protected List<String> segmentsFileLabels = new ArrayList<>();

    @Argument(
            doc="List of columns in either segment file that should be reported in the output file.  If the column header exists in both, it will have an append.",
            fullName = COLUMNS_OF_INTEREST_LONG_NAME, shortName = COLUMNS_OF_INTEREST_SHORT_NAME, minElements = 1
    )
    protected Set<String> columnsOfInterest = new HashSet<>();

    @Argument(
            doc="Output tsv file with union'ed segments",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    protected File outputFile;

    @Override
    public void traverse() {
        try {

            final VersatileAnnotatedRegionParser parser = new VersatileAnnotatedRegionParser();
            final List<SimpleAnnotatedGenomicRegion> segments1 = parser.readAnnotatedRegions(segmentsFile.get(0), columnsOfInterest);
            final List<SimpleAnnotatedGenomicRegion> segments2 = parser.readAnnotatedRegions(segmentsFile.get(1), columnsOfInterest);
            final List<List<SimpleAnnotatedGenomicRegion>> inputSegmentsLists = Lists.newArrayList(segments1, segments2);

            if (segmentsFileLabels.size() == 0) {
                segmentsFileLabels = IntStream.range(0, inputSegmentsLists.size()).mapToObj(i -> String.valueOf(i+1)).collect(Collectors.toList());
            }

            // Check to see if we should warn the user that one or more columns of interest were not seen in any input file.
            final Set<String> allSeenAnnotations = Sets.union(segments1.get(0).getAnnotations().keySet(), segments2.get(0).getAnnotations().keySet());
            final Set<String> unusedColumnsOfInterest = Sets.difference(columnsOfInterest, allSeenAnnotations);
            if (unusedColumnsOfInterest.size() > 0) {
                final List<String> missingColumns = new ArrayList<>(unusedColumnsOfInterest);
                logger.warn("Some columns of interest specified by the user were not seen in any input files: " + StringUtils.join(missingColumns, ", "));
            }

            // Create a map of input to output headers.  I.e. the annotation in the segment1 to the output that should be written in the final file.
            //  This assumes that the keys in each entry of the list is the same.
            // TODO: If we want to support more than two segment files, this is the only bit that requires thinking.  Once this is solved, then this logic can go into a utility class.
            final Set<String> intersectingAnnotations = Sets.intersection(segments1.get(0).getAnnotations().keySet(), segments2.get(0).getAnnotations().keySet());

            // Create the obvious mappings that are identity then tack on new annotations for conflicts.
            // These are mappings that take the header name from the segment files and map to an output header to avoid conflicts.
            final Map<String, String> input1ToOutputHeaderMap = Utils.stream(segments1.get(0).getAnnotations().keySet().iterator()).filter(a -> !intersectingAnnotations.contains(a))
                    .collect(Collectors.toMap(Function.identity(), Function.identity()));
            Utils.stream(intersectingAnnotations.iterator()).forEach(a -> input1ToOutputHeaderMap.put(a, a + "_" + segmentsFileLabels.get(0)));

            final Map<String, String> input2ToOutputHeaderMap = Utils.stream(segments2.get(0).getAnnotations().keySet().iterator()).filter(a -> !intersectingAnnotations.contains(a))
                    .collect(Collectors.toMap(Function.identity(), Function.identity()));
            Utils.stream(intersectingAnnotations.iterator()).forEach(a -> input2ToOutputHeaderMap.put(a, a + "_" + segmentsFileLabels.get(1)));

            // Use the best available sequence dictionary, but if that is not available (i.e. null),
            //   create one from the segment files.

            final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary() != null ?
                    getBestAvailableSequenceDictionary() :
                    createSAMSequenceDictionaryFromIntervals(inputSegmentsLists.stream()
                            .flatMap(List::stream)
                            .collect(Collectors.toList()));

            final List<SimpleAnnotatedGenomicRegion> finalList = annotateUnionedIntervals(segments1, segments2,
                    sequenceDictionary, Lists.newArrayList(input1ToOutputHeaderMap, input2ToOutputHeaderMap));

            // TODO: Capture sample names in the comments if possible.
            // TODO:  Allow choice in output names of interval headers.
            final VersatileAnnotatedRegionParser writer = new VersatileAnnotatedRegionParser();
            writer.writeAnnotatedRegionsAsTsv(finalList, outputFile, Collections.emptyList(),
                    "contig", "start", "end");

        } catch (final IOException ioe) {
            throw new UserException.BadInput("Could not parse input file", ioe);
        }
    }

    /**
     *  Create a union of segments1 and segments2 as described in {@link IntervalUtils::unionIntervals} and annotate
     *   the union with annotations of interest in segments1 and segments2.
     *
     * @param segments1 a list of simple annotated regions
     * @param segments2 a list of simple annotated regions
     * @param dictionary a sequence dictionary
     * @param inputToOutputHeaderMaps a list of maps (of length 2) for segments1 and segments2 respectively.
     *                                Each maps the annotation name as it appears in the segment list to the output annotation name it should get.
     *                                This is required typically to avoid conflicts in the output annotation names.
     *
     * @return a list of simple annotated regions that is a union of segments1 and segments2 and contains all the annotations specified
     *  in inputToOutputHeaderMaps.
     */
    private List<SimpleAnnotatedGenomicRegion> annotateUnionedIntervals(final List<SimpleAnnotatedGenomicRegion> segments1, final List<SimpleAnnotatedGenomicRegion> segments2,
                                                                        final SAMSequenceDictionary dictionary, final List<Map<String, String>> inputToOutputHeaderMaps) {

        final List<List<SimpleAnnotatedGenomicRegion>> segmentLists = Lists.newArrayList(segments1, segments2);

        // We assume that the union'ed intervals are sorted.
        // TODO: This needs to be sorted if we want to support more than two sets of regions at once.
        final List<Locatable> unionIntervals = IntervalUtils.unionIntervals(
                segmentLists.get(0).stream().map(SimpleAnnotatedGenomicRegion::getInterval)
                        .collect(Collectors.toList()),
                segmentLists.get(1).stream().map(SimpleAnnotatedGenomicRegion::getInterval)
                        .collect(Collectors.toList()));

        // Create a list of maps where each entry corresponds to union'ed intervals to the regions in segmentList_i
        final List<Map<Locatable, List<SimpleAnnotatedGenomicRegion>>> unionIntervalsToSegmentsMaps = segmentLists.stream()
                .map(segs -> IntervalUtils.createOverlapMap(unionIntervals, segs, dictionary)).collect(Collectors.toList());

        final List<SimpleAnnotatedGenomicRegion> result = new ArrayList<>();

        for (final Locatable interval: unionIntervals) {
            final Map<String, String> intervalAnnotationMap = new HashMap<>();

            for (int i = 0; i < unionIntervalsToSegmentsMaps.size(); i ++) {
                final Map<Locatable, List<SimpleAnnotatedGenomicRegion>> unionIntervalsToSegmentsMap = unionIntervalsToSegmentsMaps.get(i);
                final Map<String, String> inputToOutputHeaderMap =  inputToOutputHeaderMaps.get(i);
                final List<SimpleAnnotatedGenomicRegion> matchingSegments = unionIntervalsToSegmentsMap.get(interval);
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

    private SAMSequenceDictionary createSAMSequenceDictionaryFromIntervals(final List<Locatable> locatables) {
        final Map<String, Integer> contigsToLengthMap = Utils.stream(locatables.stream().map(Locatable::getContig).collect(Collectors.toSet()).iterator())
                .collect(Collectors.toMap(Function.identity(), l -> 0));
        locatables.forEach(l -> {
            if (contigsToLengthMap.get(l.getContig()) < l.getEnd())
                contigsToLengthMap.put(l.getContig(), l.getEnd());
        });

        return new SAMSequenceDictionary(contigsToLengthMap.entrySet().stream()
                .map(e -> new SAMSequenceRecord(e.getKey(), e.getValue())).collect(Collectors.toList()));
    }
}
