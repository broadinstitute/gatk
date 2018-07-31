package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.*;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        oneLineSummary = "Combine the breakpoints of two segment files and annotate the resulting intervals with chosen columns from each file.",
        summary = "Combine the breakpoints of two segment files while preserving annotations.\n" +
                "This tool will load all segments into RAM.\n"+
                "Comments lines start with '#' or the data can be prepended by a SAM File Header (not both).\n" +
                "Output file comments will be the first segment file comments concatenated with the second file comments.\n" +
                "SAMFileHeaders in the input seg files are supported and will be merged in the output.\n" +
                "Output seg file will have a SAMFileHeader, even if the inputs did not.  If neither input file has " +
                  "a SAMFileHeader, then a reference must be specified.",
        programGroup = CopyNumberProgramGroup.class)
@ExperimentalFeature
public class CombineSegmentBreakpoints extends GATKTool {

    private static final Logger logger = LogManager.getLogger(CombineSegmentBreakpoints.class);

    public static final String COLUMNS_OF_INTEREST_LONG_NAME = "columns-of-interest";
    public static final String LABELS_LONG_NAME = "labels";

    @Argument(
            doc = "Input segment files -- must be specified twice, but order does not matter.",
            fullName = CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME,
            maxElements = 2,
            minElements = 2
    )
    private List<File> segmentFiles = new ArrayList<>();

    @Argument(
            doc = "Input segment file labels -- these will appear as suffixes in case of collisions.  The specification order must correspond to the input segment files.",
            fullName = LABELS_LONG_NAME,
            maxElements = 2,
            minElements = 2,
            optional = true
    )
    private List<String> segmentFileLabels = Lists.newArrayList("1", "2");

    @Argument(
            doc = "List of columns in either segment file that should be reported in the output file.  If the column header exists in both, it will have the appropriate label appended as a suffix.",
            fullName = COLUMNS_OF_INTEREST_LONG_NAME,
            minElements = 1
    )
    private Set<String> columnsOfInterest = new HashSet<>();

    @Argument(
            doc = "Output TSV file with combined segment breakpoints",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private File outputFile;

    @Override
    public void traverse() {

        final AnnotatedIntervalCollection annotatedIntervalCollection1 = AnnotatedIntervalCollection.create(segmentFiles.get(0).toPath(), columnsOfInterest);
        final List<AnnotatedInterval> segments1 = annotatedIntervalCollection1.getRecords();

        final AnnotatedIntervalCollection annotatedIntervalCollection2 = AnnotatedIntervalCollection.create(segmentFiles.get(1).toPath(), columnsOfInterest);
        final List<AnnotatedInterval> segments2 = annotatedIntervalCollection2.getRecords();

        final SAMFileHeader outputSamFileHeader = createOutputSamFileHeader(annotatedIntervalCollection1, annotatedIntervalCollection2);

        // Check to see if we should warn the user that one or more columns of interest were not seen in any input file.
        final Set<String> allSeenAnnotations = Sets.union(new HashSet<>(annotatedIntervalCollection1.getAnnotations()),
                new HashSet<>(annotatedIntervalCollection2.getAnnotations()));
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

        final List<AnnotatedInterval> finalList = annotateCombinedIntervals(segments1, segments2,
                Arrays.asList(input1ToOutputHeaderMap, input2ToOutputHeaderMap), outputSamFileHeader.getSequenceDictionary(),
                l -> progressMeter.update(l));

        final List<String> finalAnnotations = Lists.newArrayList(finalList.get(0).getAnnotations().keySet());
        finalAnnotations.sort(String::compareTo);

        final AnnotatedIntervalCollection finalCollection =
                AnnotatedIntervalCollection.create(finalList, outputSamFileHeader, finalAnnotations);
        finalCollection.write(outputFile);
    }

    /**
     * Do all the processing to create a SAMFileHeader for the given AnnotatedIntervalCollections.
     * Throws an exception if unable to generate a header with a non-empty sequence dictionary.
     *
     * @param annotatedIntervalCollection1 Never {@code null}
     * @param annotatedIntervalCollection2 Never {@code null}
     * @return Merged SAMFileHeader for the two annotated interval collections.  Never {@code null}.
     */
    private SAMFileHeader createOutputSamFileHeader(final AnnotatedIntervalCollection annotatedIntervalCollection1, final AnnotatedIntervalCollection annotatedIntervalCollection2) {
        Utils.nonNull(annotatedIntervalCollection1);
        Utils.nonNull(annotatedIntervalCollection2);

        final SAMFileHeader samFileHeader1 = annotatedIntervalCollection1.getSamFileHeader();
        final SAMFileHeader samFileHeader2 = annotatedIntervalCollection2.getSamFileHeader();

        assertSequenceDictionaryCompatibility(samFileHeader1, samFileHeader2);
        warnIfOnlyOneSequenceDictionarySpecified(samFileHeader1, samFileHeader2);

        // Since this tool will sort the segments by coordinate, this will be coordinate sorted, regardless of the input.
        final SamFileHeaderMerger samFileHeaderMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate,
                Arrays.asList(samFileHeader1, samFileHeader2),true);

        final SAMFileHeader outputSamFileHeader = samFileHeaderMerger.getMergedHeader();
        if ((outputSamFileHeader.getSequenceDictionary().getReferenceLength() == 0) && (getBestAvailableSequenceDictionary() == null)) {
            throw new UserException.BadInput("Cannot assemble a reference dictionary.  In order to use this tool, one " +
                    "of the following conditions must be satisfied:  1)  One or both input files have a SAM File header " +
                    "... 2)  A reference is provided (-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME + ")");
        }
        if (outputSamFileHeader.getSequenceDictionary().getReferenceLength() == 0) {
            outputSamFileHeader.setSequenceDictionary(getBestAvailableSequenceDictionary());
        }
        return outputSamFileHeader;
    }

    private void warnIfOnlyOneSequenceDictionarySpecified(SAMFileHeader samFileHeader1, SAMFileHeader samFileHeader2) {
        final boolean isSamFileHeader1Empty = samFileHeader1.getSequenceDictionary().getReferenceLength() == 0;
        final boolean isSamFileHeader2Empty = samFileHeader2.getSequenceDictionary().getReferenceLength() == 0;

        if (isSamFileHeader1Empty ^ isSamFileHeader2Empty) {
            logger.warn("One of the input files does not have a SAMFileHeader with a sequence dictionary.  Assuming that the specified SAMFileHeader dictionary is applicable to both input collections.");
        }
    }

    private static void assertSequenceDictionaryCompatibility(SAMFileHeader samFileHeader1, SAMFileHeader samFileHeader2) {
        final SequenceDictionaryUtils.SequenceDictionaryCompatibility compatibilityResult = SequenceDictionaryUtils.compareDictionaries(samFileHeader1.getSequenceDictionary(), samFileHeader2.getSequenceDictionary(), true);
        switch ( compatibilityResult ) {
            case UNEQUAL_COMMON_CONTIGS:
                throw new UserException.BadInput("Input files had common contigs with different lengths in the sequence dictionaries.  Were these segment files generated with the same reference?");
            case NON_CANONICAL_HUMAN_ORDER:
                logger.warn("Input files contain sequence dictionaries from human reference, but the sorting is not canonical.  Downstream errors are possible.");
                return;
            case OUT_OF_ORDER:
                throw new UserException.BadInput("Input files have different sequence dictionary ordering.  The risk of errors downstream is too high to continue.");

            // These are okay scenarios for this tool.
            case SUPERSET:
            case COMMON_SUBSET:
            case IDENTICAL:
            case NO_COMMON_CONTIGS:
            case DIFFERENT_INDICES:
                return;

            default:
                logger.warn("Received an unrecognized compatibility result: " + compatibilityResult.toString() + ".  Continuing, but errors may result downstream.");
        }
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
    private List<AnnotatedInterval> annotateCombinedIntervals(final List<AnnotatedInterval> segments1, final List<AnnotatedInterval> segments2,
                                                              final List<Map<String, String>> inputToOutputHeaderMaps,
                                                              final SAMSequenceDictionary dictionary,
                                                              final Consumer<Locatable> progressUpdater) {

        final List<List<AnnotatedInterval>> segmentLists = Arrays.asList(segments1, segments2);

        final List<Locatable> combinedIntervals = IntervalUtils.combineAndSortBreakpoints(
                segmentLists.get(0).stream().map(AnnotatedInterval::getInterval)
                        .collect(Collectors.toList()),
                segmentLists.get(1).stream().map(AnnotatedInterval::getInterval)
                        .collect(Collectors.toList()), dictionary);

        logger.info("Creating map of overlaps...");
        // Create a list of maps where each entry corresponds to overlapping intervals of the combined intervals.
        final List<Map<Locatable, List<AnnotatedInterval>>> combinedIntervalsToSegmentsMaps = segmentLists.stream()
                .map(segs -> IntervalUtils.createOverlapMap(combinedIntervals, segs, dictionary)).collect(Collectors.toList());

        final List<AnnotatedInterval> result = new ArrayList<>();

        logger.info("Annotating...");
        for (final Locatable interval: combinedIntervals) {
            final SortedMap<String, String> intervalAnnotationMap = new TreeMap<>();

            for (int i = 0; i < combinedIntervalsToSegmentsMaps.size(); i ++) {
                final Map<Locatable, List<AnnotatedInterval>> combinedIntervalsToSegmentsMap = combinedIntervalsToSegmentsMaps.get(i);
                final Map<String, String> inputToOutputHeaderMap =  inputToOutputHeaderMaps.get(i);
                final List<AnnotatedInterval> matchingSegments = combinedIntervalsToSegmentsMap.get(interval);
                inputToOutputHeaderMap.forEach((k,v) -> {
                                if (matchingSegments.size() > 1) {
                                    logger.warn(interval + " had more than one segment: " + matchingSegments + " only annotating with the first match.");
                                }
                                if (matchingSegments.size() >= 1) {
                                    intervalAnnotationMap.put(v, matchingSegments.get(0).getAnnotationValueOrDefault(k, ""));
                                } else {
                                    intervalAnnotationMap.put(v, "");
                                }
                    });
            }

            result.add(new AnnotatedInterval(new SimpleInterval(interval), intervalAnnotationMap));
            progressUpdater.accept(interval);
        }
        return result;
    }
}