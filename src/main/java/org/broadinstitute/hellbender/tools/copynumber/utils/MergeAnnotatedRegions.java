package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;

import java.io.File;
import java.util.List;

@CommandLineProgramProperties(
        oneLineSummary = "(EXPERIMENTAL) Merge annotated genomic regions based entirely on contig and annotation value.",
        summary = "(EXPERIMENTAL) Merge annotated genomic regions based entirely on contig and annotation value.  Column header order is not guaranteed to be preserved. Reference is required and will superseded any sequence dictionary in the given seg/region files.  Sequence dictionary is optional on the input file, but will be included on the output.  Conflicts of annotations is resolved by sorting the values and inserting a delimiter.  THIS TOOL IS TOTALLY UNSUPPORTED",
        programGroup = CopyNumberProgramGroup.class)

@ExperimentalFeature
public class MergeAnnotatedRegions extends GATKTool {
    @Argument(
            doc = "Input segment file or annotated segment file -- touching segments will be merged in the output.",
            fullName = CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME
    )
    private File segmentFile;

    @Argument(
            doc = "Output TSV file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private File outputFile;

    final static String DEFAULT_SEPARATOR = "__";

    @Override
    public void traverse() {

        // Load the seg file.  Not done with the collection class, in order to keep the sequence dictionary optional on the input.
        final AnnotatedIntervalCollection annotatedIntervalCollection = AnnotatedIntervalCollection.create(segmentFile.toPath(), null);
        final List<AnnotatedInterval> initialSegments = annotatedIntervalCollection.getRecords();

        // Sort the locatables
        final List<AnnotatedInterval> segments = IntervalUtils.sortLocatablesBySequenceDictionary(initialSegments,
                getBestAvailableSequenceDictionary());

        // Perform the actual merging
        // For each segment, see if we have more than one overlap.  If we do, merge to create a new segment.
        final List<AnnotatedInterval> finalSegments = AnnotatedIntervalUtils.mergeRegions(segments,
                getBestAvailableSequenceDictionary(), DEFAULT_SEPARATOR, l -> progressMeter.update(l));

        final AnnotatedIntervalCollection collection =
                AnnotatedIntervalCollection.create(finalSegments, annotatedIntervalCollection.getSamFileHeader(),
                        annotatedIntervalCollection.getAnnotations());
        collection.write(outputFile);
    }
}
