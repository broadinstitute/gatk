package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegionCollection;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegionUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;

import java.io.File;
import java.util.List;

@CommandLineProgramProperties(
        oneLineSummary = "(EXPERIMENTAL) Merge annotated genomic regions based entirely on contig and annotation value.",
        summary = "(EXPERIMENTAL) Merge annotated genomic regions based entirely on contig and annotation value.  Column header order is not guaranteed to be preserved. Reference is required and will superseded any sequence dictionary in the given seg/region files.  Sequence dictionary is optional on the input file, but will be included on the output.  Conflicts of annotations is resolved by sorting the values and inserting a delimiter.  THIS TOOL IS TOTALLY UNSUPPORTED",
        programGroup = CopyNumberProgramGroup.class)

@BetaFeature
public class MergeAnnotatedRegions extends GATKTool {
    @Argument(
            doc = "Input segment file or annotated segment file -- touching segments will be merged in the output.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME
    )
    private File segmentFile;

    @Argument(
            doc = "Output TSV file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private File outputFile;

    final static String DEFAULT_SEPARATOR = "__";

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void traverse() {

        // Load the seg file.  Not done with the collection class, in order to keep the sequence dictionary optional on the input.
        final List<SimpleAnnotatedGenomicRegion> initialSegments = SimpleAnnotatedGenomicRegion.readAnnotatedRegions(segmentFile);

        // Sort the locatables
        final List<SimpleAnnotatedGenomicRegion> segments = IntervalUtils.sortLocatablesBySequenceDictionary(initialSegments,
                getBestAvailableSequenceDictionary());

        // Perform the actual merging
        // For each segment, see if we have more than one overlap.  If we do, merge to create a new segment.
        final List<SimpleAnnotatedGenomicRegion> finalSegments = SimpleAnnotatedGenomicRegionUtils.mergeRegions(segments,
                getBestAvailableSequenceDictionary(), DEFAULT_SEPARATOR, l -> progressMeter.update(l));

        final SimpleAnnotatedGenomicRegionCollection collection = SimpleAnnotatedGenomicRegionCollection.create(finalSegments,
                getBestAvailableSequenceDictionary(), finalSegments.get(0).getAnnotationNames());
        collection.write(outputFile);
    }
}
