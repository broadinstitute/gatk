package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegionCollection;
import org.broadinstitute.hellbender.tools.copynumber.utils.germlinetagging.SimpleGermlineTagger;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(
        oneLineSummary = "(EXPERIMENTAL) Do a simplistic tagging of germline events in a tumor segment file.",
        summary = "(EXPERIMENTAL) A tool for tagging possible germline events in a tumor segment file.  The algorithm used is very simplistic.  Segments called as amplified or deleted in the normal are matched by breakpoints (+/- some padding).",
        programGroup = CopyNumberProgramGroup.class)
@BetaFeature
public class TagGermlineEvents extends GATKTool{

    final static public String MATCHED_NORMAL_SEGMENT_FILE_SHORT_NAME = "CMNSeg";
    final static public String MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME = "CalledMatchedNormalSegFile";
    final static public int DEFAULT_GERMLINE_TAG_PADDING_IN_BP = 1000;
    final static public String GERMLINE_TAG_HEADER = "POSSIBLE_GERMLINE";

    @Argument(
            doc = "Input tumor (called) segment file -- the output will be a copy of this file with the additional information.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME
    )
    private File tumorSegmentFile;

    @Argument(
            doc = "Matched normal called segment file.",
            fullName = MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME,
            shortName = MATCHED_NORMAL_SEGMENT_FILE_SHORT_NAME
    )
    private File calledNormalSegmentFile;

    @Argument(
            doc = "Output TSV file identical to the tumor segment file, but with additional germline tag column (" + GERMLINE_TAG_HEADER + ").",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private File outputFile;

    @Argument(
            doc = "Amount of padding to allow a breakpoint match.",
            shortName = "pbp",
            fullName = "paddingInBp")
    private int paddingInBp = DEFAULT_GERMLINE_TAG_PADDING_IN_BP;


    @Override
    public void traverse() {
        final SimpleAnnotatedGenomicRegionCollection tumorSimpleAnnotatedGenomicRegionCollection = SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(tumorSegmentFile);
        final List<SimpleAnnotatedGenomicRegion> initialTumorSegments = tumorSimpleAnnotatedGenomicRegionCollection.getRecords();
        final List<SimpleAnnotatedGenomicRegion> initialNormalSegments = SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(calledNormalSegmentFile).getRecords();

        final String callHeader = CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn.CALL.toString();

        final List<SimpleAnnotatedGenomicRegion> tumorSegments = SimpleGermlineTagger.tagTumorSegmentsWithGermlineActivity(
                initialTumorSegments, initialNormalSegments, callHeader, getBestAvailableSequenceDictionary(),
                GERMLINE_TAG_HEADER, paddingInBp);

        final SimpleAnnotatedGenomicRegionCollection finalCollection =
                SimpleAnnotatedGenomicRegionCollection.createCollectionFromExistingCollection(tumorSegments,
                        tumorSimpleAnnotatedGenomicRegionCollection,
                        new ArrayList<>(tumorSegments.get(0).getAnnotations().keySet()));
        finalCollection.write(outputFile);
    }

    @Override
    public boolean requiresReference() {
        return true;
    }
}
