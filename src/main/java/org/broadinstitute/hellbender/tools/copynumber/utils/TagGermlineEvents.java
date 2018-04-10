package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.utils.germlinetagging.SimpleGermlineTagger;

import java.io.File;
import java.util.List;

@CommandLineProgramProperties(
        oneLineSummary = "(EXPERIMENTAL) Do a simplistic tagging of germline events in a tumor segment file.",
        summary = "(EXPERIMENTAL) A tool for tagging possible germline events in a tumor segment file.  The algorithm used is very simplistic.  Segments called as amplified or deleted in the normal are matched by breakpoints (+/- some padding).",
        programGroup = CopyNumberProgramGroup.class)
@ExperimentalFeature
public class TagGermlineEvents extends GATKTool{

    final static public String MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME = "called-matched-normal-seg-file";
    final static public int DEFAULT_GERMLINE_TAG_PADDING_IN_BP = 1000;
    final static public String PADDING_IN_BP_LONG_NAME = "endpoint-padding";
    final static public String GERMLINE_TAG_HEADER = "POSSIBLE_GERMLINE";
    final static public String INPUT_CALL_HEADER = "input-call-header";

    @Argument(
            doc = "Input tumor (called) segment file -- the output will be a copy of this file with the additional information.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME
    )
    private File tumorSegmentFile;

    @Argument(
            doc = "Matched normal called segment file.",
            fullName = MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME
    )
    private File calledNormalSegmentFile;

    @Argument(
            doc = "Output TSV file identical to the tumor segment file, but with additional germline tag column (" + GERMLINE_TAG_HEADER + ").",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    private File outputFile;

    @Argument(
            doc = "Amount of padding, in bases, to allow a breakpoint match.",
            fullName = PADDING_IN_BP_LONG_NAME)
    private int paddingInBp = DEFAULT_GERMLINE_TAG_PADDING_IN_BP;

    @Argument(
            doc = "Column header to use for the call (amplified, neutral, or deleted).",
            fullName = INPUT_CALL_HEADER, optional = true)
    private String callColumnHeader = "CALL";

    @Override
    public void traverse() {
        final AnnotatedIntervalCollection tumorAnnotatedIntervalCollection = AnnotatedIntervalCollection.create(tumorSegmentFile.toPath(), null);
        final List<AnnotatedInterval> initialTumorSegments = tumorAnnotatedIntervalCollection.getRecords();
        final List<AnnotatedInterval> initialNormalSegments = AnnotatedIntervalCollection.create(calledNormalSegmentFile.toPath(), null).getRecords();

        final List<AnnotatedInterval> tumorSegments = SimpleGermlineTagger.tagTumorSegmentsWithGermlineActivity(
                initialTumorSegments, initialNormalSegments, callColumnHeader, getBestAvailableSequenceDictionary(),
                GERMLINE_TAG_HEADER, paddingInBp);

        final List<String> finalAnnotations = tumorAnnotatedIntervalCollection.getAnnotations();
        finalAnnotations.add(GERMLINE_TAG_HEADER);
        finalAnnotations.sort(String::compareTo);

        final AnnotatedIntervalCollection finalCollection =
                AnnotatedIntervalCollection.create(tumorSegments,
                        tumorAnnotatedIntervalCollection.getSamFileHeader(), finalAnnotations);
        finalCollection.write(outputFile);
    }

    @Override
    public boolean requiresReference() {
        return true;
    }
}
