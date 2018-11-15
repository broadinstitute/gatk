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
        summary = "(EXPERIMENTAL) A tool for tagging possible germline events in a tumor segment file.  The algorithm used is very simplistic.  Segments called as amplified or deleted in the normal are matched by breakpoints (+/- some padding) or reciprocal overlap.",
        programGroup = CopyNumberProgramGroup.class)
@ExperimentalFeature
public class TagGermlineEvents extends GATKTool{

    final static public String MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME = "called-matched-normal-seg-file";
    final static public int DEFAULT_GERMLINE_TAG_PADDING_IN_BP = 1000;
    final static public String PADDING_IN_BP_LONG_NAME = "endpoint-padding";
    final static public String GERMLINE_TAG_HEADER = "POSSIBLE_GERMLINE";
    final static public String INPUT_CALL_HEADER = "input-call-header";
    final static public String HI_MAF_HEADER = "input-high-maf-header";
    final static public String LO_MAF_HEADER = "input-low-maf-header";
    final static public String DO_NOT_CHECK_CNLOH = "skip-normal-cnloh";

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
            fullName = INPUT_CALL_HEADER,
            optional = true)
    private String callColumnHeader = "CALL";

    @Argument(
            doc = "Reciprocal threshold to match a segment (when endpoints do not match).  Over this threshold will cause a tumor segment to be tagged.",
            fullName = "reciprocal-threshold",
            optional = true)
    private double reciprocalThreshold = 0.75;

    @Argument(
            doc = "Maximum maf to match a segment by the minor allelic fraction (MAF) when the call in the normal is copy neutral.  Under this threshold will allow a tumor segment to be tagged.",
            fullName = "maf-max-threshold",
            minValue = 0.0, maxValue = 0.5,
            optional = true)
    private double mafMaxThreshold = 0.47;

    @Argument(
            doc = "Column header to use for the high MAF estimate.",
            fullName = HI_MAF_HEADER,
            optional = true)
    private String mafHiColumnHeader = "MINOR_ALLELE_FRACTION_POSTERIOR_90";

    @Argument(
            doc = "Column header to use for the low MAF estimate.",
            fullName = LO_MAF_HEADER,
            optional = true)
    private String mafLoColumnHeader = "MINOR_ALLELE_FRACTION_POSTERIOR_10";

    @Argument(
            doc = "Disable attempt to find CNLoH events in the matched normal.",
            fullName = DO_NOT_CHECK_CNLOH,
            optional = true)
    private Boolean isNoCnLoHCheck = false;  // Default to enable the check.


    @Override
    public void traverse() {
        final AnnotatedIntervalCollection tumorAnnotatedIntervalCollection = AnnotatedIntervalCollection.create(tumorSegmentFile.toPath(), null);
        final List<AnnotatedInterval> initialTumorSegments = tumorAnnotatedIntervalCollection.getRecords();
        final List<AnnotatedInterval> initialNormalSegments = AnnotatedIntervalCollection.create(calledNormalSegmentFile.toPath(), null).getRecords();

        List<AnnotatedInterval> tumorSegments = null;
        if (isNoCnLoHCheck) {
            tumorSegments = SimpleGermlineTagger.tagTumorSegmentsWithGermlineActivity(
                    initialTumorSegments, initialNormalSegments, callColumnHeader, getBestAvailableSequenceDictionary(),
                    GERMLINE_TAG_HEADER, paddingInBp, reciprocalThreshold);
        } else {
            tumorSegments = SimpleGermlineTagger.tagTumorSegmentsWithGermlineActivity(
                    initialTumorSegments, initialNormalSegments, callColumnHeader, getBestAvailableSequenceDictionary(),
                    GERMLINE_TAG_HEADER, paddingInBp, reciprocalThreshold, mafMaxThreshold, mafHiColumnHeader,
                    mafLoColumnHeader);
        }

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
