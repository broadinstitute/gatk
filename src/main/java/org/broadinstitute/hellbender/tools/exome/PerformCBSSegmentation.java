package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExomeAnalysisProgramGroup;
import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;

@CommandLineProgramProperties(
        summary = "Segment genomic data into regions of constant copy-ratio",
        oneLineSummary = "Segment genomic data into regions of constant copy-ratio",
        programGroup = ExomeAnalysisProgramGroup.class
)
public final class PerformCBSSegmentation extends CommandLineProgram {

    public static final String SAMPLE_NAME_LONG_NAME = "sampleName";
    public static final String SAMPLE_NAME_SHORT_NAME = "S";

    public static final String TARGETS_FILE_LONG_NAME = "targets";
    public static final String TARGETS_FILE_SHORT_NAME = "T";

    public static final String SEGMENT_FILE_LONG_NAME = StandardArgumentDefinitions.OUTPUT_LONG_NAME;
    public static final String SEGMENT_FILE_SHORT_NAME = StandardArgumentDefinitions.OUTPUT_SHORT_NAME;

    public static final String MIN_LOG_VALUE_LONG_NAME = "minLogValue";
    public static final String MIN_LOG_VALUE_SHORT_NAME = "M";

    @Argument(
            doc = "Name of the sample being segmented",
            shortName = SAMPLE_NAME_SHORT_NAME,
            fullName = SAMPLE_NAME_LONG_NAME,
            optional = false
    )
    protected String sampleName;

    @Argument(
            doc = "Genomic targets file",
            shortName = TARGETS_FILE_SHORT_NAME,
            fullName =  TARGETS_FILE_LONG_NAME,
            optional = false
    )
    protected String tangentFile;

    @Argument(
            doc = "Full path to the outputted segment file",
            shortName = SEGMENT_FILE_SHORT_NAME,
            fullName = SEGMENT_FILE_LONG_NAME,
            optional = false
    )
    protected String outFile;

    @Argument(
            doc = "Any values under this threshold will be set to it",
            shortName = MIN_LOG_VALUE_SHORT_NAME,
            fullName = MIN_LOG_VALUE_LONG_NAME,
            optional = true
    )
    protected float minLogValue;

    @Override
    protected Object doWork() {
        applySegmentation(sampleName, tangentFile, outFile, minLogValue);
        return "Success";
    }

    private void applySegmentation(String sampleName, String tangentFile, String outFile, float minLogValue){
        RCBSSegmenter.writeSegmentFile(sampleName, tangentFile, outFile, minLogValue);
    }
}
