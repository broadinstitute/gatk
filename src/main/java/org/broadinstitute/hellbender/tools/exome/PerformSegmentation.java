package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;

@CommandLineProgramProperties(
        summary = "Segment genomic data into regions of constant copy-ratio",
        oneLineSummary = "Segment genomic data into regions of constant copy-ratio",
        programGroup = CopyNumberProgramGroup.class
)
public final class PerformSegmentation extends CommandLineProgram {

    @Argument(
            doc = "Name of the sample being segmented",
            fullName = ExomeStandardArgumentDefinitions.SAMPLE_LONG_NAME,
            optional = false
    )
    protected String sampleName;

    @Argument(
            doc = "Genomic targets file",
            shortName = ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME,
            optional = false
    )
    protected String tangentFile;

    @Argument(
            doc = "Full path to the outputted segment file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false
    )
    protected String outFile;

    @Argument(
            doc = "If input data has had a log2 transform applied",
            shortName = ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.LOG2_LONG_NAME,
            optional = true
    )
    protected Boolean log = false;

    @Override
    protected Object doWork() {
        applySegmentation(sampleName, tangentFile, outFile);
        return "Success";
    }

    private void applySegmentation(String sampleName, String tangentFile, String outFile){
        RCBSSegmenter.writeSegmentFile(sampleName, tangentFile, outFile, log);
    }
}
