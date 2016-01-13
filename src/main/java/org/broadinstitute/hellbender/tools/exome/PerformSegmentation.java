package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Segment genomic data into regions of constant copy-ratio",
        oneLineSummary = "Segment genomic data into regions of constant copy-ratio",
        programGroup = CopyNumberProgramGroup.class
)
public final class PerformSegmentation extends CommandLineProgram {

    public static final String TARGET_WEIGHT_FILE_LONG_NAME= "targetWeights";
    public static final String TARGET_WEIGHT_FILE_SHORT_NAME = "tw";

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

    @Argument(
            doc = "File with target weights.  This is the 1/var(post-projected targets for each normal).  " +
                    "Listed one value per line in plain text.  Values of zero or less, Nan, Inf, and -Inf are not " +
                    "acceptable.  Must have the same number of values as there are in the tangentFile.",
            shortName = TARGET_WEIGHT_FILE_SHORT_NAME,
            fullName = TARGET_WEIGHT_FILE_LONG_NAME,
            optional = true
    )
    protected File weightFile = null;

    @Override
    protected Object doWork() {
        applySegmentation(sampleName, tangentFile, outFile);
        return "Success";
    }

    private void applySegmentation(String sampleName, String tangentFile, String outFile) {
        RCBSSegmenter.writeSegmentFile(sampleName, tangentFile, outFile, log, 2, weightFile);
    }
}
