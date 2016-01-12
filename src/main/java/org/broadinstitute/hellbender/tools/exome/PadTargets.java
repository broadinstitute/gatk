package org.broadinstitute.hellbender.tools.exome;


import htsjdk.tribble.bed.BEDFeature;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Creates a new target BED file with targets extended on both sides by the specified number of bases.  IMPORTANT:  This tool will only preserve contig, start, end, and name columns.",
        oneLineSummary = "Create a new target file with padded targets.",
        programGroup = CopyNumberProgramGroup.class
)
public final class PadTargets extends CommandLineProgram {

    protected static final String TARGET_FILE_FULL_NAME = "targets";
    protected static final String TARGET_FILE_SHORT_NAME = "T";

    protected static final String PADDING_SHORT_NAME = "p";
    protected static final String PADDING_FULL_NAME = "padding";

    @Argument(
            doc = "File containing the targets (BED) for padding.  Should have no existing overlap.",
            shortName = TARGET_FILE_SHORT_NAME,
            fullName = TARGET_FILE_FULL_NAME,
            optional = true
    )
    protected File targetFile = null;

    @Argument(
            doc = "Output target BED file name.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false
    )
    protected File outFile = null;

    @Argument(
            doc = "Amount of bases to be added on either side of all targets.",
            shortName = PADDING_SHORT_NAME,
            fullName = PADDING_FULL_NAME,
            optional = true
    )
    protected Integer padding = 250;

    @Override
    protected Object doWork() {

        final TargetCollection<? extends BEDFeature> inputTargetCollection = TargetUtils.readTargetFile(targetFile);

        final TargetCollection<Target> paddedTargetCollection = TargetPadder.padTargetsFromBEDFeatures(inputTargetCollection, padding);

        TargetCoverageUtils.writeTargetsAsBed(outFile, paddedTargetCollection.targets());

        return "SUCCESS";
    }
}
