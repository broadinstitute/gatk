package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.programgroups.ReferenceProgramGroup;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

/**
 * Simply creates the reference index image file.
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Creates the image file for use by BwaMemAligner",
        oneLineSummary = "Creates the image file for use by BwaMemAligner",
        programGroup = ReferenceProgramGroup.class
)
public final class BwaMemIndexImageCreator extends CommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "Input reference fasta file location.")
    private String referenceFastaLoc = null;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output reference index image file (ending in \".img\").",
            optional = true)
    private String referenceIndexImageOutputLoc = null;

    @Override
    protected final Object doWork() {
        if (referenceIndexImageOutputLoc == null) {
            referenceIndexImageOutputLoc = referenceFastaLoc + ".img";
        }
        BwaMemIndex.createIndexImageFromFastaFile(referenceFastaLoc, referenceIndexImageOutputLoc);
        return null;
    }
}
