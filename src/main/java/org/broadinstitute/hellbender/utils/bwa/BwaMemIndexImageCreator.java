package org.broadinstitute.hellbender.utils.bwa;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.BwaMemUtilitiesProgramGroup;


/**
 * Simply creates the reference index image file.
 *
 */
@CommandLineProgramProperties(
        summary = "Counts filtered reads at het sites for allele specific expression estimate",
        oneLineSummary = "Generates table of filtered base counts at het sites for allele specific expression",
        programGroup = BwaMemUtilitiesProgramGroup.class
)
public final class BwaMemIndexImageCreator extends CommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
              shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
              doc = "Output file (if not provided, defaults to STDOUT)",
              optional = false)
    private String referenceFastaLoc = null;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              doc = "Output file (if not provided, defaults to STDOUT)",
              optional = false)
    private String referenceIndexImageOutputLoc = null;

    @Override
    protected final Object doWork() {
        BwaMemIndex.createIndexImage(referenceFastaLoc, referenceIndexImageOutputLoc);
        return null;
    }
}
