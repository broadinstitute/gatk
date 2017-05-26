package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.util.List;

/**
 * Extend target intervals on either side without overlapping consecutive intervals.
 *
 * <p>Input and output target files are both in the format described in {@link TargetWriter}.</p>
 * <p>Any overlapping padding is split equally between upstream and downstream targets.</p>
 *
 * <h3>Example</h3>
 *
 * <p>For whole exome sequencing (WES) targets: </p>
 *
 * <pre>
 * java -Xmx1g -jar $gatk_jar PadTargets \
 *   --targets targets.tsv \
 *   --output targets.padded.tsv
 * </pre>
 *
 * <p>
 *     By default, the tool pads targets 250 basepairs on either side of each target.
 *     This amount of padding was empirically determined to increase sensitivity for CNV calling
 *     for Illumina short-read whole exome sequencing (WES) data similar to TCGA Project data.
 * </p>
 */
@CommandLineProgramProperties(
        summary = "Creates a new target file with targets extended on both sides by the specified number of bases.  IMPORTANT:  This tool will only preserve contig, start, end, and name columns.",
        oneLineSummary = "Create a new target file with padded targets",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class PadTargets extends CommandLineProgram {

    protected static final String TARGET_FILE_FULL_NAME = ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME;
    protected static final String TARGET_FILE_SHORT_NAME = ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME;

    protected static final String PADDING_SHORT_NAME = "p";
    protected static final String PADDING_FULL_NAME = "padding";

    @Argument(
            doc = "File containing the targets for padding.  Should have no existing overlap.",
            fullName = TARGET_FILE_FULL_NAME,
            shortName = TARGET_FILE_SHORT_NAME,
            optional = false
    )
    protected File targetFile = null;

    @Argument(
            doc = "Output target file name.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected File outFile = null;

    @Argument(
            doc = "Amount of bases to be added on either side of all targets.",
            fullName = PADDING_FULL_NAME,
            shortName = PADDING_SHORT_NAME,
            optional = true
    )
    protected Integer padding = 250;

    @Override
    protected Object doWork() {
        if (padding < 0) {
            throw new UserException.BadInput("Padding parameter must be >= 0");
        }

        final List<Target> inputTargets = TargetTableReader.readTargetFile(targetFile);
        final List<Target> paddedTargets = TargetPadder.padTargets(inputTargets, padding);
        TargetWriter.writeTargetsToFile(outFile, paddedTargets);
        return "SUCCESS";
    }
}
