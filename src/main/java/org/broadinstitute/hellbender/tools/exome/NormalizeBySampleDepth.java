package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.IOException;

/**
 * Calculates basic statistics for coverage per targets and per sample.
 * <p>
 *   This tool takes a read counts file and outputs a read counts file in which coverage at each target
 *   and sample is divided by the average coverage of that sample.
 * </p>
 * <p>
 *   The input and output format for the coverage file is described in {@link ReadCountCollectionUtils}.
 * </p>
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Normalize read counts by samples' average coverage.",
        oneLineSummary = "Normalize read counts by samples' average coverage",
        programGroup = CopyNumberProgramGroup.class
)
public class NormalizeBySampleDepth extends CommandLineProgram {
    public static final String WEIGHTED_AVERAGE_LONG_NAME = "weightedAverage";
    public static final String WEIGHTED_AVERAGE_SHORT_NAME = "w";

    @Argument(
            doc = "Input coverage file",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            optional = false
    )
    protected File inputFile;

    @Argument(
            doc = "Output coverage file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected File outputFile;

    @Argument(
            doc = "Weight average with weights proportional to target sizes",
            fullName = WEIGHTED_AVERAGE_LONG_NAME,
            shortName = WEIGHTED_AVERAGE_SHORT_NAME,
            optional = true
    )
    protected boolean weightByTargetSize = false;

    @Override
    public Object doWork() {
        final ReadCountCollection inputCounts;
        try {
            inputCounts = ReadCountCollectionUtils.parse(inputFile);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, "cannot reach or read the input target coverage file");
        }

        final ReadCountCollection normalizedCounts;
        try {
            normalizedCounts = inputCounts.normalizeByColumnAverages(weightByTargetSize);
        } catch (final IllegalArgumentException e) {
            throw new UserException.BadInput("Weighting by target size requested but input read normalizedCounts lack target intervals.");
        }

        try {
            ReadCountCollectionUtils.write(outputFile, normalizedCounts);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "cannot  write to the given output file");
        }

        return "SUCCESS";
    }
}
