package org.broadinstitute.hellbender.tools.exome.coveragestats;

import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.ReadCountRecord;
import org.broadinstitute.hellbender.tools.exome.ReadCountsReader;
import org.broadinstitute.hellbender.tools.exome.Target;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Calculates basic statistics for coverage per targets and per sample.
 * <p>
 *   This tool generates two output files that contains tables for summary stats
 *   per target ({@link #targetOutputFile} argument) and
 *   per sample ({@link #sampleOutputFile} argument).
 * </p>
 * <p>
 *   The user can indicate only one of the two output files. In this case the other type of
 *   output is omitted.
 * </p>
 * <p>
 *   No output file specified by the user will result in a {@link UserException} to be thrown.
 * </p>
 * <p>
 *   The input format for the coverage file is described in {@link ReadCountCollectionUtils}.
 * </p>
 * <p>
 *   The output format for the per-target output file is described in {@link TargetCoverageStatsWriter}
 *   whereas the output format for the per-sample output file is described in {@link SampleCoverageStatsWriter}.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Calculates the mean and variance of the coverage per target and sample",
        oneLineSummary = "Calculate statistics of the coverage per target and sample",
        programGroup = CopyNumberProgramGroup.class
)
public final class CalculateCoverageStats extends CommandLineProgram {

    /**
     * Full name for the {@link #sampleOutputFile} argument.
     */
    public static final String SAMPLE_OUTPUT_FILE_FULL_NAME = "sampleOutputFile";

    /**
     * Short name for the {@link #sampleOutputFile} argument.
     */
    public static final String SAMPLE_OUTPUT_FILE_SHORT_NAME = "so";

    /**
     * Full name for the {@link #targetOutputFile} argument.
     */
    public static final String TARGET_OUTPUT_FILE_FULL_NAME = "targetOutputFile";

    /**
     * Short name for the {@link #targetOutputFile} argument.
     */
    public static final String TARGET_OUTPUT_FILE_SHORT_NAME = "to";

    @Argument(
        doc = "Input target coverage file",
        fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
        shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
        optional = false
    )
    protected File inputFile;

    @Argument(
        doc = "Output for per sample statistics",
        fullName = SAMPLE_OUTPUT_FILE_FULL_NAME,
        shortName = SAMPLE_OUTPUT_FILE_SHORT_NAME,
        optional = true
    )
    protected File sampleOutputFile;

    @Argument(
        doc = "Output for per target statistics",
        fullName = TARGET_OUTPUT_FILE_FULL_NAME,
        shortName = TARGET_OUTPUT_FILE_SHORT_NAME,
        optional = true
    )
    protected File targetOutputFile;

    /**
     * Holds a list of sample names in the same order they appear in the input coverage count columns.
     */
    private List<String> sampleNames;

    /**
     * Holds the sum of all coverage value per sample. The order in the array corresponds to the sample
     * order in {@link #sampleNames}.
     */
    private double[] sampleSum;

    /**
     * Hodls the sum of all squared coverage values per sample. The order in the array corresponds to the sample
     * order in {@link #sampleNames}.
     */
    private double[] sampleSquaredSum;

    /**
     * Number of targets so far processed.
     */
    private long targetCount;

    /**
     * Writer for the per-target output table.
     */
    private TargetCoverageStatsWriter targetWriter;

    @Override
    public Object doWork() {
        try (final ReadCountsReader reader = new ReadCountsReader(inputFile)) {
            onTraversalStart();
            reader.stream().forEach(this::processCoverage);
            onTraversalDone();
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(inputFile, ex);
        }
        return null;
    }

    private void onTraversalStart() {
        if (sampleOutputFile == null && targetOutputFile == null) {
            throw new UserException.BadArgumentValue(TARGET_OUTPUT_FILE_FULL_NAME, "you must indicate an output target file name or an output sample file");
        }
        final Pair<Boolean, List<String>> inputCoverageContent = checkInputCoverageContent();
        onTraversalStartForSampleOutput(inputCoverageContent);
        onTraversalStartForTargetOutput(inputCoverageContent);
    }

    private void onTraversalStartForTargetOutput(Pair<Boolean, List<String>> inputCoverageContent) {
        if (targetOutputFile != null) {
            try {
                targetWriter = new TargetCoverageStatsWriter(targetOutputFile, inputCoverageContent.getFirst());
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(targetOutputFile, ex);
            }
        }
    }

    private void onTraversalStartForSampleOutput(Pair<Boolean, List<String>> inputCoverageContent) {
        if (sampleOutputFile != null) {
            sampleNames = inputCoverageContent.getSecond();
            sampleSum = new double[sampleNames.size()];
            sampleSquaredSum = new double[sampleNames.size()];
        }
    }

    private Pair<Boolean, List<String>> checkInputCoverageContent() {
        try (final ReadCountsReader reader = new ReadCountsReader(inputFile)) {
            final List<String> countColumnNames = reader.getCountColumnNames();
            if (countColumnNames.isEmpty()) {
                throw new UserException.BadInput("the input coverage does not contain any coverage column");
            }
            return new Pair<>(reader.hasTargetIntervals(), countColumnNames);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(inputFile);
        }
    }

    private void processCoverage(final ReadCountRecord coverage) {
        final double[] values = coverage.getDoubleCounts();
        final Target target = coverage.getTarget();
        processCoverageForSampleOutput(values);
        processCoverageForTargetOutput(values, target);
        targetCount++;
    }

    private void processCoverageForSampleOutput(final double[] coverage) {
        if (sampleOutputFile != null) {
            for (int i = 0; i < coverage.length; i++) {
                final double value = coverage[i];
                sampleSum[i] += value;
                sampleSquaredSum[i] += value * value;
            }
        }
    }

    private void processCoverageForTargetOutput(final double[] coverage, final Target target) {
        if (targetOutputFile != null) {
            try {
                targetWriter.writeRecord(TargetCoverageStats.fromCoverage(target, coverage));
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(targetOutputFile, ex);
            }
        }
    }

    public Object onTraversalDone() {
        onTraversalDoneForTargetOutput();
        onTraversalDoneForSampleOutput();
        if (targetCount == 0) {
            throw new UserException.BadInput("the input coverage does not contain any targets");
        }
        return null;
    }

    private void onTraversalDoneForSampleOutput() {
        if (sampleOutputFile != null) {
            try (final SampleCoverageStatsWriter writer = new SampleCoverageStatsWriter(sampleOutputFile)) {
                for (int i = 0; i < sampleNames.size(); i++) {

                    writer.writeRecord(SampleCoverageStats.fromSums(sampleNames.get(i), targetCount,
                            sampleSum[i], sampleSquaredSum[i]));
                }
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(sampleOutputFile, ex);
            }
        }
    }

    private void onTraversalDoneForTargetOutput() {
        if (targetWriter != null) {
            try {
                targetWriter.close();
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(targetOutputFile, ex);
            }
        }
    }
}
