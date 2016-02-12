package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.DefaultRealMatrixPreservingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

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
        summary = "Normalizes read counts by samples' average coverage.",
        oneLineSummary = "Normalizes read counts by samples' average coverage.",
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
            doc = "Weighted average with weights proportional to target sizes",
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

        final RealMatrix counts = inputCounts.counts();

        //kludgy way to deep copy input targets list
        final List<Target> targets = inputCounts.targets().stream()
                .map(t -> new Target(t.getName(), t.getInterval())).collect(Collectors.toList());
        final List<String> columnNames = new ArrayList<>(inputCounts.columnNames());

        final int[] targetSizes;
        if (weightByTargetSize) {
            try {
                targetSizes = targets.stream().mapToInt(t -> length(t)).toArray();
            } catch (final IllegalStateException e) {
                throw new UserException.BadInput("Weighting by target size requested but input read counts lack target intervals.");
            }
        } else {
            targetSizes = new int[targets.size()];
            Arrays.fill(targetSizes, 1);
        }
        final int totalTargetSize = Arrays.stream(targetSizes).sum();

        final double[] sampleWeightedAverages = new double[counts.getColumnDimension()]; // elements initialized to 0.0

        counts.walkInOptimizedOrder(new DefaultRealMatrixPreservingVisitor() {
            @Override
            public void visit(final int target, final int sample, final double coverage) {
                sampleWeightedAverages[sample] += targetSizes[target] * coverage / totalTargetSize;
            }
        });

        counts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int target, final int sample, final double coverage) {
                return coverage / sampleWeightedAverages[sample];
            }
        });

        final ReadCountCollection normalizedCounts = new ReadCountCollection(SetUniqueList.setUniqueList(targets),
                SetUniqueList.setUniqueList(columnNames), counts);

        try {
            ReadCountCollectionUtils.write(outputFile, normalizedCounts);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "cannot  write to the given output file");
        }

        return "SUCCESS";
    }

    private static int length(final Target target) { return target.getEnd() - target.getStart() + 1; }
}
