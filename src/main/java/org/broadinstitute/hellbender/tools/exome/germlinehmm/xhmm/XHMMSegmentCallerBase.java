package org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm;

import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm.XHMMArgumentCollection;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm.XHMMModel;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

/**
 * Parent class for those tools that make CNV segment calls based on a {@link XHMMModel}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public abstract class XHMMSegmentCallerBase extends CommandLineProgram {

    public static final String ZSCORE_DIMENSION_FULL_NAME = "standardizeBy";
    public static final String ZSCORE_DIMENSION_SHORT_NAME = "standardizeBy";

    @ArgumentCollection
    protected XHMMArgumentCollection modelArguments = new XHMMArgumentCollection();

    /**
     * Z-scores transformation dimensions.
     */
    public enum ZScoreDimension {

        /**
         * Input coverage must be transformed in z-scores sample by sample.
         */
        SAMPLE,

        /**
         * Input coverage must be transformed in z-scores target by target.
         */
        TARGET,

        /**
         * No transformation must be performed.
         */
        NONE
    }

    @Argument(
            doc = "In what dimension apply the zscore transformation",
            fullName = ZSCORE_DIMENSION_FULL_NAME,
            shortName = ZSCORE_DIMENSION_SHORT_NAME,
            optional = true
    )
    protected ZScoreDimension zscoreDimension = ZScoreDimension.NONE;

    @Argument(
            doc = "Input counts",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    protected File inputFile;

    @Argument(
            doc = "Output file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected File outputFile;

    @ArgumentCollection
    protected final TargetArgumentCollection targetArguments = new TargetArgumentCollection(() -> inputFile);

    /* HMM results */
    protected List<ForwardBackwardAlgorithm.Result<XHMMEmissionData, Target, CopyNumberTriState>>
            sampleForwardBackwardResults;
    protected List<List<CopyNumberTriState>> sampleBestPaths;

    @Override
    protected final Object doWork() {
        final XHMMModel model = modelArguments.createModel();
        final TargetCollection<Target> targets = targetArguments.readTargetCollection(false);
        final ReadCountCollection inputCounts = readAndSortCountsByTargetCoordinates(targets);
        applyZScoreTransformation(inputCounts);
        checkForMissingTargetsInInputCounts(targets, inputCounts);
        runForwardBackwardAndViterbi(model, targets, inputCounts);
        openOutput(outputFile, model, targets, inputCounts);
        makeCalls(model, targets, inputCounts);
        closeOutput(outputFile);
        return "SUCCESS";
    }

    /**
     * Run forward-backward algorithm and Viterbi algorithm on each sample
     *
     * @param model an instance of {@link XHMMModel}
     * @param targets input target collection
     * @param inputCounts input read count collection
     */
    private void runForwardBackwardAndViterbi(final XHMMModel model, final TargetCollection<Target> targets,
                                              final ReadCountCollection inputCounts) {
        sampleForwardBackwardResults = new ArrayList<>(inputCounts.columnNames().size());
        sampleBestPaths = new ArrayList<>(inputCounts.columnNames().size());
        for (int sampleIndex = 0; sampleIndex < inputCounts.columnNames().size(); sampleIndex++) {
            final List<XHMMEmissionData> emissionData = DoubleStream.of(inputCounts.counts().getColumn(sampleIndex))
                    .mapToObj(XHMMEmissionData::new)
                    .collect(Collectors.toList());
            final ForwardBackwardAlgorithm.Result<XHMMEmissionData, Target, CopyNumberTriState> fbResult =
                    ForwardBackwardAlgorithm.apply(emissionData, targets.targets(), model);
            final List<CopyNumberTriState> bestPath = ViterbiAlgorithm.apply(emissionData, targets.targets(), model);
            sampleForwardBackwardResults.add(fbResult);
            sampleBestPaths.add(bestPath);
        }
    }

    /**
     * Opens the output where to dump the calls.
     *
     * @throws UserException.CouldNotCreateOutputFile if something went wrong creating/opening or writing in the file.
     * @throws UserException.BadInput if there inputs (coverage or targets) are somehow invalid.
     * @throws IllegalArgumentException if any of the parameters is {@code null}.
     */
    protected abstract void openOutput(final File outputFile, final XHMMModel model,
                                       final TargetCollection<Target> targets, final ReadCountCollection inputCounts);

    /**
     * Closes the output.
     *
     * @throws UserException.CouldNotCreateOutputFile if something went wrong.
     */
    protected abstract void closeOutput(final File outputFile);

    /**
     * Does the actual CN segment calling, outputting results to the output file.
     *
     * @param model the CN segment model.
     * @param targets the target under analysis.
     * @param inputCounts the input coverage values already processed.
     * @throws IllegalArgumentException if any of the parameters is {@code null} or otherwise invalid.
     * @throws IllegalStateException if the tool is not ready for calling (e.g. the output file was not properly opened).
     * @throws UserException if anything went wrong due to something under the user's control.
     * @throws GATKException if anything went wrong due to something outside the user's control.
     */
    protected abstract void makeCalls(final XHMMModel model, final TargetCollection<Target> targets,
                                      final ReadCountCollection inputCounts);

    /**
     * Transform read counts to z-scores.
     * @param inputCounts the input read-counts, modified in-situ.
     */
    private void applyZScoreTransformation(final ReadCountCollection inputCounts) {
        final RealMatrix counts = inputCounts.counts();
        switch (zscoreDimension) {
            case SAMPLE:
                standardizeBySample(counts);
                break;
            case TARGET:
                standardizeByTarget(counts);
                break;
        }
    }

    /**
     * Standardize read counts (per-target).
     * Note: modification is done in-place.
     *
     * @param counts original read counts
     */
    private void standardizeByTarget(final RealMatrix counts) {
        final double[] rowMeans = GATKProtectedMathUtils.rowMeans(counts);
        final double[] rowStdDev = GATKProtectedMathUtils.rowStdDevs(counts);

        counts.walkInColumnOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return (value - rowMeans[row]) / rowStdDev[row];
            }

        });
    }

    /**
     * Standardize read counts (per-sample).
     * Note: modification is done in-place.
     *
     * @param counts original read counts
     */
    private void standardizeBySample(final RealMatrix counts) {
        final double[] columnMeans = GATKProtectedMathUtils.columnMeans(counts);
        final double[] columnStdDev = GATKProtectedMathUtils.columnStdDevs(counts);

        counts.walkInColumnOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return (value - columnMeans[column]) / columnStdDev[column];
            }

        });
    }

    /**
     * Read the input read counts table on the specified targets and order them with respect to the
     * provided list of targets.
     *
     * @param targets input target collection
     * @return sorted read counts
     */
    private ReadCountCollection readAndSortCountsByTargetCoordinates(final TargetCollection<Target> targets) {

        // read the input read-counts.
        final ReadCountCollection originalReadCounts = readInputCounts(targets);

        // Get the input read count target positions in the input targets collection.
        final int[] sortedPositions = originalReadCounts.targets().stream()
                .mapToInt(t -> targets.index(t.getName())).toArray();

        // Check whether the input is in the same order as the targets in the input target collection.
        // If not we rearrange them accordingly.
        for (int i = 1, previousIndex = -1; i < sortedPositions.length; i++) {
            if (sortedPositions[i] == -1) {
                return originalReadCounts.arrangeTargets(targets.targets());
            } else {
                if (previousIndex >= sortedPositions[i]) {
                    return originalReadCounts.arrangeTargets(targets.targets());
                }
                previousIndex = sortedPositions[i];
            }
        }
        return originalReadCounts;
    }

    /**
     * Parse the input read count table on specified targets (missing targets will be ignored)
     *
     * @param targets input target collection
     * @return read counts
     */
    private ReadCountCollection readInputCounts(final TargetCollection<Target> targets) {
        try {
            return ReadCountCollectionUtils.parse(inputFile, targets, true);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(inputFile, ex);
        }
    }

    /**
     * Checks the input read count table against input target list for possibly missing targets and generate
     * warning messages accordingly.
     *
     * @param targets input target collection
     * @param inputCounts input read count collection
     */
    private void checkForMissingTargetsInInputCounts(final TargetCollection<Target> targets, final ReadCountCollection inputCounts) {
        if (inputCounts.targets().size() != targets.targetCount()) {
            logger.warn(String.format("The input counts are missing some of the targets required (%d): e.g. %s",
                    inputCounts.targets().size() - targets.targetCount(), targets.targets().stream()
                            .filter(t -> !inputCounts.targets().contains(t))
                            .limit(5)));
        }
    }
}
