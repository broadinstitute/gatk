package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Parent class for those tools that perform calls CNV segment based on a {@link CopyNumberTriStateHiddenMarkovModel}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public abstract class CopyNumberTriStateSegmentCaller extends CommandLineProgram {

    public static final String ZSCORE_DIMENSION_FULL_NAME = "standardizeBy";
    public static final String ZSCORE_DIMENSION_SHORT_NAME = "standardizeBy";

    /**
     * Threshold used to determine best way to calculate log(1- exp(a))
     * based on https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
     */
    private static final double LN_1_M_EXP_THRESHOLD = - Math.log(2);

    /**
     * Maximum tolerated log probability value.
     * <p>Values between this
     * constant and 0.0 as considered as a 0.0.
     * </p>
     * Values above this threshold
     * are considered to indicate a log probability calculation problem that is
     * worth to be reported as an error or warning.
     */
    private static final double MAX_LOG_PROB = 1e-3;

    /**
     * Cached value of ln(10)^-1
     */
    protected static final double INV_LN_10 = 1.0 / Math.log(10);

    /**
     * Maximum reportable output quality score; higher quality scores will
     * capped to this value.
     */
    public static final double MAX_QUAL_SCORE = 199.99999999999;

    public static final double PHRED_SCORE_PRECISION = .1;

    public static final double LOG10_PROB_PRECISION = PHRED_SCORE_PRECISION * .1;

    @ArgumentCollection
    protected CopyNumberTriStateHiddenMarkovModelArgumentCollection modelArguments =
            new CopyNumberTriStateHiddenMarkovModelArgumentCollection();

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

    @Override
    protected final Object doWork() {
        final CopyNumberTriStateHiddenMarkovModel model = modelArguments.createModel();
        final TargetCollection<Target> targets = targetArguments.readTargetCollection(false);
        final ReadCountCollection inputCounts = readAndSortCountsByTargetCoordinates(targets);
        applyZScoreTransformation(inputCounts);
        checkForMissingTargetsInInputCounts(targets, inputCounts);
        openOutput(outputFile, model, targets, inputCounts);
        makeCalls(model, targets, inputCounts);
        closeOutput(outputFile);
        return "SUCCESS";
    }

    /**
     * Opens the output where to dump the calls.
     *
     * @throws UserException.CouldNotCreateOutputFile if something went wrong creating/opening or writing in the file.
     * @throws UserException.BadInput if there inputs (coverage or targets) are somehow invalid.
     * @throws IllegalArgumentException if any of the parameters is {@code null}.
     */
    protected abstract void openOutput(final File outputFile, final CopyNumberTriStateHiddenMarkovModel model,
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
    protected abstract void makeCalls(final CopyNumberTriStateHiddenMarkovModel model,
                                        final TargetCollection<Target> targets,
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

    private ReadCountCollection readInputCounts(final TargetCollection<Target> targets) {
        try {
            return ReadCountCollectionUtils.parse(inputFile, targets, true);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(inputFile, ex);
        }
    }

    private void checkForMissingTargetsInInputCounts(final TargetCollection<Target> targets, final ReadCountCollection inputCounts) {
        if (inputCounts.targets().size() != targets.targetCount()) {
            logger.warn(String.format("The input counts are missing some of the targets required (%d): e.g. %s",
                    inputCounts.targets().size() - targets.targetCount(), targets.targets().stream()
                            .filter(t -> !inputCounts.targets().contains(t))
                            .limit(5)));
        }
    }

    protected final double logSomeProbability(final int firstTarget, final int length, final CopyNumberTriState call,
                                      final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        // trivial case when the segment has length 0.
        if (length == 0) {
            return 0;
        } else {
            final Set<CopyNumberTriState> otherStates = EnumSet.complementOf(EnumSet.of(call));
            final List<Set<CopyNumberTriState>> otherStatesConstraints = Collections.nCopies(length, otherStates);
            final double logOtherStates = fbResult.logConstrainedProbability(firstTarget, otherStatesConstraints);
            return logProbComplement(logOtherStates);
        }
    }

    /**
     * Calculates the probability of a hidden state switch at the beginning of the segment
     * from a different state to the call state.
     */
    protected final double logStartProbability(final int firstTarget,
                                         final CopyNumberTriState call,
                                         final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        if (firstTarget == 0) {
            return fbResult.logProbability(firstTarget, call);
        } else {
            final EnumSet<CopyNumberTriState> onlyCallState = EnumSet.of(call);
            final EnumSet<CopyNumberTriState> allStatesExceptCall = EnumSet.complementOf(onlyCallState);
            return fbResult.logConstrainedProbability(firstTarget - 1, Arrays.asList(allStatesExceptCall, onlyCallState));
        }
    }

    /**
     * Calculates the probability of a hidden state switch at the end of the segment
     * from the call state to a different state.
     */
    protected final double logEndProbability(final int lastTarget,
                                       final CopyNumberTriState call,
                                       final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        if (lastTarget + 1 == fbResult.positions().size()) {
            return fbResult.logProbability(lastTarget, call);
        } else {
            final EnumSet<CopyNumberTriState> onlyCallState = EnumSet.of(call);
            final EnumSet<CopyNumberTriState> allStatesExceptCall = EnumSet.complementOf(onlyCallState);
            return fbResult.logConstrainedProbability(lastTarget, Arrays.asList(onlyCallState, allStatesExceptCall));
        }
    }

    /**
     * Calculates the complement of a log probability.
     *
     * <p>
     *     With complement of {@code x} we mean: {@code log(1-log(x))}.
     * </p>
     * @param x the input log probability.
     * @return {@code log(1-log(x))}
     */
    private double logProbComplement(final double x) {
        return x >= LN_1_M_EXP_THRESHOLD
                ? Math.log(-Math.expm1(x))
                : Math.log1p(-Math.exp(x));
    }

    /**
     * Transform a log scaled probability (x) into the Phred scaled
     * equivalent or its complement (1-x) Phred scaled equivalent.
     * <p>
     *     This method tolerates probabilities slightly larger than 1.0
     *     (> 0.0 in log scale) which may occur occasionally due to
     *     float point calculation rounding.
     * </p>
     * <p>
     *     The value returned is a phred score capped by {@link #MAX_QUAL_SCORE}.
     * </p>
     *
     * @param rawLogProb the probability.
     * @param complement whether to return the direct Phred transformation ({@code false})
     *                    or its complenent ({@code true)}.
     * @return a values between 0 and {@link #MAX_QUAL_SCORE}.
     * @throws GATKException if {@code rawLogProb} is larger than {@link #MAX_LOG_PROB}.
     */
    protected final double logProbToPhredScore(final double rawLogProb, final boolean complement) {
        if (rawLogProb > MAX_LOG_PROB) {
            throw new GATKException(
                    String.format("numerical instability problem: the log-probability is too large: %g > 0.0 (with maximum tolerance %g)", rawLogProb, MAX_LOG_PROB));
        }
        // make sure that log probs are less than 1 in linear scale.
        // there are cases in that they are just over 0.0 due to float point precision.
        final double logProbEqOrLessThan0 = Math.min(0.0, rawLogProb);

        // Accurate way to calculate log(1-exp(a))
        // based on https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
        final double finalLogProb = complement
                ? logProbComplement(logProbEqOrLessThan0)
                : logProbEqOrLessThan0;

        final double absoluteQualScore = QualityUtils.phredScaleLog10ErrorRate(finalLogProb * INV_LN_10);
        final double exactValue = Math.min(MAX_QUAL_SCORE, absoluteQualScore);
        // We round the value to the required precession.
        return roundPhred(exactValue);
    }

    protected static double roundPhred(final double value) {
        return Math.round(value / PHRED_SCORE_PRECISION) * PHRED_SCORE_PRECISION;
    }

    protected static double roundLog10Prob(final double value) {
        return Math.round(value / LOG10_PROB_PRECISION) * LOG10_PROB_PRECISION;
    }
}
