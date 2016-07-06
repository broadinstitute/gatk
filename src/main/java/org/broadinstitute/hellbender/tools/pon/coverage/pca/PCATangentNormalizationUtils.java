package org.broadinstitute.hellbender.tools.pon.coverage.pca;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.linalg.DenseMatrix;
import org.apache.spark.mllib.linalg.Matrix;
import org.apache.spark.mllib.linalg.distributed.RowMatrix;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.pon.coverage.CaseToPoNTargetMapper;
import org.broadinstitute.hellbender.tools.pon.coverage.CoveragePanelOfNormals;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.spark.SparkConverter;

import java.util.List;
import java.util.stream.IntStream;

/**
 * Utility class for package-private methods for performing tangent normalization (and related operations).
 *
 * Currently, only supports tangent normalization in the reduced hyperplane, not the logNormal hyperplane
 */
public final class PCATangentNormalizationUtils {
    private static final Logger logger = LogManager.getLogger(PCATangentNormalizationUtils.class);

    /**
     * Minimum target normalized and column centered count possible.
     *
     * <p>
     *     It must be small yet greater than 0 to avoid -Inf problems in the calculations.
     * </p>
     */
    public static final double EPSILON = 1E-9;
    private static final double INV_LN2 = GATKProtectedMathUtils.INV_LN2;
    private static final double LOG_2_EPSILON = Math.log(EPSILON) * INV_LN2;

    private static final int TN_NUM_SLICES_SPARK = 50;

    private PCATangentNormalizationUtils() {}

    /**
     * Target-factor normalizes a {@link RealMatrix} in-place given target factors..
     */
    static void factorNormalize(final RealMatrix input, final double[] targetFactors) {
        Utils.nonNull(input, "Input matrix cannot be null.");
        Utils.nonNull(targetFactors, "Target factors cannot be null.");
        Utils.validateArg(targetFactors.length == input.getRowDimension(),
                "Number of target factors does not correspond to the number of rows.");
        // Divide all counts by the target factor for the row.
        input.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return value / targetFactors[row];
            }
        });
    }

    /**
     *  Do the full tangent normalization process given a {@link PCACoveragePoN} and a proportional-coverage profile.
     *
     *  This includes:
     *   <ul><li>normalization by target factors (optional)</li>
     *   <li>projection of the normalized coverage profile into the hyperplane from the PoN</li>
     *   </ul>
     *
     * @param pon -- never {@code null}
     * @param profile -- never {@code null}.  Must contain data for at least one sample.
     * @param ctx spark context.  Use {@code null} if no context is available
     * @param doFactorNormalization if true, perform factor normalization (set to false for normalizing normals in a PoN that have already been factor normalized)
     * @return never {@code null}
     */
    static PCATangentNormalizationResult tangentNormalize(final PCACoveragePoN pon,
                                                          final ReadCountCollection profile,
                                                          final boolean doFactorNormalization,
                                                          final JavaSparkContext ctx) {
        Utils.nonNull(pon, "PoN cannot be null.");
        Utils.nonNull(profile, "Proportional coverages cannot be null.");
        ParamUtils.isPositive(profile.columnNames().size(), "Column names cannot be an empty list.");

        //normals stored in a PoN may already be factor normalized
        final ReadCountCollection factorNormalizedCoverage = doFactorNormalization ? mapTargetsToPoNAndFactorNormalize(profile, pon) : profile;

        return tangentNormalize(factorNormalizedCoverage, pon.getPanelTargetNames(), pon.getReducedPanelCounts(), pon.getReducedPanelPInverseCounts(), ctx);
    }

    /**
     * Project all of the normals used to create the PoN into the reduced panel using the raw PoN data.
     * This is required to calculate target variances for PoN initialization
     * (since target variances are considered part of {@link PCACoveragePoN}, this cannot be done simply using
     * {@link CoveragePanelOfNormals#normalizeNormalsInPoN} since the PoN is not fully initialized).
     */
    static PCATangentNormalizationResult tangentNormalizeNormalsInPoN(final ReadCountCollection normalizedCounts,
                                                                      final List<String> panelTargetNames,
                                                                      final RealMatrix reducedPanelCounts,
                                                                      final RealMatrix reducedPanelPInvCounts,
                                                                      final JavaSparkContext ctx) {
        Utils.nonNull(normalizedCounts, "Normalized counts cannot be null.");
        Utils.nonNull(panelTargetNames, "Panel target names cannot be null.");
        Utils.nonNull(reducedPanelCounts, "Reduced-panel counts cannot be null.");
        Utils.nonNull(reducedPanelPInvCounts, "Reduced-panel pseudoinverse cannot be null.");
        // For each sample in the PoN, tangent normalize against the qc reduced PoN.
        logger.info("Tangent normalizing the normals (normalized by target factors) ...");
        return tangentNormalize(normalizedCounts, panelTargetNames, reducedPanelCounts, reducedPanelPInvCounts, ctx);
    }

    /**
     * Returns a target-factor-normalized {@link ReadCountCollection} given a {@link PCACoveragePoN}..
     */
    private static ReadCountCollection mapTargetsToPoNAndFactorNormalize(final ReadCountCollection input, final PCACoveragePoN pon) {
        final CaseToPoNTargetMapper targetMapper = new CaseToPoNTargetMapper(input.targets(), pon.getTargetNames());
        final RealMatrix inputCounts = targetMapper.fromCaseToPoNCounts(input.counts());
        factorNormalize(inputCounts, pon.getTargetFactors());   //factor normalize in-place
        return targetMapper.fromPoNtoCaseCountCollection(inputCounts, input.columnNames());
    }

    /**
     * Tangent normalize given the raw PoN data.  Non-Spark or Spark implementation automatically chosen.
     */
    private static PCATangentNormalizationResult tangentNormalize(final ReadCountCollection targetFactorNormalizedCounts,
                                                                  final List<String> panelTargetNames,
                                                                  final RealMatrix reducedPanelCounts,
                                                                  final RealMatrix reducedPanelPInvCounts,
                                                                  final JavaSparkContext ctx) {
        final CaseToPoNTargetMapper targetMapper = new CaseToPoNTargetMapper(targetFactorNormalizedCounts.targets(), panelTargetNames);

        // The input counts with rows (targets) sorted so that they match the PoN's order.
        final RealMatrix tangentNormalizationRawInputCounts = targetMapper.fromCaseToPoNCounts(targetFactorNormalizedCounts.counts());

        // We prepare the counts for tangent normalization.
        final RealMatrix tangentNormalizationInputCounts = composeTangentNormalizationInputMatrix(tangentNormalizationRawInputCounts);

        if (ctx == null) {
            return tangentNormalizeNonSpark(targetFactorNormalizedCounts, reducedPanelCounts, reducedPanelPInvCounts, targetMapper, tangentNormalizationInputCounts);
        } else {
            return tangentNormalizeSpark(targetFactorNormalizedCounts, reducedPanelCounts, reducedPanelPInvCounts, targetMapper, tangentNormalizationInputCounts, ctx);
        }
    }

    /**
     * Tangent normalize given the raw PoN data without using Spark.
     */
    private static PCATangentNormalizationResult tangentNormalizeNonSpark(final ReadCountCollection targetFactorNormalizedCounts,
                                                                          final RealMatrix reducedPanelCounts,
                                                                          final RealMatrix reducedPanelPInvCounts,
                                                                          final CaseToPoNTargetMapper targetMapper,
                                                                          final RealMatrix tangentNormalizationInputCounts) {
        // Calculate the beta-hats for the input read count columns (samples).
        logger.info("Calculating beta hats...");
        final RealMatrix tangentBetaHats = calculateBetaHats(reducedPanelPInvCounts, tangentNormalizationInputCounts, EPSILON);

        // Actual tangent normalization step.
        logger.info("Performing actual tangent normalization (" + tangentNormalizationInputCounts.getColumnDimension() + " columns)...");
        final RealMatrix tangentNormalizedCounts = tangentNormalize(reducedPanelCounts, tangentNormalizationInputCounts, tangentBetaHats);

        // Output the tangent normalized counts.
        logger.info("Post-processing tangent normalization results...");
        final ReadCountCollection tangentNormalized = targetMapper.fromPoNtoCaseCountCollection(
                tangentNormalizedCounts, targetFactorNormalizedCounts.columnNames());
        final ReadCountCollection preTangentNormalized = targetMapper.fromPoNtoCaseCountCollection(
                tangentNormalizationInputCounts, targetFactorNormalizedCounts.columnNames());

        return new PCATangentNormalizationResult(tangentNormalized, preTangentNormalized, tangentBetaHats, targetFactorNormalizedCounts);
    }

    /**
     * Calculate the beta-hats that best fit case read counts given the panel of normals.
     *
     * @param normalsPseudoinverse the log-normalized or reduced-panel pseudoinverse from a panel of normals
     * @param input a {@code TxS} matrix where {@code T} is the number of targets and {@code S} the number of count groups (e.g. case samples).
     * @return never {@code null} an {@code NxS} matrix, where N is the number of samples in
     *  the panel and S the original name of count groups.
     */
    @VisibleForTesting
    public static RealMatrix calculateBetaHats(final RealMatrix normalsPseudoinverse,
                                               final RealMatrix input,
                                               final double epsilon) {
        Utils.nonNull(normalsPseudoinverse, "Normals inverse matrix cannot be null.");
        Utils.nonNull(input, "Input counts cannot be null.");
        Utils.validateArg(epsilon > 0, String.format("Invalid epsilon value, must be > 0: %f", epsilon));
        final double targetThreshold = (Math.log(epsilon) / Math.log(2)) + 1;

        // copy case samples in order to mask targets in-place and mask (set to zero) targets with coverage below threshold
        final RealMatrix maskedInput = input.copy();
        maskedInput.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return value > targetThreshold ? value : 0;
            }
        });

        return normalsPseudoinverse.multiply(maskedInput);
    }

    /**
     * Applies tangent normalization.
     * <p>
     *     The input row order should match the panel's target order.
     * </p>
     *
     * @param normals the log-normalized or reduced-panel counts from a panel of normals
     * @param input the input counts to normalize. This matrix is TxS where T is the number of targets
     *              and S the number of count columns.
     * @param betaHats the beta-hats for the projection to use for the normalization. This matrix
     *                  is NxS where N is the number of samples in the panel of choice and S is the number of count columns.
     * @return never {@code null}.
     */
    private static RealMatrix tangentNormalize(final RealMatrix normals,
                                               final RealMatrix input,
                                               final RealMatrix betaHats) {
        Utils.validateArg(input.getColumnDimension() == betaHats.getColumnDimension(),
                String.format("the input count column count (%d) does not match the number of columns in the beta-hats (%d)",
                        input.getColumnDimension(), betaHats.getColumnDimension()));
        Utils.validateArg(normals.getColumnDimension() == betaHats.getRowDimension(),
                String.format("beta-hats component count (%d) does not match the number of samples in the PoN (%d)",
                        normals.getRowDimension(), normals.getColumnDimension()));
        final RealMatrix projection = normals.multiply(betaHats);
        return input.subtract(projection);
    }

    /**
     * Tangent normalize given the raw PoN data using Spark:  the code here is a little more complex for optimization purposes.
     *
     *  Please see notes in docs/PoN ...
     *
     *  Ahat^T = (C^T P^T) A^T
     *  Therefore, C^T is the RowMatrix
     *
     *  pinv: P
     *  panel: A
     *  projection: Ahat
     *  cases: C
     *  betahat: C^T P^T
     *  tangentNormalizedCounts: C - Ahat
     */
    private static PCATangentNormalizationResult tangentNormalizeSpark(final ReadCountCollection targetFactorNormalizedCounts,
                                                                       final RealMatrix reducedPanelCounts,
                                                                       final RealMatrix reducedPanelPInvCounts,
                                                                       final CaseToPoNTargetMapper targetMapper,
                                                                       final RealMatrix tangentNormalizationInputCounts,
                                                                       final JavaSparkContext ctx) {
        // Make the C^T a distributed matrix (RowMatrix)
        final RowMatrix caseTDistMat = SparkConverter.convertRealMatrixToSparkRowMatrix(
                ctx, tangentNormalizationInputCounts.transpose(), TN_NUM_SLICES_SPARK);

        // Spark local matrices (transposed)
        final Matrix pinvTLocalMat = new DenseMatrix(
                reducedPanelPInvCounts.getRowDimension(), reducedPanelPInvCounts.getColumnDimension(),
                Doubles.concat(reducedPanelPInvCounts.getData()), true).transpose();
        final Matrix panelTLocalMat = new DenseMatrix(
                reducedPanelCounts.getRowDimension(), reducedPanelCounts.getColumnDimension(),
                Doubles.concat(reducedPanelCounts.getData()), true).transpose();

        // Calculate the projection transpose in a distributed matrix, then convert to Apache Commons matrix (not transposed)
        final RowMatrix betahatDistMat = caseTDistMat.multiply(pinvTLocalMat);
        final RowMatrix projectionTDistMat = betahatDistMat.multiply(panelTLocalMat);
        final RealMatrix projection = SparkConverter.convertSparkRowMatrixToRealMatrix(
                projectionTDistMat, tangentNormalizationInputCounts.transpose().getRowDimension()).transpose();

        // Subtract the projection from the cases
        final RealMatrix tangentNormalizedCounts = tangentNormalizationInputCounts.subtract(projection);

        // Construct the result object and return it with the correct targets.
        final ReadCountCollection tangentNormalized = targetMapper.fromPoNtoCaseCountCollection(
                tangentNormalizedCounts, targetFactorNormalizedCounts.columnNames());
        final ReadCountCollection preTangentNormalized = targetMapper.fromPoNtoCaseCountCollection(
                tangentNormalizationInputCounts, targetFactorNormalizedCounts.columnNames());
        final RealMatrix tangentBetaHats = SparkConverter.convertSparkRowMatrixToRealMatrix(
                betahatDistMat, tangentNormalizedCounts.getColumnDimension());

        return new PCATangentNormalizationResult(tangentNormalized, preTangentNormalized, tangentBetaHats.transpose(), targetFactorNormalizedCounts);
    }

    /**
     * Prepares the data to perform tangent normalization.
     * <p>
     * This is done by count group or column:
     *   <ol>
     *     </li>we divide counts by the column mean,</li>
     *     </li>then we transform value to their log_2,</li>
     *     </li>and finally we center them around the median.</li>
     *   </ol>
     * </p>
     *
     * @param matrix input matrix.
     * @return never {@code null}.
     */
    private static RealMatrix composeTangentNormalizationInputMatrix(final RealMatrix matrix) {
        final RealMatrix result = matrix.copy();

        // step 1: divide by column means and log_2 transform
        final double[] columnMeans = GATKProtectedMathUtils.columnMeans(matrix);
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return truncatedLog2(value / columnMeans[column]);
            }
        });

        // step 2: subtract column medians
        final double[] columnMedians = IntStream.range(0, matrix.getColumnDimension())
                .mapToDouble(c -> new Median().evaluate(result.getColumn(c))).toArray();
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return value - columnMedians[column];
            }
        });

        return result;
    }

    private static double truncatedLog2(final double x) {
        return x < EPSILON ? LOG_2_EPSILON : Math.log(x) * INV_LN2;
    }
}
