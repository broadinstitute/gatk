package org.broadinstitute.hellbender.tools.exome;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Doubles;
import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.linalg.DenseMatrix;
import org.apache.spark.mllib.linalg.Matrix;
import org.apache.spark.mllib.linalg.distributed.RowMatrix;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hdf5.PoN;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.spark.SparkConverter;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Utility class for performing tangent normalization (and related operations).
 *
 * Currently, only supports tangent normalization in the reduced hyperplane, not the logNormal hyperplane
 */
public final class TangentNormalizer {
    private static final Logger logger = LogManager.getLogger(TangentNormalizer.class);

    /**
     * Minimum target normalized and column centered count possible.
     *
     * <p>
     *     It must be small yet greater than 0 to avoid -Inf problems in the calculations.
     * </p>
     */
    public static double EPSILON = Math.pow(10, -10);
    public static double LOG_2_EPSILON = Math.log(EPSILON) / Math.log(2.0);

    /**
     * Cached inverse of the natural logarithm of 2.
     */
    private static double INV_LN2 = 1.0 / Math.log(2.0);

    private static final int TN_NUM_SLICES_SPARK = 50;

    private TangentNormalizer() { }

    /**
     * Tangent normalize a coverage profile.
     *
     * <p>Notes about the Spark tangent normalization can be found in docs/PoN/</p>
     *
     * @param pon Not {@code null}
     * @param targetFactorNormalizedCounts ReadCountCollection of counts that have already been normalized fully
     *                                     (typically, including the target factor normalization).  I.e. a coverage profile
     *                                     The column names should be intact.  Not {@code null}
     *                                     See {@link TangentNormalizer::createCoverageProfile}
     * @return never {@code null}
     */
    private static TangentNormalizationResult tangentNormalize(final PoN pon, final ReadCountCollection targetFactorNormalizedCounts, JavaSparkContext ctx) {

        Utils.nonNull(pon, "PoN cannot be null.");
        Utils.nonNull(targetFactorNormalizedCounts, "targetFactorNormalizedCounts cannot be null.");
        Utils.nonNull(targetFactorNormalizedCounts.columnNames(), "targetFactorNormalizedCounts column names cannot be null.");
        ParamUtils.isPositive(targetFactorNormalizedCounts.columnNames().size(), "targetFactorNormalizedCounts column names cannot be an empty list.");

        final Case2PoNTargetMapper targetMapper = new Case2PoNTargetMapper(targetFactorNormalizedCounts.targets(), pon.getPanelTargetNames());

        // The input counts with rows (targets) sorted so that they match the PoN's order.
        final RealMatrix tangentNormalizationRawInputCounts = targetMapper.fromCaseToPoNCounts(targetFactorNormalizedCounts.counts());

        // We prepare the counts for tangent normalization.
        final RealMatrix tangentNormalizationInputCounts = composeTangentNormalizationInputMatrix(tangentNormalizationRawInputCounts);

        if (ctx == null) {
            // Calculate the beta-hats for the input read count columns (samples).
            logger.info("Calculating beta hats...");
            final RealMatrix tangentBetaHats = pon.betaHats(tangentNormalizationInputCounts, true, EPSILON);

            // Actual tangent normalization step.
            logger.info("Performing actual tangent normalization (" + tangentNormalizationInputCounts.getColumnDimension() + " columns)...");
            final RealMatrix tangentNormalizedCounts = pon.tangentNormalization(tangentNormalizationInputCounts, tangentBetaHats, true);

            // Output the tangent normalized counts.
            logger.info("Post-processing tangent normalization results...");
            final ReadCountCollection tangentNormalized = targetMapper.fromPoNtoCaseCountCollection(tangentNormalizedCounts, targetFactorNormalizedCounts.columnNames());
            final ReadCountCollection preTangentNormalized = targetMapper.fromPoNtoCaseCountCollection(tangentNormalizationInputCounts, targetFactorNormalizedCounts.columnNames());

            return new TangentNormalizationResult(tangentNormalized, preTangentNormalized, tangentBetaHats, targetFactorNormalizedCounts);
        } else {

            /*
            Using Spark:  the code here is a little more complex for optimization purposes.

            Please see notes in docs/PoN ...

            Ahat^T = (C^T P^T) A^T
            Therefore, C^T is the RowMatrix

            pinv: P
            panel: A
            projection: Ahat
            cases: C
            betahat: C^T P^T
            tangentNormalizedCounts: C - Ahat
             */
            final RealMatrix pinv = pon.getReducedPanelPInverseCounts();
            final RealMatrix panel = pon.getReducedPanelCounts();

            // Make the C^T a distributed matrix (RowMatrix)
            final RowMatrix caseTDistMat = SparkConverter.convertRealMatrixToSparkRowMatrix(ctx, tangentNormalizationInputCounts.transpose(), TN_NUM_SLICES_SPARK);

            // Spark local matrices (transposed)
            final Matrix pinvTLocalMat = new DenseMatrix(pinv.getRowDimension(), pinv.getColumnDimension(), Doubles.concat(pinv.getData()), true).transpose();
            final Matrix panelTLocalMat = new DenseMatrix(panel.getRowDimension(), panel.getColumnDimension(), Doubles.concat(panel.getData()), true).transpose();

            // Calculate the projection transpose in a distributed matrix, then convert to Apache Commons matrix (not transposed)
            final RowMatrix betahatDistMat = caseTDistMat.multiply(pinvTLocalMat);
            final RowMatrix projectionTDistMat = betahatDistMat.multiply(panelTLocalMat);
            final RealMatrix projection = SparkConverter.convertSparkRowMatrixToRealMatrix(projectionTDistMat, tangentNormalizationInputCounts.transpose().getRowDimension()).transpose();

            // Subtract the cases from the projection
            final RealMatrix tangentNormalizedCounts = tangentNormalizationInputCounts.subtract(projection);

            // Construct the result object and return it with the correct targets.
            final ReadCountCollection tangentNormalized = targetMapper.fromPoNtoCaseCountCollection(tangentNormalizedCounts, targetFactorNormalizedCounts.columnNames());
            final ReadCountCollection preTangentNormalized = targetMapper.fromPoNtoCaseCountCollection(tangentNormalizationInputCounts, targetFactorNormalizedCounts.columnNames());
            final RealMatrix tangentBetaHats = SparkConverter.convertSparkRowMatrixToRealMatrix(betahatDistMat, tangentNormalizedCounts.getColumnDimension());
            return new TangentNormalizationResult(tangentNormalized, preTangentNormalized, tangentBetaHats.transpose(), targetFactorNormalizedCounts);
        }
    }

    /**
     *  Do the full tangent normalization process given proportional coverage data.
     *
     *  This includes:
     *   <ul><li>normalization by target factors</li>
     *   <li>projection of the normalized coverage profile into the hyperplane from the PoN</li>
     *   </ul>
     *
     * @param pon -- never {@code null}
     * @param pcov -- never {@code null}.  Must contain data for at least one sample.
     * @param ctx spark context.  Use {@code null} if no context is available
     * @return never {@code null}
     */
    public static TangentNormalizationResult tangentNormalizePcov(final PoN pon, final ReadCountCollection pcov, final JavaSparkContext ctx) {
        Utils.nonNull(pon, "PoN cannot be null.");
        Utils.nonNull(pcov, "input pcov read counts cannot be null when creating a coverage profile.");
        ParamUtils.isPositive(pcov.columnNames().size(), "input cov profile column names cannot be an empty list.");
        final ReadCountCollection coverageProfile = createCoverageProfile(pon, pcov);
        return TangentNormalizer.tangentNormalize(pon, coverageProfile, ctx);
    }

    /**
     * {@link TangentNormalizer :: tangentNormalizePcov} with no spark context
     *
     * @param pon never {@code null}.  Must contain data for at least one sample.
     * @param pcov never {@code null}.  Must contain data for at least one sample.
     * @return never {@code null}
     */
    public static TangentNormalizationResult tangentNormalizePcov(final PoN pon, final ReadCountCollection pcov) {
        return TangentNormalizer.tangentNormalizePcov(pon, pcov, null);
    }

    private static ReadCountCollection createCoverageProfile(final PoN pon, final ReadCountCollection inputReadCounts) {
        Utils.nonNull(pon, "PoN cannot be null.");
        Utils.nonNull(inputReadCounts, "input read counts cannot be null when creating a coverage profile.");
        ParamUtils.isPositive(inputReadCounts.columnNames().size(), "inputReadCounts column names cannot be an empty list.");
        final Case2PoNTargetMapper targetMapper = new Case2PoNTargetMapper(inputReadCounts.targets(), pon.getTargetNames());
        final RealMatrix inputCounts = targetMapper.fromCaseToPoNCounts(inputReadCounts.counts());
        final RealMatrix targetNormalizedCounts = pon.factorNormalization(inputCounts);

        return targetMapper.fromPoNtoCaseCountCollection(targetNormalizedCounts, inputReadCounts.columnNames());
    }

    /**
     * Project all of the normals used to create the PoN into the reduced panel.
     *
     * Same as {@link TangentNormalizer :: tangentNormalizeNormalsInPoN} but will not attempt to use spark context.
     * @param pon Not {@code null}.
     * @return never {@code null}.  Result will contain multiple columns in each of the ReadCountCollection attributes.
     */
    public static TangentNormalizationResult tangentNormalizeNormalsInPoN(final PoN pon) {
        return tangentNormalizeNormalsInPoN(pon, null);
    }

    /**
     * Project all of the normals used to create the PoN into the reduced panel.
     *
     * @param pon Not {@code null}.
     * @param ctx Spark context,  {@code null} indicates no spark context available, which is supported.
     * @return never {@code null}.  Result will contain multiple columns in each of the ReadCountCollection attributes.
     */
    public static TangentNormalizationResult tangentNormalizeNormalsInPoN(final PoN pon, final JavaSparkContext ctx) {
        // Get the list of sample names in a modifiable List
        final List<String> sampleNamesCopy = new ArrayList<>(pon.getSampleNames());

        logger.info("Loading normalized counts...");
        final ReadCountCollection coverageProfile = new ReadCountCollection(SetUniqueList.setUniqueList(pon.getTargets()), SetUniqueList.setUniqueList(sampleNamesCopy), pon.getNormalizedCounts());

        logger.info("Tangent normalizing the normals (normalized by target factors) ...");

        // For each sample in the PoN, tangent normalize against the qc reduced PoN.
        return TangentNormalizer.tangentNormalize(pon, coverageProfile, ctx);
    }

    /**
     * Prepares the data to perform tangent factor normalization.
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
    @VisibleForTesting
    static RealMatrix composeTangentNormalizationInputMatrix(final RealMatrix matrix) {
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
