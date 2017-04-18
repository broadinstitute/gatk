package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

/**
 * This class takes an instance of {@link IntegerCopyNumberTransitionMatrix}, the integer copy number
 * transition matrix on a one-base basis, and computes/caches the log of the integral powers of this matrix
 * to obtain the transition probability over longer intervals (and on the same contig).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberTransitionProbabilityCache implements Serializable {

    private static final long serialVersionUID = 4424349483293880309L;

    private final Map<Integer, LogTransitionProbabilityMatrix> cache = new ConcurrentHashMap<>();

    private final int maxCopyNumber;

    private final double[] eigs;
    private final RealMatrix V;
    private final RealMatrix invV;

    private final double[] logStationaryProbabilities;

    public IntegerCopyNumberTransitionProbabilityCache(@Nonnull final IntegerCopyNumberTransitionMatrix transitionMatrixData) {
        Utils.nonNull(transitionMatrixData, "The integer copy number transition matrix data must be non-null");
        maxCopyNumber = transitionMatrixData.getMaxCopyNumber();
        final LogTransitionProbabilityMatrix T = new LogTransitionProbabilityMatrix(transitionMatrixData.getTransitionMatrix());

        /* perform eigen decomposition */
        final EigenDecomposition ed = new EigenDecomposition(T.getMatrix());
        eigs = ed.getRealEigenvalues();
        V = ed.getV();
        invV = new LUDecomposition(ed.getV()).getSolver().getInverse();

        /* calculate the stationary probabilities */
        logStationaryProbabilities = calculateLogStationaryProbabilities(transitionMatrixData);
        // Arrays.stream(logStationaryProbabilities).mapToObj(Double::toString).forEach(System.out::println);
    }

    /**
     * Calculate log stationary probabilities
     *
     * This is done by unpadding the transition matrix if necessary and finding the stationary distribution via
     * LU decomposition
     *
     * @param transitionMatrixData an instance of {@link IntegerCopyNumberTransitionMatrix}
     * @return an array of doubles
     */
    private double[] calculateLogStationaryProbabilities(@Nonnull final IntegerCopyNumberTransitionMatrix transitionMatrixData) {
        final RealMatrix paddedTransitionMatrix = transitionMatrixData.getTransitionMatrix();
        final RealMatrix unpaddedTransitionMatrix;
        if (transitionMatrixData.getPadding() > 0) {
            unpaddedTransitionMatrix = paddedTransitionMatrix.getSubMatrix(
                    0, paddedTransitionMatrix.getRowDimension() - transitionMatrixData.getPadding() - 1,
                    0, paddedTransitionMatrix.getColumnDimension() - transitionMatrixData.getPadding() - 1);
        } else {
            unpaddedTransitionMatrix = paddedTransitionMatrix;
        }

        /* find the stationary distribution */
        final int dim = unpaddedTransitionMatrix.getColumnDimension();
        final RealMatrix coeff = unpaddedTransitionMatrix.subtract(MatrixUtils.createRealIdentityMatrix(dim));
        coeff.setRow(dim - 1, IntStream.range(0, dim).mapToDouble(i -> 1.0).toArray());
        final RealVector vec = new ArrayRealVector(dim);
        vec.setEntry(dim - 1, 1.0);
        final RealVector unpaddedStationaryDistribution = new LUDecomposition(coeff).getSolver().solve(vec);
        final RealVector paddedStationaryDistribution;
        if (transitionMatrixData.getPadding() > 0) {
            paddedStationaryDistribution = unpaddedStationaryDistribution.append(
                    new ArrayRealVector(transitionMatrixData.getPadding()));
        } else {
            paddedStationaryDistribution = unpaddedStationaryDistribution;
        }
        return Arrays.stream(paddedStationaryDistribution.toArray())
                .map(FastMath::abs) /* some negative values may result from round-off errors */
                .map(FastMath::log)
                .toArray();
    }

    /**
     * Calculates the log transition probability matrix for a given {@param distance} between targets
     *
     * @param distance distance between targets
     * @return an instance of {@link LogTransitionProbabilityMatrix}
     */
    private LogTransitionProbabilityMatrix get(final int distance) {
        if (!cache.containsKey(distance)) {
            calculateMatrix(distance);
        }
        return cache.get(distance);
    }

   /**
     * Caches an instance of log transition matrix for a given distance
     *
     * @param distance distance between targets
     */
    public void cacheLogTransitionMatrix(final int distance) {
        if (!cache.containsKey(distance)) {
            calculateMatrix(distance);
        }
    }

    public void clearCache() {
        cache.clear();
    }

    /**
     * Get the log transition probability between different states
     *
     * @param distance distance between targets
     * @param to destination state
     * @param from departure state
     * @return a double value between (-Infinity) and 0 (inclusive)
     */
    public double logTransitionProbability(final int distance, final IntegerCopyNumberState to, final IntegerCopyNumberState from) {
        /* TODO github/gatk-protected issue #853 -- may remove in the future for performance gains */
        Utils.nonNull(to, "The destination state must be non-null");
        Utils.nonNull(from, "The departure state must be non-null");
        return get(distance).get(to, from);
    }

    public double logStationaryProbability(final IntegerCopyNumberState state) {
        Utils.nonNull(state, "The integer copy number state must be non-null");
        return logStationaryProbabilities[state.getCopyNumber()];
    }

    protected RealMatrix getTransitionProbabilityMatrix(final int distance) {
        final RealMatrix result = get(distance).getMatrix().copy();
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int to, final int from, final double logProbability) {
                return Math.exp(logProbability);
            }
        });
        return result;
    }

    public RealVector getStationaryProbabilityVector() {
        return new ArrayRealVector(Arrays.stream(logStationaryProbabilities).map(FastMath::exp).toArray());
    }

    /**
     * Calculate T^d, where T is the per-base transition matrix and d is the distance in bases, and caches the result
     * @param distance distance between targets (in bases)
     */
    private void calculateMatrix(final int distance) {
        final RealMatrix eigsExponentiated = new DiagonalMatrix(Arrays.stream(eigs)
                .map(x -> Math.pow(x, distance)).toArray());

        /* calculate matrix in probability space, then modify in-place with an element-by-element log */
        final RealMatrix matrixInLogSpace = V.multiply(eigsExponentiated.multiply(invV));
        matrixInLogSpace.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int to, final int from, final double probability) {
                /**
                 * Due to round-off errors, some of the entries can be close to -0.0 rather than +0.0,
                 * resulting in (NaN) rather than (-Infinity). The absolutely value fixes the problem.
                 */
                return Math.log(Math.abs(probability));
            }
        });

        final LogTransitionProbabilityMatrix result = new LogTransitionProbabilityMatrix(matrixInLogSpace);
        cache.put(distance, result);
    }

    public int getMaxCopyNumber() {
        return maxCopyNumber;
    }

    /**
     * Wrapper for a {@link RealMatrix} that is indexed by {@link IntegerCopyNumberState}
     *
     * Note: this class is made public for Serialization issues. Otherwise, we get the following exception
     * when running on GCS:
     *
     * java.lang.IllegalAccessException: Class com.twitter.chill.Instantiators$$anonfun$normalJava$1
     * can not access a member of class org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriStateTransitionProbabilityCache
     * $LogTransitionProbabilityMatrix with modifiers "public"
     *
     */
    public static class LogTransitionProbabilityMatrix implements Serializable {

        private static final long serialVersionUID = -4747127132799560991L;

        private final RealMatrix matrix;

        public LogTransitionProbabilityMatrix(final RealMatrix matrix) {
            if (!matrix.isSquare()) {
                throw new IllegalArgumentException("The copy number transition matrix must be square");
            }
            ParamUtils.isPositive(matrix.getColumnDimension(), "The row/col dimension of the transition matrix must be positive");
            this.matrix = matrix;
        }

        /**
         * Get the log transition probability from an integer copy number state to another
         *
         * @param to the destination integer copy number state
         * @param from the departure integer copy number state
         * @throws OutOfRangeException if either the departure or destination integer copy number states are out of bounds
         * @return a double value
         */
        public double get(final IntegerCopyNumberState to, final IntegerCopyNumberState from) {
            return matrix.getEntry(to.getCopyNumber(), from.getCopyNumber());
        }

        public RealMatrix getMatrix() {
            return matrix;
        }
    }
}
