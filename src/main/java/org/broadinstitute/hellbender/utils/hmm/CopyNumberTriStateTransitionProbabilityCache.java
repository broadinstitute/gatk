package org.broadinstitute.hellbender.utils.hmm;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.linear.*;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * store transition matrices for different genomic distances
 *
 * Created by davidben on 3/10/16.
 */
public class CopyNumberTriStateTransitionProbabilityCache {

    private final Map<Integer, LogTransitionProbabilityMatrix> cache = new HashMap<>();
    private final double eventContinueProbability;

    // diagonalize the per-base transition matrix T for fast exponentiation: T = V*Lambda*inv(V), hence T^d = V*Lambda^d*inv(V)
    private final LogTransitionProbabilityMatrix T;
    private final double[] Lambda;  // Lambda (the matrix) = diag(Lambda (the array))
    private final RealMatrix V;
    private final RealMatrix invV;

    public CopyNumberTriStateTransitionProbabilityCache(final double meanEventSize, final double eventStartProbability) {
        ParamUtils.isPositive(meanEventSize, "Mean event size must be positive");
        ParamUtils.inRange(eventStartProbability, 0, 1, "Event start probabiloity must be between 0 and 1");
        eventContinueProbability = Math.exp(-1 / meanEventSize);
        T = new LogTransitionProbabilityMatrix();
        T.set(CopyNumberTriState.DELETION, CopyNumberTriState.DELETION, eventContinueProbability);
        T.set(CopyNumberTriState.DUPLICATION, CopyNumberTriState.DUPLICATION, eventContinueProbability);

        T.set(CopyNumberTriState.DELETION, CopyNumberTriState.NEUTRAL, eventStartProbability);
        T.set(CopyNumberTriState.DUPLICATION, CopyNumberTriState.NEUTRAL, eventStartProbability);
        T.set(CopyNumberTriState.NEUTRAL, CopyNumberTriState.NEUTRAL, 1 - 2 * eventStartProbability);

        T.set(CopyNumberTriState.DELETION, CopyNumberTriState.DUPLICATION, 0);
        T.set(CopyNumberTriState.DUPLICATION, CopyNumberTriState.DELETION, 0);

        T.set(CopyNumberTriState.NEUTRAL, CopyNumberTriState.DELETION, 1 - eventContinueProbability);
        T.set(CopyNumberTriState.NEUTRAL, CopyNumberTriState.DUPLICATION, 1 - eventContinueProbability);


        final EigenDecomposition ed = new EigenDecomposition(T.getMatrix());
        Lambda = ed.getRealEigenvalues();
        V = ed.getV();
        invV = new LUDecomposition(ed.getV()).getSolver().getInverse();
    }

    private LogTransitionProbabilityMatrix get(final int distance) {
        if (!cache.containsKey(distance)) {
            calculateMatrix(distance);
        }
        return cache.get(distance);
    }

    public double logProbability(final int distance, final CopyNumberTriState to, final CopyNumberTriState from) {
        return get(distance).get(to, from);
    }

    @VisibleForTesting
    RealMatrix getAsMatrixInProbabilitySpace(final int distance) {
        final RealMatrix result = get(distance).getMatrix().copy();
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int to, final int from, final double logProbability) {
                return Math.exp(logProbability);
            }
        });
        return result;
    }

    // calculate T^d, where T is the per-base transition matrix and d is the distance in bases
    private void calculateMatrix(final int distance) {
        final RealMatrix LambdaExponentiated = new DiagonalMatrix(Arrays.stream(Lambda).map(x -> Math.pow(x, distance)).toArray());

        //calculate matrix in probability space, then modify in-place with an element-by-element log
        final RealMatrix matrixInLogSpace = V.multiply(LambdaExponentiated.multiply(invV));

        matrixInLogSpace.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int to, final int from, final double probability) {
                return Math.log(probability);
            }
        });

        final LogTransitionProbabilityMatrix result = new LogTransitionProbabilityMatrix(matrixInLogSpace);
        cache.put(distance, result);
    }

    // wrapper for a real matrix that is indexed by the CopyNumberTriState enum
    private static class LogTransitionProbabilityMatrix {
        static final int size = CopyNumberTriState.values().length;
        private final RealMatrix matrix;

        public LogTransitionProbabilityMatrix() {
            matrix = new Array2DRowRealMatrix(size, size);
        }

        public LogTransitionProbabilityMatrix(final RealMatrix matrix) {
            this.matrix = matrix;
        }

        public void set(final CopyNumberTriState to, final CopyNumberTriState from, final double value) {
            matrix.setEntry(to.ordinal(), from.ordinal(), value);
        }

        public double get(final CopyNumberTriState to, final CopyNumberTriState from) {
            return matrix.getEntry(to.ordinal(), from.ordinal());
        }

        public RealMatrix getMatrix() {
            return matrix;
        }
    }
}
