package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.linear.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * This class calculates and returns transition matrices for different genomic distances with caching under the hood.
 * It is thread-safe for both retrieval and updates, and serializable.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class CopyNumberTriStateTransitionProbabilityCache implements Serializable {

    private static final long serialVersionUID = -2870206010566207323L;

    private final Map<Integer, LogTransitionProbabilityMatrix> cache = new ConcurrentHashMap<>();

    private final double[] Lambda;  // Lambda (the matrix) = diag(Lambda (the array))
    private final RealMatrix V;
    private final RealMatrix invV;

    /**
     * Public constructor.
     *
     * @param meanEventSize mean size of CNV events (in units of bases)
     * @param eventStartProbability the per-base probability of making a transition from one copy number state to another
     */
    public CopyNumberTriStateTransitionProbabilityCache(final double meanEventSize, final double eventStartProbability) {
        ParamUtils.isPositive(meanEventSize, "Mean event size must be positive");
        ParamUtils.inRange(eventStartProbability, 0, 1, "Event start probability must be between 0 and 1");

        final double eventContinueProbability = Math.exp(-1 / meanEventSize);
        final LogTransitionProbabilityMatrix T = new LogTransitionProbabilityMatrix();

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

    /**
     * Get the log transition probability between different states
     *
     * @param distance distance between targets
     * @param to destination state
     * @param from departure state
     * @return a double value between (-Infinity) and 0 (inclusive)
     */
    public double logProbability(final int distance, final CopyNumberTriState to, final CopyNumberTriState from) {
        /* TODO github/gatk-protected issue #853 -- may remove in the future for performance gains */
        Utils.nonNull(to, "The destination state must be non-null");
        Utils.nonNull(from, "The departure state must be non-null");
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

    /**
     * Calculate T^d, where T is the per-base transition matrix and d is the distance in bases, and caches the result
     * @param distance distance between targets (in bases)
     */
    private void calculateMatrix(final int distance) {
        final RealMatrix LambdaExponentiated = new DiagonalMatrix(Arrays.stream(Lambda)
                .map(x -> Math.pow(x, distance)).toArray());

        // calculate matrix in probability space, then modify in-place with an element-by-element log
        final RealMatrix matrixInLogSpace = V.multiply(LambdaExponentiated.multiply(invV));
        matrixInLogSpace.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int to, final int from, final double probability) {
                /**
                 * Due to round-off errors, some of the entries can be close to -0.0 rather than +0.0,
                 * resulting in (NaN) rather than (-Infinity). The absolute value fixes the problem.
                 */
                return Math.log(Math.abs(probability));
            }
        });

        final LogTransitionProbabilityMatrix result = new LogTransitionProbabilityMatrix(matrixInLogSpace);
        cache.put(distance, result);
    }

    /**
     * Wrapper for a {@link RealMatrix} that is indexed by the CopyNumberTriState enum
     *
     * @implNote This class is made public for Serialization issues. Otherwise, we get the following exception
     * when running on GCS:
     *
     * java.lang.IllegalAccessException: Class com.twitter.chill.Instantiators$$anonfun$normalJava$1
     * can not access a member of class org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriStateTransitionProbabilityCache
     * $LogTransitionProbabilityMatrix with modifiers "public"
     *
     */
    public static class LogTransitionProbabilityMatrix implements Serializable {

        private static final long serialVersionUID = -9072224697693405187L;

        final static int size = CopyNumberTriState.values().length;
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