package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.util.Arrays;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class IntegerCopyNumberTransitionProbabilityCacheUnitTest extends BaseTest {

    private static final double EPSILON = 1e-6;

    /**
     * Some distances
     */
    private static final int[] DISTANCES = {1, 2, 5, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000,
        1000000000};

    private static final int[] PADDINGS = {0, 1, 2, 5, 10};

    /**
     * Some normalized transition matrices of different sizes
     */
    private static final RealMatrix[] TRANSITION_MATRICES = {
            new Array2DRowRealMatrix(new double[][] {
                    {0.1, 0.5, 0.9},
                    {0.6, 0.2, 0.05},
                    {0.3, 0.3, 0.05}}),
            new Array2DRowRealMatrix(new double[][] {
                    {0.3, 0.5},
                    {0.7, 0.5}}),
            new Array2DRowRealMatrix(new double[][] {
                    {0.1, 0.5, 0.9, 0.1},
                    {0.3, 0.2, 0.05, 0.2},
                    {0.6, 0.1, 0.05, 0.1},
                    {0.0, 0.2, 0.0, 0.6}})
    };

    @Test
    public void testBasicSoundness() {
        for (final RealMatrix transitionMatrix : TRANSITION_MATRICES) {
            final IntegerCopyNumberTransitionProbabilityCache cache = new IntegerCopyNumberTransitionProbabilityCache(
                    new IntegerCopyNumberTransitionMatrixData(transitionMatrix, 0));
            for (final int dist : DISTANCES) {
                final RealMatrix transitionMatrixExponentiated = cache.getTransitionProbabilityMatrix(dist);

                /* assert positivity */
                Assert.assertTrue(Arrays.stream(transitionMatrixExponentiated.getData())
                        .flatMapToDouble(Arrays::stream)
                        .allMatch(d -> d >= 0));

                /* assert conservation of probability */
                for (int c = 0; c < transitionMatrix.getColumnDimension(); c++) {
                    Assert.assertEquals(Arrays.stream(transitionMatrixExponentiated.getColumn(c)).sum(), 1.0, EPSILON);
                }

                /* assert correctness, T(2*d) = T(d)*T(d) */
                assertEqualMatrices(cache.getTransitionProbabilityMatrix(2*dist),
                        transitionMatrixExponentiated.multiply(transitionMatrixExponentiated));
            }

            /* assert loss of initial state over long distances, i.e. all columns must be equal */
            final RealMatrix longRangeTransitionMatrix = cache.getTransitionProbabilityMatrix(Integer.MAX_VALUE);
            final double[] firstColumn = longRangeTransitionMatrix.getColumn(0);
            final RealMatrix syntheticLongRangeTransitionMatrix = new Array2DRowRealMatrix(firstColumn.length,
                    firstColumn.length);
            for (int i = 0; i < firstColumn.length; i++) {
                syntheticLongRangeTransitionMatrix.setColumn(i, firstColumn);
            }
            assertEqualMatrices(longRangeTransitionMatrix, syntheticLongRangeTransitionMatrix);

            final double[] stationary = cache.getStationaryProbabilityVector().toArray();
            ArrayAsserts.assertArrayEquals(stationary, firstColumn, EPSILON);
        }
    }

    @Test
    public void testStationaryProbabilities() {
        for (final int padding : PADDINGS) {
            for (final RealMatrix transitionMatrix : TRANSITION_MATRICES) {
                final IntegerCopyNumberTransitionProbabilityCache cache = new IntegerCopyNumberTransitionProbabilityCache(
                        new IntegerCopyNumberTransitionMatrixData(transitionMatrix, padding));
                final double[] stationary = cache.getStationaryProbabilityVector().toArray();
                final double[] stationaryFromMultiplication = cache.getTransitionProbabilityMatrix(Integer.MAX_VALUE)
                        .getColumn(0);
                ArrayAsserts.assertArrayEquals(stationary, stationaryFromMultiplication, EPSILON);
                Assert.assertEquals(Arrays.stream(stationary).sum(), 1.0, EPSILON);
            }
        }

    }

    private static void assertEqualMatrices(final RealMatrix mat1, final RealMatrix mat2) {
        Assert.assertEquals(mat1.subtract(mat2).getNorm(), 0, EPSILON);
    }

}
