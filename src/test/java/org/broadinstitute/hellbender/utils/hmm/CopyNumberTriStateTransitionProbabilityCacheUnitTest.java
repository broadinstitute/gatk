package org.broadinstitute.hellbender.utils.hmm;

import org.apache.commons.math3.linear.DefaultRealMatrixPreservingVisitor;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Created by davidben on 3/10/16.
 */
public class CopyNumberTriStateTransitionProbabilityCacheUnitTest {
    private static final int[] DISTANCES = new int[] {10, (int) 1e3, (int) 1e6, (int) 1e8};
    private static final int HUGE_DISTANCE = Integer.MAX_VALUE;
    private static final double EPSILON = 1e-6;

    //test various properties of a transition matrix
    @Test(dataProvider = "meanEventSizeAndEventStartProbability")
    public void markovianPropertiesTest(final double meanEventSize, final double eventStartProbability) {
        final CopyNumberTriStateTransitionProbabilityCache cache =
                new CopyNumberTriStateTransitionProbabilityCache(meanEventSize, eventStartProbability);

        for (final int d : DISTANCES) {
            //check symmetries -- these are part of the model and need not be true in the future
            final RealMatrix transitionMatrix = cache.getAsMatrixInProbabilitySpace(d);
            assertSymmetries(transitionMatrix);

            //check that columns sums equal 1
            for (int column = 0; column < transitionMatrix.getColumnDimension(); column++) {
                Assert.assertEquals(MathUtils.sum(transitionMatrix.getColumn(column)), 1, EPSILON);
            }

            //check that all elements are positive
            transitionMatrix.walkInOptimizedOrder(new DefaultRealMatrixPreservingVisitor() {
                @Override
                public void visit(int row, int column, double value) {
                    Assert.assertTrue(value >= 0);
                }
            });

            //check that T(2d) = T(d)*T(d)
            assertEqualMatrices(cache.getAsMatrixInProbabilitySpace(2*d), transitionMatrix.multiply(transitionMatrix));

            //check that the largest eigenvalue of the transition matrix is 1 (this corresponds to the asymptotic stationary state)
            Assert.assertEquals(MathUtils.arrayMax(new EigenDecomposition(transitionMatrix).getRealEigenvalues()), 1, EPSILON);
        }

        // check that at long distances memory of the initial state is lost and all initial distributions tend toward
        // the same asymptotic stationary distribution.  That is, all columns of the large-distance transition matrix are equal
        final RealMatrix asymptoticMatrix = cache.getAsMatrixInProbabilitySpace(HUGE_DISTANCE);
        for (int column = 1; column < asymptoticMatrix.getColumnDimension(); column++) {
            Assert.assertEquals(asymptoticMatrix.getColumnVector(0).subtract(asymptoticMatrix.getColumnVector(column)).getL1Norm(), 0, EPSILON);
        }
    }

    @DataProvider(name="meanEventSizeAndEventStartProbability")
    public Object[][] meanEventSizeAndEventStartProbability() {
        return new Object[][] {
                {1e3, 1e-8},
                {1e4, 1e-8},
                {1e5, 1e-6},
                {1e5, 1e-8},
                {1e3, 1e-5},
                {1e6, 1e-8}
        };
    }

    private static void assertEqualMatrices(final RealMatrix mat1, final RealMatrix mat2) {
        Assert.assertEquals(mat1.subtract(mat2).getNorm(), 0, EPSILON);
    }

    private void assertSymmetries(final RealMatrix transitionMatrix) {
        final int deletionIndex = CopyNumberTriState.DELETION.ordinal();
        final int duplicationIndex = CopyNumberTriState.DUPLICATION.ordinal();
        final int neutralIndex = CopyNumberTriState.NEUTRAL.ordinal();
        //neutral to deletion == neutral to duplication
        Assert.assertEquals(transitionMatrix.getEntry(deletionIndex, neutralIndex),
                transitionMatrix.getEntry(duplicationIndex, neutralIndex), 1e-10);

        // deletion to deletion == duplication to duplication
        Assert.assertEquals(transitionMatrix.getEntry(deletionIndex, deletionIndex),
                transitionMatrix.getEntry(duplicationIndex, duplicationIndex), 1e-10);

        // deletion to neutral == duplication to neutral
        Assert.assertEquals(transitionMatrix.getEntry(neutralIndex, deletionIndex),
                transitionMatrix.getEntry(neutralIndex, duplicationIndex), 1e-10);
    }
}