package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Assert;
import org.nd4j.linalg.factory.Nd4j;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link FourierLinearOperatorNDArray}
 *
 * TODO github/gatk-protected issue #843 -- create a unit test by generating sinusoidal test data and making sure
 * the filter works as intended
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */

public class FourierLinearOperatorNDArrayUnitTest extends BaseTest {

    private static final double EPS = 1e-8;

    @DataProvider(name = "testDataWithoutZeroPadding")
    public Object[][] getTestDataWithoutZeroPadding() {
        return new Object[][] {
            {4, /* dimension */
                new double[] {1.0, 1.0, 1.0}, /* Fourier factors */
                new double[] {1.0, 2.0, 3.0, 4.0}, /* input vector */
                new double[] {1.0, 2.0, 3.0, 4.0}}, /* transformed */
            {4,
                new double[] {0.0, 1.0, 1.0},
                new double[] {1.0, 2.0, 3.0, 4.0},
                new double[] {-1.5, -0.5, 0.5, 1.5}},
            {4,
                new double[] {-1.0, 3.0, 4.0},
                new double[] {1.0, 2.0, 3.0, 4.0},
                new double[] {-7.5, -3.5, -1.5, 2.5}},
            {5,
                new double[] {1.0, 1.0, 1.0},
                new double[] {1.0, 2.0, 3.0, 4.0, 5.0},
                new double[] {1.0, 2.0, 3.0, 4.0, 5.0}},
            {5,
                new double[] {0.0, 1.0, 1.0},
                new double[] {1.0, 2.0, 3.0, 4.0, 5.0},
                new double[] {-2.0, -1.0, 0.0, 1.0, 2.0}},
            {5,
                new double[] {1.0, 2.0, 1.0},
                new double[] {1.0, 2.0, 3.0, 4.0, 5.0},
                new double[] {0, 3.81966011e-01, 3.0, 5.61803399e+00, 6.0}}
        };
    }

    @DataProvider(name = "testDataWithZeroPadding")
    public Object[][] getTestDataWithZeroPadding() {
        return new Object[][] {
                {4, /* dimension */
                        new double[] {1.0, 1.0, 1.0}, /* Fourier factors */
                        new double[] {1.0, 2.0, 3.0, 4.0}, /* input vector */
                        new double[] {1.0, 2.0, 3.0, 4.0}}, /* transformed */
                {4,
                        new double[] {0.0, 1.0, 1.0},
                        new double[] {1.0, 2.0, 3.0, 4.0},
                        new double[] {-1.5, -0.5, 0.5, 1.5}},
                {4,
                        new double[] {-1.0, 3.0, 4.0},
                        new double[] {1.0, 2.0, 3.0, 4.0},
                        new double[] {-7.5, -3.5, -1.5, 2.5}},
                {5,
                        new double[] {1.0, 1.0, 1.0},
                        new double[] {1.0, 2.0, 3.0, 4.0, 5.0},
                        new double[] {1.0, 2.0, 3.0, 4.0, 5.0}},
                {5,
                        new double[] {1.0, 1.0, 1.0},
                        new double[] {1.0, 2.0, 3.0, 4.0, 5.0},
                        new double[] {1.0, 2.0, 3.0, 4.0, 5.0}},
                {5,
                        new double[] {0.0, 1.0, 1.0},
                        new double[] {1.0, 2.0, 3.0, 4.0, 5.0},
                        new double[] {1.0 - 15.0/8.0, 2.0 - 15.0/8.0, 3.0 - 15.0/8.0, 4.0 - 15.0/8.0, 5.0 - 15.0/8.0}},
                {5,
                        new double[] {1.0, 2.0, 1.0},
                        new double[] {1.0, 2.0, 3.0, 4.0, 5.0},
                        new double[] {0.39644661, 1.82322330, 4.06066017, 6.73743687, 7.10355339}}
        };
    }

    /**
     * The main test routine
     *
     * @param dim dimension of the operator
     * @param fourierFactors Fourier factors
     * @param x input data
     * @param y expected output data
     */
    @Test(dataProvider = "testDataWithoutZeroPadding")
    public void performTestWithoutZeroPadding(final int dim, final double[] fourierFactors, final double[] x, final double[] y) {
        FourierLinearOperatorNDArray linOp = new FourierLinearOperatorNDArray(dim, fourierFactors, false);
        final double[] yCalc = linOp.operate(Nd4j.create(x)).data().asDouble();
        Assert.assertArrayEquals(y, yCalc, EPS);
    }

    /**
     * The main test routine
     *
     * @param dim dimension of the operator
     * @param fourierFactors Fourier factors
     * @param x input data
     * @param y expected output data
     */
    @Test(dataProvider = "testDataWithZeroPadding")
    public void performTestWithZeroPadding(final int dim, final double[] fourierFactors, final double[] x, final double[] y) {
        FourierLinearOperatorNDArray linOp = new FourierLinearOperatorNDArray(dim, fourierFactors, true);
        final double[] yCalc = linOp.operate(Nd4j.create(x)).data().asDouble();
        Assert.assertArrayEquals(y, yCalc, EPS);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadDimension_0() {
        /* negative dimensions not allowed */
        new FourierLinearOperatorNDArray(-1, new double[5], false);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadDimension_1() {
        /* dimension >= 2 */
        new FourierLinearOperatorNDArray(1, new double[5], false);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadDimension_2() {
        /* fourierFactors.length = floor(dimension/2) + 1 */
        new FourierLinearOperatorNDArray(15, new double[3], false);
    }
}