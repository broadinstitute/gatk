package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Assert;
import org.nd4j.linalg.factory.Nd4j;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link FourierLinearOperatorNDArray}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */

public class FourierLinearOperatorNDArrayUnitTest extends BaseTest {

    @DataProvider(name = "testDataWithoutZeroPadding")
    public Object[][] getTestDataWithoutZeroPadding() {
        return new Object[][] {
            {4,
                new double[] {1.0, 1.0, 1.0},
                new double[] {1.0, 2.0, 3.0, 4.0},
                new double[] {1.0, 2.0, 3.0, 4.0}},
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
                {4,
                        new double[] {1.0, 1.0, 1.0},
                        new double[] {1.0, 2.0, 3.0, 4.0},
                        new double[] {1.0, 2.0, 3.0, 4.0}},
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
     * @param dim dimension of the operator
     * @param fourierFacts Fourier factors
     * @param x input data
     * @param y expected output data
     */
    @Test(dataProvider = "testDataWithoutZeroPadding")
    public void performTestWithoutZeroPadding(final int dim, final double[] fourierFacts, final double[] x, final double[] y) {
        FourierLinearOperatorNDArray linOp = new FourierLinearOperatorNDArray(dim, fourierFacts, false);
        final double[] yCalc = linOp.operate(Nd4j.create(x)).data().asDouble();
        Assert.assertArrayEquals("", y, yCalc, 1e-8);
    }

    /**
     * The main test routine
     * @param dim dimension of the operator
     * @param fourierFacts Fourier factors
     * @param x input data
     * @param y expected output data
     */
    @Test(dataProvider = "testDataWithZeroPadding")
    public void performTestWithZeroPadding(final int dim, final double[] fourierFacts, final double[] x, final double[] y) {
        FourierLinearOperatorNDArray linOp = new FourierLinearOperatorNDArray(dim, fourierFacts, true);
        final double[] yCalc = linOp.operate(Nd4j.create(x)).data().asDouble();
        Assert.assertArrayEquals("", y, yCalc, 1e-8);
    }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadDimension_0() {
        /* negative dimensions not allowed */
        FourierLinearOperatorNDArray linop = new FourierLinearOperatorNDArray(-1, new double[5], false);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadDimension_1() {
        /* dimension >= 2 */
        FourierLinearOperatorNDArray linop = new FourierLinearOperatorNDArray(1, new double[5], false);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadDimension_2() {
        /* fourierFactors.length = floor(dimension/2) + 1 */
        FourierLinearOperatorNDArray linop = new FourierLinearOperatorNDArray(15, new double[3], false);
    }

}