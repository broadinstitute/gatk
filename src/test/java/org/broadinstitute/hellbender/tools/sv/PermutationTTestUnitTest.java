package org.broadinstitute.hellbender.tools.sv;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class PermutationTTestUnitTest {

    @DataProvider(name = "testData")
    public Object[][] getTestData() {
        return new Object[][]{
                {
                    new double[]{},
                    new double[]{},
                    PermutationTTest.Alternative.GREATER,
                    Double.NaN,
                    Double.NaN
                },
                {
                        new double[]{0},
                        new double[]{},
                        PermutationTTest.Alternative.GREATER,
                        Double.NaN,
                        Double.NaN
                },
                {
                        new double[]{},
                        new double[]{0},
                        PermutationTTest.Alternative.GREATER,
                        Double.NaN,
                        Double.NaN
                },
                {
                        new double[]{0},
                        new double[]{0},
                        PermutationTTest.Alternative.GREATER,
                        Double.NaN,
                        Double.NaN
                },
                {
                        new double[]{1.},
                        new double[]{1.},
                        PermutationTTest.Alternative.GREATER,
                        Double.NaN,
                        Double.NaN
                },
                {
                        new double[]{0, 0},
                        new double[]{0, 0},
                        PermutationTTest.Alternative.GREATER,
                        Double.NaN,
                        Double.NaN
                },
                {
                        new double[]{0, 1},
                        new double[]{0, 1},
                        PermutationTTest.Alternative.GREATER,
                        0.,
                        0.5
                },
                {
                        new double[]{0, 1},
                        new double[]{2, 3},
                        PermutationTTest.Alternative.GREATER,
                        -1.5491933384829666,
                        0.9393323748207589
                },
                {
                        new double[]{2, 3, 4},
                        new double[]{3, 4, 5},
                        PermutationTTest.Alternative.GREATER,
                        -1.1677484162422844,
                        0.8785458695478379
                },
                {
                        new double[]{2, 3, 4},
                        new double[]{3, 4, 5},
                        PermutationTTest.Alternative.TWO_SIDED,
                        -1.1677484162422844,
                        0.24290826090432427
                },
        };
    }

    // Tests GREATER or TWO_SIDED; if GREATER, automatically tests LESS as well
    @Test(dataProvider = "testData")
    public void test(final double[] x, final double[] y, final PermutationTTest.Alternative alternative, final double expectedStat, final double expectedP) {
        final PermutationTTest.PermTestResult result = PermutationTTest.test(x, y, alternative);
        SVTestUtils.assertFloatWithinTolerance(result.statistic(), expectedStat);
        SVTestUtils.assertFloatWithinTolerance(result.pValue(), expectedP);
        if (alternative == PermutationTTest.Alternative.GREATER) {
            final PermutationTTest.PermTestResult resultLess = PermutationTTest.test(x, y, PermutationTTest.Alternative.LESS);
            SVTestUtils.assertFloatWithinTolerance(resultLess.statistic(), expectedStat);
            SVTestUtils.assertFloatWithinTolerance(resultLess.pValue(), 1 - expectedP);
        }
    }
}