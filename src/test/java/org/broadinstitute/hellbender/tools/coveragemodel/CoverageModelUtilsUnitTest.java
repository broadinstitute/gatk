package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link CoverageModelUtils}
 *
 * TODO github/gatk-protected issue #843
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class CoverageModelUtilsUnitTest extends BaseTest {

    @Test
    public void testEstimateReadDepthDensityFromPoissonModel() {
        Assert.assertEquals(CoverageModelUtils.estimateReadDepthFromPoissonModel(
                new double[] {10, 10, 10}, /* read counts */
                new double[] {1, 1, 1}, /* multiplicative biases */
                new int[] {2, 2, 2}, /* ploidies */
                new double[] {1, 1, 1}, /* copy ratio */
                new int[] {1, 1, 1}), /* mask */
                (FastMath.sqrt(401) - 1) / 4, 1e-12);

        Assert.assertEquals(CoverageModelUtils.estimateReadDepthFromPoissonModel(
                new double[] {10, 20, 30}, /* read counts */
                new double[] {1, 1, 1}, /* multiplicative biases */
                new int[] {2, 2, 2}, /* ploidies */
                new double[] {1, 1, 1}, /* copy ratio */
                new int[] {1, 1, 1}), /* mask */
                (FastMath.sqrt(8 * 1400 / 6.0 + 1) - 1) / 4, 1e-12);

        /* the zero ploidy entry must be masked */
        Assert.assertEquals(CoverageModelUtils.estimateReadDepthFromPoissonModel(
                new double[] {10, 20, 30}, /* read counts */
                new double[] {1, 1, 4567}, /* multiplicative biases */
                new int[] {2, 2, 0}, /* ploidies */
                new double[] {1, 1, 1234}, /* copy ratio */
                new int[] {1, 1, 1}), /* mask */
                (FastMath.sqrt(1001) - 1) / 4, 1e-12);
    }

    @Test(expectedExceptions = IllegalArgumentException.class, dataProvider = "estimateReadDepthExceptionsDataProvider")
    public void testEstimateReadDepthFromPoissonModelExceptions(final double[] readCounts,
                                                                final double[] multiplicativeBiases,
                                                                final int[] ploidies,
                                                                final double[] copyRatio,
                                                                final int[] mask) {
        CoverageModelUtils.estimateReadDepthFromPoissonModel(readCounts, multiplicativeBiases,
                ploidies, copyRatio, mask);
    }

    @DataProvider(name = "estimateReadDepthExceptionsDataProvider")
    public Object[][] estimateReadDepthExceptionsDataProvider() {
        return new Object[][] {{
                /* all masked */
                new double[] {10, 20, 30}, /* read counts */
                new double[] {1, 1, 4567}, /* multiplicative biases */
                new int[] {2, 2, 2}, /* ploidies */
                new double[] {1, 1, 1234}, /* copy ratio */
                new int[] {0, 0, 0} /* mask */
        }, {    /* all post-masked */
                new double[] {10, 20, 30}, /* read counts */
                new double[] {0.0, 1, 4567}, /* multiplicative biases */
                new int[] {2, 0, 0}, /* ploidies */
                new double[] {1, 1, 1234}, /* copy ratio */
                new int[] {1, 1, 1} /* mask */
        }, {    /* bad lengths */
                new double[] {10, 20, 30, 40}, /* read counts */
                new double[] {1, 1, 4567}, /* multiplicative biases */
                new int[] {2, 2, 0}, /* ploidies */
                new double[] {1, 1, 1234}, /* copy ratio */
                new int[] {1, 1, 1} /* mask */
        }};
    }
}