package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link CoverageModelEMWorkspaceMathUtils}
 *
 * TODO github/gatk-protected issue #843
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class CoverageModelEMWorkspaceMathUtilsUnitTest extends BaseTest {
    private static final double EPS = 1e-12;

    @Test
    public void testLogDet() {
        assertLogDetCorrectness(Nd4j.create(new double[][]{
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1}}), 0.0);
        assertLogDetCorrectness(Nd4j.create(new double[][]{
                {1, 0, 0},
                {0, -1, 0},
                {0, 0, 1}}), 0.0);
        assertLogDetCorrectness(Nd4j.create(new double[][]{
                {1, 4, 20},
                {0, -2, -10},
                {0, 0, 3}}), FastMath.log(6));
        assertLogDetCorrectness(Nd4j.create(new double[][]{
                {1, 10},
                {4, 5}}), FastMath.log(35));
    }

    public void assertLogDetCorrectness(final INDArray mat, final double expected) {
        final double logDet = CoverageModelEMWorkspaceMathUtils.logdet(mat);
        Assert.assertEquals(logDet, expected, EPS);
    }
}
