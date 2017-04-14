package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.buffer.DataBuffer;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * This test is to ensure that the DType is correctly set to Double
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class Nd4jUnitTest extends BaseTest {
    private static final double EPS = 1e-10;

    @Test
    public void testDType() {
        Assert.assertTrue(Nd4j.dataType().equals(DataBuffer.Type.DOUBLE), "Data type for Nd4j must be set to" +
                " double; otherwise, coverage model EM algorithm will not function properly");
    }

    /**
     * A test for implementing a multiplication like:
     *
     *      X_{b c} = \sum_{a} W_{a b c} v_{a}
     *
     * using matrix products and successive reshapes.
     */
    @Test
    public void testMulTensorVector() {
        /* generate random data */
        final int A = 5;
        final int B = 6;
        final int C = 7;
        final INDArray W = Nd4j.rand(new int[] {A, B, C});
        final INDArray v = Nd4j.rand(new int[] {A, 1});

        /* result using reshapes and matrix products */
        final INDArray X = W.reshape(new int[] {A, B*C}).transpose().mmul(v).reshape(new int[] {B, C});

        /* check against brute force result */
        for (int b = 0; b < B; b++) {
            for (int c = 0; c < C; c++) {
                double prod = 0;
                for (int a = 0; a < A; a++) {
                    prod += W.get(NDArrayIndex.point(a), NDArrayIndex.point(b), NDArrayIndex.point(c)).getDouble(0) *
                            v.getDouble(a);
                }
                Assert.assertEquals(X.getScalar(b, c).getDouble(0), prod, EPS);
            }
        }
    }

    /**
     * A test for implementing a multiplication like:
     *
     *      X_{a} = \sum_{b, c} W_{a b c} V_{b c}
     */
    @Test
    public void testMulTensorMatrix() {
        /* generate random data */
        final int A = 5;
        final int B = 6;
        final int C = 7;
        final INDArray W = Nd4j.rand(new int[] {A, B, C});
        final INDArray V = Nd4j.rand(new int[] {B, C});

        /* result using reshapes and matrix products */
        final INDArray X = W.reshape(new int[] {A, B*C}).mmul(V.reshape(new int[] {B*C, 1}));

        /* check against brute force result */
        for (int a = 0; a < A; a++) {
            double prod = 0;
            for (int b = 0; b < B; b++) {
                for (int c = 0; c < C; c++) {
                    prod += W.get(NDArrayIndex.point(a), NDArrayIndex.point(b), NDArrayIndex.point(c)).getDouble(0) *
                            V.get(NDArrayIndex.point(b), NDArrayIndex.point(c)).getDouble(0);
                }
            }
            Assert.assertEquals(X.getScalar(a).getDouble(0), prod, EPS);
        }
    }
}
