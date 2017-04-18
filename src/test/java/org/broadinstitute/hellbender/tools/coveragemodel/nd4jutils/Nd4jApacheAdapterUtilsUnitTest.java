package org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

/**
 * Unit tests for {@link Nd4jApacheAdapterUtils}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class Nd4jApacheAdapterUtilsUnitTest extends BaseTest {

    private static final double EPS = 1e-12;

    @Test
    public void testApacheMatrixToINDArray() {
        /* row vector edge case */
        assertApacheMatrixToINDArrayCorrectness(
                new BlockRealMatrix(new double[][] {{1.0, 2.0, 3.0}}));

        /* column vector edge case */
        assertApacheMatrixToINDArrayCorrectness(
                new BlockRealMatrix(new double[][] {{1.0}, {2.0}, {3.0}}));

        /* general matrix */
        assertApacheMatrixToINDArrayCorrectness(
                new BlockRealMatrix(new double[][] {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}));
    }

    @Test
    public void testINDArrayToApacheMatrix() {
        final INDArray rowArrCOrder = Nd4j.randn('c', new int[] {1, 4});
        final INDArray rowArrFOrder = Nd4j.randn('f', new int[] {1, 4});
        final INDArray colArrCOrder = Nd4j.randn('c', new int[] {4, 1});
        final INDArray colArrFOrder = Nd4j.randn('f', new int[] {4, 1});
        final INDArray generalCOrder = Nd4j.randn('c', new int[] {4, 5});
        final INDArray generalFOrder = Nd4j.randn('f', new int[] {4, 5});

        assertINDArrayToApacheMatrixCorrectness(rowArrCOrder);
        assertINDArrayToApacheMatrixCorrectness(rowArrFOrder);
        assertINDArrayToApacheMatrixCorrectness(colArrCOrder);
        assertINDArrayToApacheMatrixCorrectness(colArrFOrder);
        assertINDArrayToApacheMatrixCorrectness(generalCOrder);
        assertINDArrayToApacheMatrixCorrectness(generalFOrder);

        /* test on INDArray views */
        assertINDArrayToApacheMatrixCorrectness(rowArrCOrder.get(NDArrayIndex.all(), NDArrayIndex.interval(0, 3)));
        assertINDArrayToApacheMatrixCorrectness(rowArrFOrder.get(NDArrayIndex.all(), NDArrayIndex.interval(0, 3)));
        assertINDArrayToApacheMatrixCorrectness(colArrCOrder.get(NDArrayIndex.interval(0, 3), NDArrayIndex.all()));
        assertINDArrayToApacheMatrixCorrectness(colArrFOrder.get(NDArrayIndex.interval(0, 3), NDArrayIndex.all()));
        assertINDArrayToApacheMatrixCorrectness(generalCOrder.get(NDArrayIndex.interval(1, 4), NDArrayIndex.interval(2, 4)));
        assertINDArrayToApacheMatrixCorrectness(generalFOrder.get(NDArrayIndex.interval(1, 4), NDArrayIndex.interval(2, 4)));
    }

    @Test
    public void testINDArrayToApacheVector() {
        final INDArray rowArrCOrder = Nd4j.randn('c', new int[] {1, 5});
        final INDArray rowArrFOrder = Nd4j.randn('f', new int[] {1, 5});
        final INDArray colArrCOrder = Nd4j.randn('c', new int[] {5, 1});
        final INDArray colArrFOrder = Nd4j.randn('f', new int[] {5, 1});

        assertINDArrayToApacheVectorCorrectness(rowArrCOrder);
        assertINDArrayToApacheVectorCorrectness(rowArrFOrder);
        assertINDArrayToApacheVectorCorrectness(colArrCOrder);
        assertINDArrayToApacheVectorCorrectness(colArrFOrder);

        /* test on INDArray views */
        assertINDArrayToApacheVectorCorrectness(rowArrCOrder.get(NDArrayIndex.all(), NDArrayIndex.interval(2, 4)));
        assertINDArrayToApacheVectorCorrectness(rowArrFOrder.get(NDArrayIndex.all(), NDArrayIndex.interval(2, 4)));
        assertINDArrayToApacheVectorCorrectness(colArrCOrder.get(NDArrayIndex.interval(2, 4), NDArrayIndex.all()));
        assertINDArrayToApacheVectorCorrectness(colArrFOrder.get(NDArrayIndex.interval(2, 4), NDArrayIndex.all()));
    }

    public void assertINDArrayToApacheMatrixCorrectness(final INDArray arr) {
        /* brute force result */
        final RealMatrix expected = new BlockRealMatrix(arr.rows(), arr.columns());
        for (int i = 0; i < arr.rows(); i++) {
            for (int j = 0; j < arr.columns(); j++) {
                expected.setEntry(i, j, arr.getDouble(i, j));
            }
        }
        final RealMatrix result = Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(arr);
        Assert.assertEquals(result.getColumnDimension(), expected.getColumnDimension());
        Assert.assertEquals(result.getRowDimension(), expected.getRowDimension());
        for (int i = 0; i < expected.getRowDimension(); i++) {
            ArrayAsserts.assertArrayEquals(result.getData()[i], expected.getData()[i], EPS);
        }
    }

    public void assertApacheMatrixToINDArrayCorrectness(final RealMatrix mat) {
        final INDArray result = Nd4jApacheAdapterUtils.convertApacheMatrixToINDArray(mat);
        for (int i = 0; i < mat.getRowDimension(); i++) {
            for (int j = 0; j < mat.getColumnDimension(); j++) {
                Assert.assertEquals(result.getDouble(i, j), mat.getEntry(i, j), EPS);
            }
        }
    }

    public void assertINDArrayToApacheVectorCorrectness(final INDArray arr) {
        final RealVector result = Nd4jApacheAdapterUtils.convertINDArrayToApacheVector(arr);
        final RealVector expected = new ArrayRealVector(arr.dup().data().asDouble());
        Assert.assertEquals(result.getDimension(), expected.getDimension());
        ArrayAsserts.assertArrayEquals(result.toArray(), expected.toArray(), EPS);
    }
}
