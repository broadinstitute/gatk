package org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;

import static org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jIOUtils.getIntArrayFromStringRepr;
import static org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jIOUtils.getStringReprFromIntArray;

/**
 * Unit tests for {@link Nd4jIOUtils}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class Nd4jIOUtilsUnitTest extends BaseTest {

    @Test
    public void binaryDumpTest() {
        binaryReadWriteAsserter(Nd4j.rand(10, 10));
        binaryReadWriteAsserter(Nd4j.rand(1, 10));
        binaryReadWriteAsserter(Nd4j.rand(10, 1));
    }

    @Test
    public void tsvMatrixReadWriteTest() {
        tsvMatrixReadWriteAsserter(Nd4j.rand(10, 10));
        tsvMatrixReadWriteAsserter(Nd4j.rand(1, 10));
        tsvMatrixReadWriteAsserter(Nd4j.rand(10, 1));
    }

    @Test
    public void tsvTensorReadWriteTest() {
        tsvTensorReadWriteAsserter(Nd4j.rand(new int[] {10, 10}));
        tsvTensorReadWriteAsserter(Nd4j.rand(new int[] {1, 10}));
        tsvTensorReadWriteAsserter(Nd4j.rand(new int[] {10, 1}));
        tsvTensorReadWriteAsserter(Nd4j.rand(new int[] {10, 10, 10}));
        tsvTensorReadWriteAsserter(Nd4j.rand(new int[] {1, 10, 15}));
        tsvTensorReadWriteAsserter(Nd4j.rand(new int[] {5, 3, 8}));
    }

    @Test
    public void testIntArrayStringRepr() {
        assertIntArrayStringRepr(new int[] {1});
        assertIntArrayStringRepr(new int[] {1, -5, 2});
        assertIntArrayStringRepr(new int[] {4, 3});
        assertIntArrayStringRepr(new int[] {1, 2, 3, 4, 5, 6});

    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntArrayStringReprEmpty() {
        assertIntArrayStringRepr(new int[] {});
    }

    public void assertIntArrayStringRepr(final int[] indices) {
        final String repr = getStringReprFromIntArray(indices);
        final int[] reconstructed = getIntArrayFromStringRepr(repr);
        ArrayAsserts.assertArrayEquals(indices, reconstructed);
    }

    public void tsvMatrixReadWriteAsserter(final INDArray arr) {
        final File outFile = createTempFile("Nd4j_matrix_tsv_test", ".tsv");
        Nd4jIOUtils.writeNDArrayMatrixToTextFile(arr, outFile, "TEST", null, null);
        final INDArray arrLoaded = Nd4jIOUtils.readNDArrayMatrixFromTextFile(outFile);
        ArrayAsserts.assertArrayEquals(arr.shape(), arrLoaded.shape());
        ArrayAsserts.assertArrayEquals(arr.dup().data().asDouble(), arrLoaded.dup().data().asDouble(), 1e-12);
    }

    public void binaryReadWriteAsserter(final INDArray arr) {
        final File outFile = createTempFile("Nd4j_binary_dump_test", ".nd4j");
        Nd4jIOUtils.writeNDArrayToBinaryDumpFile(arr, outFile);
        final INDArray arrLoaded = Nd4jIOUtils.readNDArrayFromBinaryDumpFile(outFile);
        ArrayAsserts.assertArrayEquals(arr.shape(), arrLoaded.shape());
        ArrayAsserts.assertArrayEquals(arr.dup().data().asDouble(), arrLoaded.dup().data().asDouble(), 1e-12);
    }

    public void tsvTensorReadWriteAsserter(final INDArray arr) {
        final File outFile = createTempFile("Nd4j_tensor_tsv_test", ".tsv");
        Nd4jIOUtils.writeNDArrayTensorToTextFile(arr, outFile, "TEST", null);
        final INDArray arrLoaded = Nd4jIOUtils.readNDArrayTensorFromTextFile(outFile);
        ArrayAsserts.assertArrayEquals(arr.shape(), arrLoaded.shape());
        ArrayAsserts.assertArrayEquals(arr.dup().data().asDouble(), arrLoaded.dup().data().asDouble(), 1e-12);
    }

}
