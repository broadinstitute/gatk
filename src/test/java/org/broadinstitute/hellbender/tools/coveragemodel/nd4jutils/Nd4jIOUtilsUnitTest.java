package org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;

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
    public void tsvReadWriteTest() {
        tsvReadWriteAsserter(Nd4j.rand(10, 10));
        tsvReadWriteAsserter(Nd4j.rand(1, 10));
        tsvReadWriteAsserter(Nd4j.rand(10, 1));
    }

    public void tsvReadWriteAsserter(final INDArray arr) {
        final File outFile = createTempFile("Nd4j_tsv_test", ".tsv");
        Nd4jIOUtils.writeNDArrayToTextFile(arr, outFile, null, null);
        final INDArray arrLoaded = Nd4jIOUtils.readNDArrayFromTextFile(outFile);
        ArrayAsserts.assertArrayEquals(arr.shape(), arrLoaded.shape());
        ArrayAsserts.assertArrayEquals(arr.data().asDouble(), arrLoaded.data().asDouble(), 1e-12);
    }

    public void binaryReadWriteAsserter(final INDArray arr) {
        final File outFile = createTempFile("Nd4j_binary_dump_test", ".nd4j");
        Nd4jIOUtils.writeNDArrayToBinaryDumpFile(arr, outFile);
        final INDArray arrLoaded = Nd4jIOUtils.readNDArrayFromBinaryDumpFile(outFile);
        ArrayAsserts.assertArrayEquals(arr.shape(), arrLoaded.shape());
        ArrayAsserts.assertArrayEquals(arr.data().asDouble(), arrLoaded.data().asDouble(), 1e-12);
    }
}
