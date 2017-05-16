package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.tools.coveragemodel.cachemanager.ComputableNodeFunction;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

/**
 * Unit tests for {@link CoverageModelEMComputeBlock}
 *
 * TODO github/gatk-protected issue #843 -- unit tests to ensure that all {@link ComputableNodeFunction}
 * and all functions annotated with {@link org.broadinstitute.hellbender.tools.coveragemodel.CoverageModelEMComputeBlock.QueriesICG}
 * treat ICG node values
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class CoverageModelEMComputeBlockUnitTest extends BaseTest {

    private static final double EPSILON = 1e-12;

    @Test
    public void testReplaceMaskedEntries() {
        final double CORRECTION_VALUE = 5.0;

        /* a tainted array */
        final INDArray arr = Nd4j.create(new double[][]{
                {1, Double.NaN, 3, 18, 3, -5},
                {4, 0, Double.POSITIVE_INFINITY, 2, 4, 12},
                {Double.NaN, 8, 9, Double.NEGATIVE_INFINITY, -12, Double.MAX_VALUE}});

        /* mask on some values */
        final INDArray mask = Nd4j.create(new double[][]{
                {1, 0, 1, 1, 1, 1},
                {0, 1, 0, 1, 1, 1},
                {0, 1, 1, 0, 1, 0}});

        /* corrected array */
        final INDArray cor = Nd4j.create(new double[][]{
                {1, CORRECTION_VALUE, 3, 18, 3, -5},
                {CORRECTION_VALUE, 0, CORRECTION_VALUE, 2, 4, 12},
                {CORRECTION_VALUE, 8, 9, CORRECTION_VALUE, -12, CORRECTION_VALUE}});

        /* test on full view correction */
        assertNDArrayEquals(CoverageModelEMComputeBlock.replaceMaskedEntries(arr, mask, CORRECTION_VALUE), cor);

        /* test on partial views */
        assertNDArrayEquals(CoverageModelEMComputeBlock.replaceMaskedEntries(
                arr.get(NDArrayIndex.interval(0, 2), NDArrayIndex.all()),
                mask.get(NDArrayIndex.interval(0, 2), NDArrayIndex.all()),
                CORRECTION_VALUE),
                cor.get(NDArrayIndex.interval(0, 2), NDArrayIndex.all()));
        assertNDArrayEquals(CoverageModelEMComputeBlock.replaceMaskedEntries(
                arr.get(NDArrayIndex.interval(2, 5), NDArrayIndex.interval(1, 4)),
                mask.get(NDArrayIndex.interval(2, 5), NDArrayIndex.interval(1, 4)),
                CORRECTION_VALUE),
                cor.get(NDArrayIndex.interval(2, 5), NDArrayIndex.interval(1, 4)));
        assertNDArrayEquals(CoverageModelEMComputeBlock.replaceMaskedEntries(
                arr.get(NDArrayIndex.all(), NDArrayIndex.interval(2, 3)),
                mask.get(NDArrayIndex.all(), NDArrayIndex.interval(2, 3)),
                CORRECTION_VALUE),
                cor.get(NDArrayIndex.all(), NDArrayIndex.interval(2, 3)));
    }

    private void assertNDArrayEquals(final INDArray arr1, final INDArray arr2) {
        ArrayAsserts.assertArrayEquals(arr1.dup().data().asDouble(), arr2.dup().data().asDouble(), EPSILON);
    }
}
