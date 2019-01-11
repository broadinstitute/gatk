package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class HDF5UtilsUnitTest extends GATKBaseTest {
    private static final double DOUBLE_MATRIX_TOLERANCE = 1E-4;
    private static final int CHUNK_DIVISOR = 256;
    private static final int MAX_CHUNK_SIZE = HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX / CHUNK_DIVISOR;

    @DataProvider(name = "testCreateLargeMatrixData")
    public Object[][] dataProvider() {
        return new Object[][] {
                new Object[] {100, 100000},
                new Object[] {CHUNK_DIVISOR + 16, MAX_CHUNK_SIZE},
                new Object[] {MAX_CHUNK_SIZE, CHUNK_DIVISOR + 16}
        };
    }

    //this test requires a large amount of heap space by design, so we disable it
    @Test(dataProvider = "testCreateLargeMatrixData", enabled = false)
    public void testCreateLargeMatrix(final int numRows,
                                      final int numColumns) {
        final double mean = 0.;
        final double sigma = 1.;
        final String matrixPath = "/test/matrix";

        final RealMatrix largeMatrix = createMatrixOfGaussianValues(numRows, numColumns, mean, sigma);
        final File tempOutputHD5 = IOUtils.createTempFile("large-matrix-", ".hd5");
        try (final HDF5File hdf5File = new HDF5File(tempOutputHD5, HDF5File.OpenMode.CREATE)) {
            HDF5Utils.writeChunkedDoubleMatrix(hdf5File, matrixPath, largeMatrix.getData(), MAX_CHUNK_SIZE);
        }

        try (final HDF5File hdf5FileForReading = new HDF5File(tempOutputHD5, HDF5File.OpenMode.READ_ONLY)) {
            final double[][] result = HDF5Utils.readChunkedDoubleMatrix(hdf5FileForReading, matrixPath);
            final RealMatrix resultAsRealMatrix = new Array2DRowRealMatrix(result, false);
            Assert.assertTrue(resultAsRealMatrix.getRowDimension() == numRows);
            Assert.assertTrue(resultAsRealMatrix.getColumnDimension() == numColumns);
            assertEqualsMatrix(resultAsRealMatrix, largeMatrix, DOUBLE_MATRIX_TOLERANCE);
        }
    }

    private static RealMatrix createMatrixOfGaussianValues(final int numRows,
                                                           final int numColumns,
                                                           final double mean,
                                                           final double sigma) {
        final RealMatrix largeMatrix = new Array2DRowRealMatrix(numRows, numColumns);
        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator();
        randomDataGenerator.reSeed(337337337);
        largeMatrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return randomDataGenerator.nextGaussian(mean, sigma);
            }
        });
        return largeMatrix;
    }

    /**
     * Test whether two matrices are equal
     * @param left never {@code null}
     * @param right never {@code null}
     */
    private static void assertEqualsMatrix(final RealMatrix left, final RealMatrix right, final double tolerance) {
        Assert.assertEquals(left.getRowDimension(), right.getRowDimension());
        Assert.assertEquals(left.getColumnDimension(), right.getColumnDimension());
        for (int i = 0; i < left.getRowDimension(); i++) {
            final double[] leftRow = left.getRow(i);
            final double[] rightRow = right.getRow(i);
            for (int j = 0; j < leftRow.length; j++) {
                Assert.assertEquals(leftRow[j], rightRow[j], tolerance);
            }
        }
    }
}
