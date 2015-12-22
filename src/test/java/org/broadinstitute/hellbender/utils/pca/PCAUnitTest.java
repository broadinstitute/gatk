package org.broadinstitute.hellbender.utils.pca;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.svd.SVDFactory;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link PCA}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class PCAUnitTest extends BaseTest {

    /**
     * Precision required when comparing eigen vector values.
     */
    public static final double EPSILON = 0.0001;

    /**
     * The test matrix.
     */
    public static final double[][] TEST_MATRIX = {
            {-1.53010281, -0.8839249, -0.007075833, 0.04151694},
            {1.57662029, -0.7109790, 0.010064938, 0.85383808},
            {-0.70763334, -0.6530172, -0.944054729, 1.14103165},
            {-1.03426795, -0.4949822, 1.289649707, -0.99308528},
            {-0.82874316, -2.0274880, 0.314147081, -0.18261399},
            {0.01674156, 3.1476897, -0.411199655, 0.24109695},
    };

    /**
     * Obtained using {@link #TEST_MATRIX} in R:
     *
     * <code>
     *     prcomp(t(TEST_MATRIX))$rotation
     * </code>
     */
    public static final double[][] TEST_EXPECTED_EIGENVECTORS = {
            {0.14501664, 0.22996366, 0.57396632, -0.41715764},
            {0.30390126, -0.50045292, -0.28479337, 0.14382405},
            {0.08699514, -0.38386777,  0.722709780,  0.53997803},
            {0.16341388, 0.72446041, -0.04526613, 0.60138021},
            {0.47639747,  0.14913586,  0.18374377, -0.38383406},
            {-0.79081439, 0.04716718, 0.17665770, -0.06878317}
    };

    /**
     * Obtained using {@link #TEST_MATRIX} in R:
     *
     * <code>
     *     prcomp(t(TEST_MATRIX))$centers
     * </code>
     */
    public static final double[] TEST_EXPECTED_CENTERS = new double[] {
            -0.5948967, 0.4323861, -0.2909184, -0.3081714, -0.6811745, 0.7485821
    };

    /**
     * Obtained using {@link #TEST_MATRIX} in R:
     *
     * <code>
     *     prcomp(t(TEST_MATRIX))$sdevs
     * </code>
     */
    public static final double[] TEST_EXPECTED_SDEVS = new double[] {
            2.035583e+00, 1.435669e+00, 1.064676e+00, 2.195049e-16
    };

    /**
     * {@link #TEST_MATRIX} principal component variances.
     */
    public static final double[] TEST_EXPECTED_VARS = DoubleStream.of(TEST_EXPECTED_SDEVS)
            .map(d -> d * d).toArray();

    @Test()
    public void testResidentSVDOnTestMatrix() {
        final RealMatrix dataMatrix = new Array2DRowRealMatrix(TEST_MATRIX);
        final List<String> variables = createVariableNames(TEST_MATRIX.length);
        final List<String> sampleNames = createSampleNames(TEST_MATRIX[0].length);
        final PCA pca = PCA.createPCA(variables, sampleNames, dataMatrix, SVDFactory::createSVD);
        final RealVector centers = pca.getCenters();
        final RealVector variances = pca.getVariances();
        final RealMatrix eigenVectors = pca.getEigenVectors();
        Assert.assertNotNull(eigenVectors);
        assertEquals(centers, new ArrayRealVector(TEST_EXPECTED_CENTERS), EPSILON);
        assertEquals(variances, new ArrayRealVector(TEST_EXPECTED_VARS), EPSILON);
        assertEqualEigenVectors(eigenVectors, new Array2DRowRealMatrix(TEST_EXPECTED_EIGENVECTORS), variances, EPSILON);
    }

    @Test()
    public void testResidentSVDOnTestMatrixWithoutVariableNorSampleNames() {
        final RealMatrix dataMatrix = new Array2DRowRealMatrix(TEST_MATRIX);
        final PCA pca = PCA.createPCA(dataMatrix, SVDFactory::createSVD);
        final RealVector centers = pca.getCenters();
        final RealVector variances = pca.getVariances();
        final RealMatrix eigenVectors = pca.getEigenVectors();
        Assert.assertNotNull(eigenVectors);
        assertEquals(centers, new ArrayRealVector(TEST_EXPECTED_CENTERS), EPSILON);
        assertEquals(variances, new ArrayRealVector(TEST_EXPECTED_VARS), EPSILON);
        assertEqualEigenVectors(eigenVectors, new Array2DRowRealMatrix(TEST_EXPECTED_EIGENVECTORS), variances, EPSILON);
    }

    private static List<String> createVariableNames(final int size) {
        return IntStream.range(0, size).mapToObj(i -> "TARGET_" + i).collect(Collectors.toList());
    }

    private static List<String> createSampleNames(final int size) {
        return IntStream.range(0, size).mapToObj(i -> "SAMPLE_" + i).collect(Collectors.toList());
    }

    @Test()
    public void testSparkSVDOnTestMatrix() {
        final RealMatrix dataMatrix = new Array2DRowRealMatrix(TEST_MATRIX);
        final JavaSparkContext sparkContext = SparkContextFactory.getSparkContext("test", Collections.emptyMap(), SparkContextFactory.DEFAULT_SPARK_MASTER);
        final List<String> variables = createVariableNames(TEST_MATRIX.length);
        final List<String> sampleNames = createSampleNames(TEST_MATRIX[0].length);
        final PCA pca = PCA.createPCA(variables, sampleNames, dataMatrix, dm -> SVDFactory.createSVD(dm, sparkContext));
        final RealVector centers = pca.getCenters();
        final RealVector variances = pca.getVariances();
        final RealMatrix eigenVectors = pca.getEigenVectors();
        Assert.assertNotNull(eigenVectors);
        assertEquals(centers, new ArrayRealVector(TEST_EXPECTED_CENTERS), EPSILON);
        assertEquals(variances, new ArrayRealVector(TEST_EXPECTED_VARS), EPSILON);
        assertEqualEigenVectors(eigenVectors, new Array2DRowRealMatrix(TEST_EXPECTED_EIGENVECTORS), variances, EPSILON);
    }

    public static void assertEquals(final RealVector v1, final RealVector v2, final double epsilon) {
        Assert.assertEquals(v1.getDimension(), v2.getDimension());
        for (int i = 0; i < v1.getDimension(); i++) {
            Assert.assertEquals(v1.getEntry(i), v2.getEntry(i), epsilon, " Differences between " + Arrays.toString(v1.toArray()) + " and " + Arrays.toString(v2.toArray()));
        }
    }

    public static void assertEqualEigenVectors(final RealMatrix m1, final RealMatrix m2, final RealVector variances, final double epsilon) {
        Assert.assertEquals(m1.getRowDimension(), m2.getRowDimension());
        Assert.assertEquals(m1.getColumnDimension(), m2.getColumnDimension());
        final boolean[] sameSign = new boolean[m1.getColumnDimension()];

        // Find out whether the eigen-vector is the same or inverted (which is equivalent).
        for (int j = 0; j < sameSign.length; j++ ) {
            for (int i = 0; i < m1.getRowDimension(); i++) {
                final double product = m1.getEntry(0, j) * m2.getEntry(0, j);
                if (product != 0) {
                    sameSign[j] = product > 0;
                    break;
                }
            }
        }
        for (int i = 0; i < m1.getRowDimension(); i++) {
            for (int j = 0; j < m1.getColumnDimension(); j++) {
                // multiplied with the variance to ignore very low variance components that might
                // have divergent axis due to computation instability.
                final double v1 = m1.getEntry(i, j) * variances.getEntry(j);
                final double v2 = m2.getEntry(i, j) * variances.getEntry(j);
                Assert.assertEquals(v1, sameSign[j] ? v2 : - v2, epsilon, " Differences between " + Arrays.toString(m1.getRow(i)) + " and " + Arrays.toString(m2.getRow(i)));
            }
        }
    }
}
