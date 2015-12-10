package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.SparkToggleCommandLineProgram;
import org.broadinstitute.hellbender.utils.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.pca.PCA;
import org.broadinstitute.hellbender.utils.pca.PCAUnitTest;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.utils.pca.PCAUnitTest.TEST_MATRIX;

/**
 * Integration tests for {@link CalculateCoverageComponents} class.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class CalculateCoverageComponentsIntegrationTest extends CommandLineProgramTest {

    @Test()
    public void testSimpleMatrixNoSpark() {
        final File inputFile = createTempFile("ccc-test",".tab");
        final File outputFile = createTempFile("ccc-test-out",".hd5");
        writeReadCounts(inputFile, TEST_MATRIX);
        runCommandLine(new String[]{"-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + SparkToggleCommandLineProgram.DISABLE_SPARK_SHORT_NAME});

        checkPCAOutputFile(outputFile);
    }

    //TODO Once we fix the problem when using Spark for SVD/PCA, we can enable the test again.
    //TODO Corresponding issue number #242.
    @Test(enabled = false)
    public void testSimpleMatrixSpark() {
        final File inputFile = createTempFile("ccc-test",".tab");
        final File outputFile = createTempFile("ccc-test-out",".hd5");
        writeReadCounts(inputFile, TEST_MATRIX);
        runCommandLine(new String[]{"-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath()});
        checkPCAOutputFile(outputFile);
    }

    private void checkPCAOutputFile(File outputFile) {
        Assert.assertTrue(outputFile.exists());
        try (final HDF5File outHdf5 = new HDF5File(outputFile, HDF5File.OpenMode.READ_ONLY)) {
            final double[][] eigenVectors = outHdf5.readDoubleMatrix(PCA.EIGEN_VECTORS_FULL_PATH);
            final double[] variances = outHdf5.readDoubleArray(PCA.VARIANCES_FULL_PATH);
            final String[] sampleNames = outHdf5.readStringArray(PCA.SAMPLES_FULL_PATH);
            final String[] targetNames = outHdf5.readStringArray(PCA.VARIABLES_FULL_PATH);
            final double[] centers = outHdf5.readDoubleArray(PCA.CENTERS_FULL_PATH);
            Assert.assertEquals(sampleNames.length, PCAUnitTest.TEST_MATRIX[0].length);
            Assert.assertEquals(targetNames.length, PCAUnitTest.TEST_MATRIX.length);
            for (int i = 0; i < sampleNames.length; i++) {
                Assert.assertEquals(sampleNames[i], "SAMPLE_" + i);
            }
            for (int i = 0; i < targetNames.length; i++) {
                Assert.assertEquals(targetNames[i], "TARGET_" + i);
            }
            PCAUnitTest.assertEquals(new ArrayRealVector(centers), new ArrayRealVector(PCAUnitTest.TEST_EXPECTED_CENTERS), PCAUnitTest.EPSILON);
            PCAUnitTest.assertEquals(new ArrayRealVector(variances), new ArrayRealVector(PCAUnitTest.TEST_EXPECTED_VARS), PCAUnitTest.EPSILON);
            PCAUnitTest.assertEqualEigenVectors(new Array2DRowRealMatrix(eigenVectors),
                    new Array2DRowRealMatrix(PCAUnitTest.TEST_EXPECTED_EIGENVECTORS),
                    new ArrayRealVector(variances), PCAUnitTest.EPSILON);
        }
    }

    private void writeReadCounts(final File outputFile, final double[][] testMatrix) {
        final List<String> sampleNames = IntStream.range(0, testMatrix[0].length)
                .mapToObj(i -> "SAMPLE_" + i)
                .collect(Collectors.toList());
        final List<Target> targets = IntStream.range(0, testMatrix.length)
                .mapToObj(i -> new Target("TARGET_" + i, new SimpleInterval("1", i * 100 + 1, i * 100 + 51)))
                .collect(Collectors.toList());

        try (final TableWriter<ReadCountRecord> writer = ReadCountCollectionUtils.writerWithIntervals(new FileWriter(outputFile), sampleNames)) {
            for (int i = 0; i < testMatrix.length; i++) {
                writer.writeRecord(new ReadCountRecord(targets.get(i), testMatrix[i]));
            }
        } catch (final IOException ex) {
            Assert.fail("problems writing in the input coverage file");
        }
    }

}
