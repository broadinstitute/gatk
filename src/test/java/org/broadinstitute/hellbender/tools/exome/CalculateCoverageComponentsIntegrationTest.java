package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.SparkToggleCommandLineProgram;
import org.broadinstitute.hellbender.utils.pca.PCA;
import org.broadinstitute.hellbender.utils.pca.PCAUnitTest;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
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
    public void testSimpleMatrixWithSampleList() throws IOException {
        final File samplesFile = createTempFile("ccc-test",".tab");
        final File inputFile = createTempFile("ccc-test",".tab");
        final File outputFile = createTempFile("ccc-test-out",".hd5");
        writeReadCounts(inputFile, TEST_MATRIX);
        writeSampleNames(samplesFile, IntStream.range(0, TEST_MATRIX[0].length).toArray());
        runCommandLine(new String[]{"-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + SparkToggleCommandLineProgram.DISABLE_SPARK_SHORT_NAME,
                "-" + CalculateCoverageComponents.SAMPLES_FILE_SHORT_NAME, samplesFile.getPath()});

        checkPCAOutputFile(outputFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testSimpleMatrixWithSampleListWithUnknownSample() throws IOException {
        final File samplesFile = createTempFile("ccc-test",".tab");
        final File inputFile = createTempFile("ccc-test",".tab");
        final File outputFile = createTempFile("ccc-test-out",".hd5");
        writeReadCounts(inputFile, TEST_MATRIX);
        writeSampleNames(samplesFile, 0, 1, 2, 3, 4, 5, 6);
        runCommandLine(new String[]{"-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + SparkToggleCommandLineProgram.DISABLE_SPARK_SHORT_NAME,
                "-" + CalculateCoverageComponents.SAMPLES_FILE_SHORT_NAME, samplesFile.getPath()});

        checkPCAOutputFile(outputFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testSimpleMatrixWithEmptySampleList() throws IOException {
        final File samplesFile = createTempFile("ccc-test",".tab");
        final File inputFile = createTempFile("ccc-test",".tab");
        final File outputFile = createTempFile("ccc-test-out",".hd5");
        writeReadCounts(inputFile, TEST_MATRIX);
        writeSampleNames(samplesFile);
        runCommandLine(new String[]{"-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + SparkToggleCommandLineProgram.DISABLE_SPARK_SHORT_NAME,
                "-" + CalculateCoverageComponents.SAMPLES_FILE_SHORT_NAME, samplesFile.getPath()});

        checkPCAOutputFile(outputFile);
    }

    @Test()
    public void testSimpleMatrixExcludingSomeSamples() throws IOException {
        final File samplesFile = createTempFile("ccc-test",".tab");
        final File inputFile = createTempFile("ccc-test",".tab");
        final File outputFile = createTempFile("ccc-test-out",".hd5");
        writeReadCounts(inputFile, TEST_MATRIX);
        writeSampleNames(samplesFile, 0, 2, 3);
        runCommandLine(new String[]{"-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + SparkToggleCommandLineProgram.DISABLE_SPARK_SHORT_NAME,
                "-" + CalculateCoverageComponents.SAMPLES_FILE_SHORT_NAME, samplesFile.getPath()});

        checkPCAOutputFileWithoutSecondSample(outputFile);
    }

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

    /**
     * Runs the {@link CalculateCoverageComponents} tool.
     * @param args the arguments to pass to the tool.
     * @return may be {@code null}.
     */
    public static Object runTool(final String ... args) {
        return new CalculateCoverageComponentsIntegrationTest().runCommandLine(args);
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

    private void checkPCAOutputFile(final File outputFile) {
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

    private void checkPCAOutputFileWithoutSecondSample(final File outputFile) {
        Assert.assertTrue(outputFile.exists());

        final double[][] TEST_EXPECTED_EIGENVECTORS = new double[][] {
                {0.2695519, 0.57468541, 0.7463344},
                {-0.4021812,-0.23574923, 0.4888022},
                {-0.3343424, 0.74282952, -0.3780054},
                {0.7407325, -0.06021966, -0.1599390},
                {0.2768058, 0.21923523, -0.1290706},
                {-0.1688228, 0.10330921, -0.1375852}
        };
        try (final HDF5File outHdf5 = new HDF5File(outputFile, HDF5File.OpenMode.READ_ONLY)) {
            final double[][] eigenVectors = outHdf5.readDoubleMatrix(PCA.EIGEN_VECTORS_FULL_PATH);
            final double[] variances = outHdf5.readDoubleArray(PCA.VARIANCES_FULL_PATH);
            final String[] sampleNames = outHdf5.readStringArray(PCA.SAMPLES_FULL_PATH);
            final String[] targetNames = outHdf5.readStringArray(PCA.VARIABLES_FULL_PATH);
            Assert.assertEquals(sampleNames.length, 3);
            Assert.assertEquals(targetNames.length, 6);
            for (int i = 0; i < sampleNames.length; i++) {
                Assert.assertEquals(sampleNames[i], (i < 1) ? "SAMPLE_" + i : "SAMPLE_" + (i + 1));
            }
            for (int i = 0; i < targetNames.length; i++) {
                Assert.assertEquals(targetNames[i], "TARGET_" + i);
            }
            PCAUnitTest.assertEqualEigenVectors(new Array2DRowRealMatrix(eigenVectors),
                    new Array2DRowRealMatrix(TEST_EXPECTED_EIGENVECTORS),
                    new ArrayRealVector(variances), PCAUnitTest.EPSILON);
        }
    }

    static void writeReadCounts(final File outputFile, final double[][] testMatrix) {
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

    private static void writeSampleNames(final File samplesFile, final int ... sampleIndexes) throws IOException {
        try (final PrintWriter writer = new PrintWriter(new FileWriter(samplesFile))) {
            writer.println(FilterSamples.OUTPUT_NAME_COLUMN);
            for (int i = 0; i < sampleIndexes.length; i++) {
                writer.println("SAMPLE_" + sampleIndexes[i]);
            }
        }
    }
}
