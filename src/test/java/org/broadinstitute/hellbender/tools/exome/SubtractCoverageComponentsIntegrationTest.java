package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.SparkToggleCommandLineProgram;
import org.broadinstitute.hellbender.utils.pca.PCA;
import org.broadinstitute.hellbender.utils.pca.PCAUnitTest;
import org.junit.AfterClass;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Integration test for {@link SubtractCoverageComponents}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class SubtractCoverageComponentsIntegrationTest extends CommandLineProgramTest {

    private static final double[][] TEST_MATRIX = PCAUnitTest.TEST_MATRIX;

    private static File TEST_PCA_INPUT_FILE;
    private static File TEST_PCA_FILE;

    @BeforeClass
    protected static void composeTestFiles() {
        TEST_PCA_INPUT_FILE = createTempFile("test-pca-input-file", ".tab");
        TEST_PCA_FILE = createTempFile("test-pca-file", ".hd5");
        CalculateCoverageComponentsIntegrationTest.writeReadCounts(TEST_PCA_INPUT_FILE, TEST_MATRIX);
        // Run the PCA tool to generate the pca file used in the tests.
        CalculateCoverageComponentsIntegrationTest.runTool("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, TEST_PCA_INPUT_FILE.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, TEST_PCA_FILE.getPath(),
                "-" + SparkToggleCommandLineProgram.DISABLE_SPARK_SHORT_NAME);
    }

    @Override
    public String getTestedClassName() {
        return SubtractCoverageComponents.class.getSimpleName();
    }

    @Test
    public void testRun() throws IOException {

        final File outputFile = createTempFile("out-file", ".tab");
        outputFile.delete();

        final String[] arguments = new String[] {
                "-" + SubtractCoverageComponents.PCA_INPUT_SHORT_NAME, TEST_PCA_FILE.getPath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, TEST_PCA_INPUT_FILE.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
        };
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.canWrite());
        Assert.assertTrue(outputFile.isFile());
        assertCanObtainOriginalValuesByReversingSubtraction(outputFile, TEST_PCA_INPUT_FILE, Integer.MAX_VALUE, 1.0);
        outputFile.delete();
    }

    @Test
    public void testNumberOfComponents() throws IOException {
        for (final int numComponentsToUse : new int[] {0, 1, 5, 10000}) {
            final File outputFile = createTempFile("out-file", ".tab");
            outputFile.delete();

            final String[] arguments = new String[]{
                    "-" + SubtractCoverageComponents.PCA_INPUT_SHORT_NAME, TEST_PCA_FILE.getPath(),
                    "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, TEST_PCA_INPUT_FILE.getPath(),
                    "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                    "-" + SubtractCoverageComponents.NUM_COMPONENTS_SHORT_NAME, Integer.toString(numComponentsToUse),
            };
            runCommandLine(arguments);
            Assert.assertTrue(outputFile.exists());
            Assert.assertTrue(outputFile.canWrite());
            Assert.assertTrue(outputFile.isFile());
            assertCanObtainOriginalValuesByReversingSubtraction(outputFile, TEST_PCA_INPUT_FILE, numComponentsToUse, 1.0);
            outputFile.delete();
        }
    }

    @Test
    public void testProportionOfVariance() throws IOException {
        for (final double proportionOfVariance : new double[] {0, 0.5, 1.0}) {
            final File outputFile = createTempFile("out-file", ".tab");
            outputFile.delete();

            final String[] arguments = new String[]{
                    "-" + SubtractCoverageComponents.PCA_INPUT_SHORT_NAME, TEST_PCA_FILE.getPath(),
                    "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, TEST_PCA_INPUT_FILE.getPath(),
                    "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                    "-" + SubtractCoverageComponents.PROPORTION_OF_VARIANCE_SHORT_NAME, Double.toString(proportionOfVariance),
            };
            runCommandLine(arguments);
            Assert.assertTrue(outputFile.exists());
            Assert.assertTrue(outputFile.canWrite());
            Assert.assertTrue(outputFile.isFile());
            assertCanObtainOriginalValuesByReversingSubtraction(outputFile, TEST_PCA_INPUT_FILE, Integer.MAX_VALUE, proportionOfVariance);
            outputFile.delete();
        }
    }

    @Test
    public void testNumberOfComponentsAndProportionOfVariance() throws IOException {
        for (final int numComponentsToUse : new int[] {0, 1, 5, 10000}) {
            for (final double proportionOfVariance : new double[]{0, 0.5, 1.0}) {
                final File outputFile = createTempFile("out-file", ".tab");
                outputFile.delete();

                final String[] arguments = new String[]{
                        "-" + SubtractCoverageComponents.PCA_INPUT_SHORT_NAME, TEST_PCA_FILE.getPath(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, TEST_PCA_INPUT_FILE.getPath(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                        "-" + SubtractCoverageComponents.NUM_COMPONENTS_SHORT_NAME, Integer.toString(numComponentsToUse),
                        "-" + SubtractCoverageComponents.PROPORTION_OF_VARIANCE_SHORT_NAME, Double.toString(proportionOfVariance),
                };
                runCommandLine(arguments);
                Assert.assertTrue(outputFile.exists());
                Assert.assertTrue(outputFile.canWrite());
                Assert.assertTrue(outputFile.isFile());
                assertCanObtainOriginalValuesByReversingSubtraction(outputFile, TEST_PCA_INPUT_FILE, numComponentsToUse, proportionOfVariance);
                outputFile.delete();
            }
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativeComponents() throws Exception {
        final File outputFile = createTempFile("out-file", ".tab");
        outputFile.delete();

        final String[] arguments = new String[]{
                "-" + SubtractCoverageComponents.PCA_INPUT_SHORT_NAME, TEST_PCA_FILE.getPath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, TEST_PCA_INPUT_FILE.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + SubtractCoverageComponents.NUM_COMPONENTS_SHORT_NAME, "-1",
        };
        runCommandLine(arguments);
        outputFile.delete();
    }

    private void assertCanObtainOriginalValuesByReversingSubtraction(final File outputFile, final File inputFile,
            final int numComponentsRequested, final double proportionOfVariance) throws IOException {

        final PCA pca;
        try (final HDF5File pcaFile = new HDF5File(TEST_PCA_FILE, HDF5File.OpenMode.READ_ONLY)) {
            pca = PCA.readHDF5(pcaFile);
        }


        final int numComponentsUsed = Math.min(numComponentsRequested,
                pca.numComponentsToAccountForVariance(proportionOfVariance));
        final RealVector centers = pca.getCenters();
        final List<String> variables = pca.getVariables();
        final List<Target> targets = variables.stream().map(Target::new).collect(Collectors.toList());
        final ReadCountCollection inputCounts = ReadCountCollectionUtils.parse(inputFile).arrangeTargets(targets);
        final ReadCountCollection outputCounts = ReadCountCollectionUtils.parse(outputFile).arrangeTargets(targets);
        final List<String> samples = inputCounts.columnNames();
        final RealMatrix eigenVectors = pca.getEigenVectors();

        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final int outputIndex = outputCounts.columnNames().indexOf(sample);
            if (outputIndex < 0) {
                Assert.fail("cannot find sample: " + sample + " in " + Arrays.toString(outputCounts.columnNames().toArray()));
            }
            final double[] inputSampleValues = inputCounts.counts().getColumn(i);
            final double[] centeredSampleValues = new double[inputSampleValues.length];
            for (int j = 0; j < centeredSampleValues.length; j++) {
                centeredSampleValues[j] = inputSampleValues[j] - centers.getEntry(j);
            }
            final double[] outputSampleValues = outputCounts.counts().getColumn(outputIndex);
            final double[] expectedSampleValues = centeredSampleValues.clone();
            final double[] revertedSampleValues = outputSampleValues.clone();
            for (int k = 0; k < Math.min(eigenVectors.getColumnDimension(), numComponentsUsed); k++) {
                final double loading = new ArrayRealVector(eigenVectors.getColumn(k)).dotProduct(new ArrayRealVector(centeredSampleValues));
                final RealVector projection = new ArrayRealVector(eigenVectors.getColumn(k)).mapMultiply(loading);
                final RealVector updatedCenteredSampleValues = new ArrayRealVector(centeredSampleValues).subtract(projection);
                System.arraycopy(updatedCenteredSampleValues.toArray(), 0, centeredSampleValues, 0, centeredSampleValues.length);
                final RealVector updatedRevertedValues = new ArrayRealVector(revertedSampleValues).add(projection);
                System.arraycopy(updatedRevertedValues.toArray(), 0, revertedSampleValues, 0, revertedSampleValues.length);
            }
            for (int k = 0; k < revertedSampleValues.length; k++) {
                Assert.assertEquals(revertedSampleValues[k], expectedSampleValues[k], 0.000001, " i and k  = " + i + " " + k);
            }
        }
    }

    @AfterClass
    protected static void disposeOfTestFile() {
        TEST_PCA_INPUT_FILE.delete();
        TEST_PCA_INPUT_FILE = null;
        TEST_PCA_FILE.delete();
        TEST_PCA_FILE = null;
    }
}
