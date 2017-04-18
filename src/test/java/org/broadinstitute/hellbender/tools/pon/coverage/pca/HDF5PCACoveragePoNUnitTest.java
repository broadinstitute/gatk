package org.broadinstitute.hellbender.tools.pon.coverage.pca;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathObjectAsserts;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;


/**
 * Unit tests that cover {@link HDF5PCACoveragePoN}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HDF5PCACoveragePoNUnitTest extends BaseTest {
    private static final File TEST_RESOURCE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private static final File TEST_PON_TARGETS = new File(TEST_RESOURCE_DIR, "test_creation_of_panel-targets.txt");
    private static final File TEST_PON_SAMPLES = new File(TEST_RESOURCE_DIR, "test_creation_of_panel-samples.txt");
    private static final File TEST_PON_TARGET_FACTORS = new File(TEST_RESOURCE_DIR, "test_creation_of_panel-target_factors.txt");
    private static final File TEST_PON_NORMALIZED_PCOV = new File(TEST_RESOURCE_DIR, "test_creation_of_panel-normalized_pcov.txt");
    private static final File TEST_PON_LOG_NORMAL_SAMPLES = TEST_PON_SAMPLES;
    private static final File TEST_PON_LOG_NORMALS = new File(TEST_RESOURCE_DIR, "test_creation_of_panel-log_normals.txt");
    private static final File TEST_PON_LOG_NORMALS_PINV = new File(TEST_RESOURCE_DIR, "test_creation_of_panel-log_normal_pinv.txt");
    private static final File TEST_PON_REDUCED_PON = TEST_PON_LOG_NORMALS;
    private static final File TEST_PON_REDUCED_PON_PINV = TEST_PON_LOG_NORMALS_PINV;
    private static final File TEST_PON = new File(TEST_RESOURCE_DIR, "test_creation_of_panel.pon");

    private static final double TEST_PON_VERSION = 1.3;

    @BeforeTest
    void loadHDF5() throws UserException {
        if (!new HDF5Library().load(null)) {  //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }
    }

    @Test
    public void testTargetNameReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PCACoveragePoN pon = new HDF5PCACoveragePoN(reader);
        final List<String> targetNames = pon.getTargetNames();
        final List<String> expected = readLines(TEST_PON_TARGETS);
        Assert.assertNotNull(targetNames);
        Assert.assertEquals(targetNames,expected);
        reader.close();
    }

    @Test
    public void testSampleNameReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PCACoveragePoN pon = new HDF5PCACoveragePoN(reader);
        final List<String> sampleNames = pon.getSampleNames();
        final List<String> expected = readLines(TEST_PON_SAMPLES);
        Assert.assertNotNull(sampleNames);
        Assert.assertEquals(sampleNames, expected);
        reader.close();
    }

    @Test(dependsOnMethods = "testSampleNameReading")
    public void testLogNormalizedSampleNameReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PCACoveragePoN pon = new HDF5PCACoveragePoN(reader);

        final List<String> expected = readLines(TEST_PON_LOG_NORMAL_SAMPLES);

        final List<String> logNormalizedSampleNames = pon.getPanelSampleNames().stream().sorted().collect(Collectors.toList());
        Assert.assertEquals(logNormalizedSampleNames, expected);
    }

    @Test(dependsOnMethods = "testTargetNameReading")
    public void testTargetFactorsReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PCACoveragePoN pon = new HDF5PCACoveragePoN(reader);
        final List<String> targets = pon.getTargetNames();
        final double[] factors = pon.getTargetFactors();
        Assert.assertEquals(factors.length, targets.size());

        final List<Double> expected = readDoubleLines(TEST_PON_TARGET_FACTORS);
        Assert.assertNotNull(factors);
        Assert.assertEquals(factors.length, expected.size());
        for (int i = 0; i < expected.size(); i++) {
            Assert.assertEquals(factors[i], expected.get(i), Math.abs(expected.get(i)) * 0.0001);
        }
        reader.close();
    }

    @Test
    public void testVersionReading() {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PCACoveragePoN pon = new HDF5PCACoveragePoN(reader);
        Assert.assertEquals(pon.getVersion(),TEST_PON_VERSION);
        reader.close();
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testSampleNameReading"})
    public void testNormalizedPcovReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PCACoveragePoN pon = new HDF5PCACoveragePoN(reader);
        final List<String> targets = pon.getTargetNames();
        final List<String> samples = pon.getSampleNames();
        final RealMatrix actual = pon.getNormalizedCounts();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), targets.size());
        Assert.assertEquals(actual.getColumnDimension(), samples.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_NORMALIZED_PCOV);
        MathObjectAsserts.assertRealMatrixEquals(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalizedSampleNameReading"})
    public void testLogNormalizedMatrixReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PCACoveragePoN pon = new HDF5PCACoveragePoN(reader);
        final List<String> targets = pon.getTargetNames();
        final List<String> samples = pon.getPanelSampleNames();
        final RealMatrix actual = pon.getLogNormalizedCounts();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), targets.size());
        Assert.assertEquals(actual.getColumnDimension(), samples.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_LOG_NORMALS);
        MathObjectAsserts.assertRealMatrixEquals(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalizedSampleNameReading"})
    public void testLogNormalizedPInvMatrixReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PCACoveragePoN pon = new HDF5PCACoveragePoN(reader);
        final List<String> targets = pon.getTargetNames();
        final List<String> samples = pon.getPanelSampleNames();
        final RealMatrix actual = pon.getLogNormalizedPInverseCounts();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), samples.size());
        Assert.assertEquals(actual.getColumnDimension(), targets.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_LOG_NORMALS_PINV);
        MathObjectAsserts.assertRealMatrixEquals(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalizedSampleNameReading"})
    public void testReducedPoNMatrixReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PCACoveragePoN pon = new HDF5PCACoveragePoN(reader);
        final List<String> targets = pon.getTargetNames();
        final List<String> samples = pon.getPanelSampleNames();
        final RealMatrix actual = pon.getReducedPanelCounts();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), targets.size());
        Assert.assertTrue(actual.getColumnDimension() <= samples.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_REDUCED_PON);
        MathObjectAsserts.assertRealMatrixEquals(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalizedSampleNameReading"})
    public void testReducedPoNPInvMatrixReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PCACoveragePoN pon = new HDF5PCACoveragePoN(reader);
        final List<String> targets = pon.getTargetNames();
        final List<String> samples = pon.getPanelSampleNames();
        final RealMatrix actual = pon.getReducedPanelPInverseCounts();
        Assert.assertNotNull(actual);
        Assert.assertTrue(actual.getRowDimension() <= samples.size());
        Assert.assertEquals(actual.getColumnDimension(), targets.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_REDUCED_PON_PINV);
        MathObjectAsserts.assertRealMatrixEquals(actual, expected);
    }

    /**
     * Reads the lines of a file into an string array.
     * @param file the input file.
     * @return never {@code null}, a list with as many elements as lines in {@code file}.
     * @throws IOException if some IO exception occurred.
     */
    private List<String> readLines(final File file) throws IOException {
        final BufferedReader reader = new BufferedReader(new FileReader(file));
        final ArrayList<String> result = new ArrayList<>(1000);
        String line;
        while ((line = reader.readLine()) != null) {
            result.add(line);
        }
        return result;
    }

    /**
     * Reads the lines of a file into a double list.
     * @param file the input file.
     * @return never {@code null}, a list with as many elements as lines in {@code file}.
     * @throws IOException if some IO exception occurred.
     * @throws NumberFormatException if some of the input lines in {@code file} cannot be converted into a {@code double}.
     */
    private List<Double> readDoubleLines(final File file) throws IOException {
        return readLines(file).stream().map(Double::valueOf).collect(Collectors.toList());
    }

    /**
     * Reads the lines of a file into a double 2-D matrix.
     * <p>
     *     Each line in the file correspond to a row in the result matrix and values are separated by blank
     *     characters (such as spaces or tabs).
     * </p>
     * @param file the input file.
     * @return never {@code null}, a matrix with as many rows as lines in the input file, and as many columns are
     *   values per input file line.
     * @throws IOException if there is some IO problem accessing {@code file}.
     * @throws NullPointerException if some value in {@code file} could not be converted to a {@code double}.
     * @throws org.apache.commons.math3.exception.NoDataException if {@code file} does not contain any data.
     * @throws org.apache.commons.math3.exception.DimensionMismatchException if {@code file} lines don't have
     *    the same length (number of columns.
     */
    private RealMatrix readDoubleMatrix(final File file) throws IOException {
        final List<String> lines = readLines(file);
        final double[][] values = lines.stream()
                .map(line -> Stream.of(line.split("\\s+")).mapToDouble(Double::valueOf).toArray())
                .collect(Collectors.toList())
                .toArray(new double[lines.size()][]);
        return new Array2DRowRealMatrix(values,false);
    }
}