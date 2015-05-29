package org.broadinstitute.hellbender.utils.hdf5;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;
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
 * Unit test that cover {@link HDF5Library}, {@link HDF5Reader} and {@link HDF5PoN}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HDF5LibraryUnitTests {


    private static File TEST_RESOURCE_DIR = new File("src/test/resources/org/broadinstitute/utils/hdf5");
    public static File TEST_PON = new File(TEST_RESOURCE_DIR,"test_creation_of_panel.pon");
    private static File TEST_PON_TARGETS = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-targets.txt");
    private static File TEST_PON_SAMPLES = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-samples.txt");
    private static File TEST_PON_TARGET_FACTORS = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-target_factors.txt");
    private static File TEST_PON_NORMALIZED_PCOV = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-normalized_pcov.txt");
    private static File TEST_PON_LOG_NORMAL_SAMPLES = TEST_PON_SAMPLES;
    private static File TEST_PON_LOG_NORMALS = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-log_normals.txt");
    private static File TEST_PON_LOG_NORMALS_PINV = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-log_normal_pinv.txt");
    private static File TEST_PON_REDUCED_PON = TEST_PON_LOG_NORMALS;
    private static File TEST_PON_REDUCED_PON_PINV = TEST_PON_LOG_NORMALS_PINV;

    @SuppressWarnings("FieldCanBeLocal")
    private static double TEST_PON_MAXIMUM_RATIO_CUTOFF = 0.0;
    @SuppressWarnings("FieldCanBeLocal")
    private static String TEST_PON_VERSION = "1.3";

    @Test(groups = "supported")
    public void testIsSupported() {
        Assert.assertTrue(HDF5Library.isSupported());
    }


    @Test(dependsOnGroups = "supported")
    public void testGetLibrary() {
        final HDF5Library library = HDF5Library.getLibrary();
        Assert.assertNotNull(library);
        Assert.assertSame(HDF5Library.getLibrary(), library);
    }

    @Test(dependsOnGroups = "supported")
    public void testOpenReadOnly() {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        reader.close();
    }

    @Test(dependsOnGroups = "supported",expectedExceptions = GATKException.class)
    public void testOpenReadOnlyOnBadFile() {
        final HDF5Reader reader = new HDF5Reader(new File("/tmp/no-file"));
        reader.close();
    }

    @Test(dependsOnGroups = "supported")
    public void testTargetNameReading() throws IOException {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targetNames = pon.targetNames();
        final List<String> expected = readLines(TEST_PON_TARGETS);
        Assert.assertNotNull(targetNames);
        Assert.assertEquals(targetNames,expected);
        reader.close();
    }

    @Test(dependsOnGroups = "supported")
    public void testSampleNameReading() throws IOException {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> sampleNames = pon.sampleNames();
        final List<String> expected = readLines(TEST_PON_SAMPLES);
        Assert.assertNotNull(sampleNames);
        Assert.assertEquals(sampleNames, expected);
        reader.close();
    }

    @Test(dependsOnMethods = "testSampleNameReading")
    public void testLogNormalSampleNameReading() throws IOException {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        final PoN pon = new HDF5PoN(reader);

        final List<String> expected = readLines(TEST_PON_LOG_NORMAL_SAMPLES);

        final List<String> logNormalSampleNames = pon.logNormalSampleNames().stream().sorted().collect(Collectors.toList());
        Assert.assertEquals(logNormalSampleNames, expected);
    }

    @Test(dependsOnGroups = "supported", dependsOnMethods = "testTargetNameReading")
    public void testTargetFactorsReading() throws IOException {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.targetNames();
        final RealMatrix factors = pon.targetFactors();
        Assert.assertEquals(factors.getRowDimension(),targets.size());
        Assert.assertEquals(factors.getColumnDimension(),1);

        final List<Double> expected = readDoubleLines(TEST_PON_TARGET_FACTORS);
        Assert.assertNotNull(factors);
        Assert.assertEquals(factors.getRowDimension(),expected.size());
        for (int i = 0; i < expected.size(); i++) {
            Assert.assertEquals(factors.getEntry(i,0),expected.get(i),Math.abs(expected.get(i)) * 0.0001);
        }
        reader.close();
    }

    @Test(dependsOnGroups = "supported")
    public void testMaximumRatioCutoffReading() {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        Assert.assertEquals(pon.maximumRatioCutoff(),TEST_PON_MAXIMUM_RATIO_CUTOFF,0.0001);
        reader.close();
    }

    @Test(dependsOnGroups = "supported")
    public void testVersionReading() {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        Assert.assertEquals(pon.version(),TEST_PON_VERSION);
        reader.close();
    }

    @Test(dependsOnGroups = "supported", dependsOnMethods = {"testTargetNameReading","testSampleNameReading"})
    public void testNormalizedPcovReading() throws IOException {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.targetNames();
        final List<String> samples = pon.sampleNames();
        final RealMatrix actual = pon.normalizedPercentCoverage();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), targets.size());
        Assert.assertEquals(actual.getColumnDimension(), samples.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_NORMALIZED_PCOV);
        assertEqualMatrix(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalSampleNameReading"})
    public void testLogNormalMatrixReading() throws IOException {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.targetNames();
        final List<String> samples = pon.logNormalSampleNames();
        final RealMatrix actual = pon.logNormals();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), targets.size());
        Assert.assertEquals(actual.getColumnDimension(), samples.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_LOG_NORMALS);
        assertEqualMatrix(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalSampleNameReading"})
    public void testLogNormalPInvMatrixReading() throws IOException {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.targetNames();
        final List<String> samples = pon.logNormalSampleNames();
        final RealMatrix actual = pon.logNormalsPseudoInverse();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), samples.size());
        Assert.assertEquals(actual.getColumnDimension(), targets.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_LOG_NORMALS_PINV);
        assertEqualMatrix(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalSampleNameReading"})
    public void testReducedPoNMatrixReading() throws IOException {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.targetNames();
        final List<String> samples = pon.logNormalSampleNames();
        final RealMatrix actual = pon.reducedPoN();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), targets.size());
        Assert.assertTrue(actual.getColumnDimension() <= samples.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_REDUCED_PON);
        assertEqualMatrix(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalSampleNameReading"})
    public void testReducedPoNPInvMatrixReading() throws IOException {
        final HDF5Reader reader = new HDF5Reader(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.targetNames();
        final List<String> samples = pon.logNormalSampleNames();
        final RealMatrix actual = pon.reducedPoNPseudoInverse();
        Assert.assertNotNull(actual);
        Assert.assertTrue(actual.getRowDimension() <= samples.size());
        Assert.assertEquals(actual.getColumnDimension(), targets.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_REDUCED_PON_PINV);
        assertEqualMatrix(actual, expected);
    }

    private void assertEqualMatrix(RealMatrix actual, RealMatrix expected) {
        Assert.assertEquals(actual.getRowDimension(), expected.getRowDimension());
        Assert.assertEquals(actual.getColumnDimension(), expected.getColumnDimension());
        for (int row = 0; row < expected.getRowDimension(); row++) {
            for (int column = 0; column < expected.getColumnDimension(); column++) {
                final double actualValue = actual.getEntry(row,column);
                final double expectedValue = expected.getEntry(row,column);
                final double epsilon = Math.min(Math.abs(actualValue),Math.abs(expectedValue)) * 0.0001;
                Assert.assertEquals(actualValue,expectedValue,epsilon);
            }
        }
    }

    /**
     * Reads the lines of a file into an string array.
     * @param file the input file.
     * @return never {@code null}, a list with as many elements as lines in {@code file}.
     * @throws IOException if some IO exception occurred.
     */
    private List<String> readLines(final File file) throws IOException {

        final BufferedReader  reader = new BufferedReader(new FileReader(file));
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
