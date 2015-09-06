package org.broadinstitute.hellbender.utils.hdf5;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
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
 * Unit test that cover {@link HDF5Library}, {@link HDF5File} and {@link HDF5PoN}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HDF5LibraryUnitTests {

    private static File TEST_RESOURCE_DIR = new File("src/test/resources/org/broadinstitute/utils/hdf5");
    private static File TEST_PON_TARGETS = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-targets.txt");
    private static File TEST_PON_SAMPLES = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-samples.txt");
    private static File TEST_PON_TARGET_FACTORS = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-target_factors.txt");
    private static File TEST_PON_NORMALIZED_PCOV = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-normalized_pcov.txt");
    private static File TEST_PON_LOG_NORMAL_SAMPLES = TEST_PON_SAMPLES;
    private static File TEST_PON_LOG_NORMALS = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-log_normals.txt");
    private static File TEST_PON_LOG_NORMALS_PINV = new File(TEST_RESOURCE_DIR,"test_creation_of_panel-log_normal_pinv.txt");
    private static File TEST_PON_REDUCED_PON = TEST_PON_LOG_NORMALS;
    private static File TEST_PON_REDUCED_PON_PINV = TEST_PON_LOG_NORMALS_PINV;
    public static File TEST_PON = new File(TEST_RESOURCE_DIR,"test_creation_of_panel.pon");

    @SuppressWarnings("FieldCanBeLocal")
    private static double TEST_PON_VERSION = 1.3;

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
        final HDF5File reader = new HDF5File(TEST_PON);
        reader.close();
    }

    @Test(dependsOnGroups = "supported",expectedExceptions = GATKException.class)
    public void testOpenReadOnlyOnBadFile() {
        final HDF5File reader = new HDF5File(new File("/tmp/no-file"));
        reader.close();
    }

    @Test(dependsOnGroups = "supported")
    public void testTargetNameReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targetNames = pon.getTargetNames();
        final List<String> expected = readLines(TEST_PON_TARGETS);
        Assert.assertNotNull(targetNames);
        Assert.assertEquals(targetNames,expected);
        reader.close();
    }

    @Test(dependsOnGroups = "supported")
    public void testSampleNameReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> sampleNames = pon.getSampleNames();
        final List<String> expected = readLines(TEST_PON_SAMPLES);
        Assert.assertNotNull(sampleNames);
        Assert.assertEquals(sampleNames, expected);
        reader.close();
    }

    @Test(dependsOnMethods = "testSampleNameReading")
    public void testLogNormalSampleNameReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PoN pon = new HDF5PoN(reader);

        final List<String> expected = readLines(TEST_PON_LOG_NORMAL_SAMPLES);

        final List<String> logNormalSampleNames = pon.getPanelSampleNames().stream().sorted().collect(Collectors.toList());
        Assert.assertEquals(logNormalSampleNames, expected);
    }

    @Test(dependsOnGroups = "supported", dependsOnMethods = "testTargetNameReading")
    public void testTargetFactorsReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.getTargetNames();
        final RealMatrix factors = pon.getTargetFactors();
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
    public void testVersionReading() {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        Assert.assertEquals(pon.getVersion(),TEST_PON_VERSION);
        reader.close();
    }

    @Test(dependsOnGroups = "supported", dependsOnMethods = {"testTargetNameReading","testSampleNameReading"})
    public void testNormalizedPcovReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.getTargetNames();
        final List<String> samples = pon.getSampleNames();
        final RealMatrix actual = pon.getNormalCounts();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), targets.size());
        Assert.assertEquals(actual.getColumnDimension(), samples.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_NORMALIZED_PCOV);
        assertEqualMatrix(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalSampleNameReading"})
    public void testLogNormalMatrixReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.getTargetNames();
        final List<String> samples = pon.getPanelSampleNames();
        final RealMatrix actual = pon.getLogNormalCounts();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), targets.size());
        Assert.assertEquals(actual.getColumnDimension(), samples.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_LOG_NORMALS);
        assertEqualMatrix(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalSampleNameReading"})
    public void testLogNormalPInvMatrixReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.getTargetNames();
        final List<String> samples = pon.getPanelSampleNames();
        final RealMatrix actual = pon.getLogNormalPInverseCounts();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), samples.size());
        Assert.assertEquals(actual.getColumnDimension(), targets.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_LOG_NORMALS_PINV);
        assertEqualMatrix(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalSampleNameReading"})
    public void testReducedPoNMatrixReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.getTargetNames();
        final List<String> samples = pon.getPanelSampleNames();
        final RealMatrix actual = pon.getReducedPanelCounts();
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getRowDimension(), targets.size());
        Assert.assertTrue(actual.getColumnDimension() <= samples.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_REDUCED_PON);
        assertEqualMatrix(actual, expected);
    }

    @Test(dependsOnMethods = {"testTargetNameReading","testLogNormalSampleNameReading"})
    public void testReducedPoNPInvMatrixReading() throws IOException {
        final HDF5File reader = new HDF5File(TEST_PON);
        final PoN pon = new HDF5PoN(reader);
        final List<String> targets = pon.getTargetNames();
        final List<String> samples = pon.getPanelSampleNames();
        final RealMatrix actual = pon.getReducedPanelPInverseCounts();
        Assert.assertNotNull(actual);
        Assert.assertTrue(actual.getRowDimension() <= samples.size());
        Assert.assertEquals(actual.getColumnDimension(), targets.size());
        final RealMatrix expected = readDoubleMatrix(TEST_PON_REDUCED_PON_PINV);
        assertEqualMatrix(actual, expected);
    }

    @Test()
    public void testCreateHDF5File() throws IOException {
        final File testFile = BaseTest.createTempFile("hdf5", ".hd5");
        testFile.delete();
        final HDF5File file = new HDF5File(testFile, HDF5File.OpenMode.CREATE);
        file.close();
    }

    @Test()
    public void testCreateGroup() throws IOException {
        final File testFile = BaseTest.createTempFile("hdf5", ".hd5");
        final HDF5File file = new HDF5File(testFile, HDF5File.OpenMode.CREATE);
        Assert.assertTrue(file.makeGroup("test-group/lola-run"));

        Assert.assertFalse(file.makeGroup("test-group"));
        Assert.assertFalse(file.makeGroup("test-group/lola-run"));
        Assert.assertTrue(file.makeGroup("test-group/peter-pan"));
        file.close();
    }

    @Test()
    public void testMakeDouble() throws IOException {
        final File testFile = File.createTempFile("hdf5", ".hd5");
        HDF5File file = new HDF5File(testFile, HDF5File.OpenMode.CREATE);
        file.makeGroup("test-group/double-group");
        Assert.assertTrue(file.makeDouble("test-group/double-group/my-double", 1.1));
        System.err.println(testFile);
        file.close();
        final long time = System.currentTimeMillis();
        Assert.assertTrue(testFile.length() > 0);
        Assert.assertTrue(testFile.lastModified() <= time);
        file = new HDF5File(testFile, HDF5File.OpenMode.READ_ONLY);
        final double theDouble = file.readDouble("test-group/double-group/my-double");
        Assert.assertEquals(theDouble, 1.1);
        file.close();
    }

    @Test()
    public void testMakeNaNDouble() throws IOException {
        final File testFile = File.createTempFile("hdf5", ".hd5");
        HDF5File file = new HDF5File(testFile, HDF5File.OpenMode.CREATE);
        file.makeGroup("test-group/double-group");
        Assert.assertTrue(file.makeDouble("test-group/double-group/my-double", Double.NaN));
        System.err.println(testFile);
        file.close();
        final long time = System.currentTimeMillis();
        Assert.assertTrue(testFile.length() > 0);
        Assert.assertTrue(testFile.lastModified() <= time);
        file = new HDF5File(testFile, HDF5File.OpenMode.READ_ONLY);
        final double theDouble = file.readDouble("test-group/double-group/my-double");
        Assert.assertTrue(Double.isNaN(theDouble));
        file.close();
    }

    @Test()
    public void testMakeDoubleArray() throws IOException {
        final File testFile = File.createTempFile("hdf5", ".hd5");
        HDF5File file = new HDF5File(testFile, HDF5File.OpenMode.CREATE);
        file.makeGroup("test-group/double-group");
        final double[] testValues = new double[] { 1.1 , -2.2, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 0.0111e10-10 };
        Assert.assertTrue(file.makeDoubleArray("test-group/double-group/my-double", testValues));
        System.err.println(testFile);
        file.close();
        final long time = System.currentTimeMillis();
        Assert.assertTrue(testFile.length() > 0);
        Assert.assertTrue(testFile.lastModified() <= time);
        file = new HDF5File(testFile, HDF5File.OpenMode.READ_ONLY);
        final double[] theDoubles = file.readDoubleArray("test-group/double-group/my-double");
        Assert.assertEquals(theDoubles, testValues.clone());
        file.close();
    }

    @Test()
    public void testReMakeDoubleArray() throws IOException {
        final File testFile = File.createTempFile("hdf5", ".hd5");
        HDF5File file = new HDF5File(testFile, HDF5File.OpenMode.CREATE);
        file.makeGroup("test-group/double-group");
        final double[] testValues1 = new double[] { 1.1 , -2.2, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 0.0111e10-10 };
        final double[] testValues2 = new double[] { 11.1 , -22.2, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1.0111e10-10 };

        Assert.assertTrue(file.makeDoubleArray("test-group/double-group/my-double", testValues1));
        System.err.println(testFile);
        file.close();
        final long time = System.currentTimeMillis();
        Assert.assertTrue(testFile.length() > 0);
        Assert.assertTrue(testFile.lastModified() <= time);
        file = new HDF5File(testFile, HDF5File.OpenMode.READ_WRITE);
        final double[] theDoubles1 = file.readDoubleArray("test-group/double-group/my-double");
        Assert.assertEquals(theDoubles1, testValues1.clone());
        Assert.assertFalse(file.makeDoubleArray("test-group/double-group/my-double", testValues2));
        final double[] theDoubles2 = file.readDoubleArray("test-group/double-group/my-double");
        Assert.assertEquals(theDoubles2, testValues2.clone());

        file.close();
    }

    @Test()
    public void testMakeDoubleMatrix() throws IOException {
        final File testFile = File.createTempFile("hdf5", ".hd5");
        HDF5File file = new HDF5File(testFile, HDF5File.OpenMode.CREATE);
        file.makeGroup("test-group/double-group");
        final double[][] testValues = new double[][] {
                new double[] { 1.1 , -2.2, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 0.0111e10-10 },
                new double[] { -1.1, 2.2, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, -0.01111e10-10 }};
        Assert.assertTrue(file.makeDoubleMatrix("test-group/double-group/my-double", testValues));
        System.err.println(testFile);
        file.close();
        final long time = System.currentTimeMillis();
        Assert.assertTrue(testFile.length() > 0);
        Assert.assertTrue(testFile.lastModified() <= time);
        file = new HDF5File(testFile, HDF5File.OpenMode.READ_ONLY);

        final double[][] theDoubles = file.readDoubleMatrix("test-group/double-group/my-double");
        Assert.assertEquals(theDoubles, testValues.clone());
        file.close();
    }

    @Test()
    public void testMakeStringArray() throws IOException {
        final File testFile = File.createTempFile("hdf5", ".hd5");
        HDF5File file = new HDF5File(testFile, HDF5File.OpenMode.CREATE);
        file.makeGroup("test-group/double-group");
        final String[] testValues = new String[] { "0", "1", "absdsd12 sdsad121 sdasadsad 1212sdasdas",
                StringUtils.repeat("x", 2000) };
        Assert.assertTrue(file.makeStringArray("test-group/double-group/my-double", testValues));
        System.err.println(testFile);
        file.close();
        FileUtils.copyFile(testFile, new File("/tmp/3.hd5"));
        final long time = System.currentTimeMillis();
        Assert.assertTrue(testFile.length() > 0);
        Assert.assertTrue(testFile.lastModified() <= time);
        file = new HDF5File(testFile, HDF5File.OpenMode.READ_ONLY);
        final String[] theStrings = file.readStringArray("test-group/double-group/my-double");
        Assert.assertEquals(theStrings, testValues.clone());
        file.close();
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
