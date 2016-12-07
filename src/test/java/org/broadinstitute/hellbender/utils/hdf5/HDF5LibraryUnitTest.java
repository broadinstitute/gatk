package org.broadinstitute.hellbender.utils.hdf5;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.tools.pon.PoNTestUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Unit tests that cover {@link HDF5Library} and {@link HDF5File}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HDF5LibraryUnitTest {
    private static final File TEST_RESOURCE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private static final File TEST_PON = new File(TEST_RESOURCE_DIR, "test_creation_of_panel.pon");

    @BeforeClass
    public void testIsSupported() {
        Assert.assertTrue(new HDF5Library().load(null));
    }

    @Test
    public void testOpenReadOnly() {
        final HDF5File reader = new HDF5File(TEST_PON);
        reader.close();
    }

    @Test(expectedExceptions = HDF5LibException.class)
    public void testOpenReadOnlyOnBadFile() {
        final HDF5File reader = new HDF5File(new File("/tmp/no-file"));
        reader.close();
    }

    @Test()
    public void testCreateHDF5File() {
        final File testFile = BaseTest.createTempFile("hdf5", ".hd5");
        testFile.delete();
        final HDF5File file = new HDF5File(testFile, HDF5File.OpenMode.CREATE);
        file.close();
    }

    @Test(dependsOnMethods = {"testCreateGroup", "testMakeDouble"})
    public void testIsPresent() {
        final File testFile = BaseTest.createTempFile("hdf5", ".hd5");
        final HDF5File file = new HDF5File(testFile, HDF5File.OpenMode.CREATE);
        Assert.assertFalse(file.isPresent("test-group"));
        Assert.assertFalse(file.isPresent("test-group/lola-run"));
        Assert.assertFalse(file.isPresent("test-group/lola-run/bernie-follows"));
        Assert.assertFalse(file.isPresent("test-group/lola-run/jill-follows"));
        file.makeGroup("test-group/lola-run");
        Assert.assertTrue(file.isPresent("test-group"));
        Assert.assertTrue(file.isPresent("test-group/lola-run"));
        Assert.assertFalse(file.isPresent("test-group/lola-run/bernie-follows"));
        Assert.assertFalse(file.isPresent("test-group/lola-run/jill-follows"));
        file.makeDouble("test-group/lola-run/bernie-follows", 0.1);
        Assert.assertTrue(file.isPresent("test-group"));
        Assert.assertTrue(file.isPresent("test-group/lola-run"));
        Assert.assertTrue(file.isPresent("test-group/lola-run/bernie-follows"));
        Assert.assertFalse(file.isPresent("test-group/lola-run/jill-follows"));
        Assert.assertFalse(file.isPresent("test-group/lola-run/bernie-follows/and-so-does-jill"));
        file.close();

        final HDF5File file2 = new HDF5File(testFile, HDF5File.OpenMode.READ_ONLY);
        Assert.assertTrue(file2.isPresent("test-group"));
        Assert.assertTrue(file2.isPresent("test-group/lola-run"));
        Assert.assertTrue(file2.isPresent("test-group/lola-run/bernie-follows"));
        Assert.assertFalse(file2.isPresent("test-group/lola-run/jill-follows"));
        Assert.assertFalse(file2.isPresent("test-group/lola-run/bernie-follows/and-so-does-jill"));
        file2.close();
    }

    @Test()
    public void testCreateGroup() {
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

    @Test
    public void testCreateLargeMatrix(){
        // Creates a large PoN of junk values and simply tests that these can be written and read.

        // Make a big, fake set of read counts.
        final int numRows = 2500000;
        final int numCols = 10;
        final double mean = 3e-7;
        final double sigma = 1e-9;

        final RealMatrix bigCounts = createMatrixOfGaussianValues(numRows, numCols, mean, sigma);
        final File tempOutputHD5 = IOUtils.createTempFile("big-ol-", ".hd5");
        final HDF5File hdf5File = new HDF5File(tempOutputHD5, HDF5File.OpenMode.CREATE);
        final String hdf5Path = "/test/m";
        hdf5File.makeDoubleMatrix(hdf5Path, bigCounts.getData());
        hdf5File.close();

        final HDF5File hdf5FileForReading = new HDF5File(tempOutputHD5, HDF5File.OpenMode.READ_ONLY);
        final double[][] result = hdf5FileForReading.readDoubleMatrix(hdf5Path);
        final RealMatrix resultAsRealMatrix = new Array2DRowRealMatrix(result);
        Assert.assertTrue(resultAsRealMatrix.getRowDimension() == numRows);
        Assert.assertTrue(resultAsRealMatrix.getColumnDimension() == numCols);
        final RealMatrix readMatrix = new Array2DRowRealMatrix(result);
        PoNTestUtils.assertEqualsMatrix(readMatrix, bigCounts, false);
    }

    private RealMatrix createMatrixOfGaussianValues(int numRows, int numCols, final double mean, final double sigma) {
        final RealMatrix bigCounts = new Array2DRowRealMatrix(numRows, numCols);
        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator();
        randomDataGenerator.reSeed(337337337);
        bigCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return randomDataGenerator.nextGaussian(mean, sigma);
            }
        });
        return bigCounts;
    }
}
