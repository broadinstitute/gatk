package org.broadinstitute.hellbender.tools.pon.coverage.pca;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.Marker;
import org.apache.logging.log4j.message.Message;
import org.apache.logging.log4j.spi.AbstractLogger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.pon.PoNTestUtils;
import org.broadinstitute.hellbender.utils.MatrixSummaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.svd.SVD;
import org.broadinstitute.hellbender.utils.svd.SVDFactory;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.OptionalInt;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Unit test for {@link HDF5PCACoveragePoNCreationUtils}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HDF5PCACoveragePoNCreationUtilsUnitTest extends BaseTest {
    private static final String TEST_DIR = "src/test/resources/org/broadinstitute/hellbender/tools/exome/";
    private static final File TEST_PCOV_FILE = new File(TEST_DIR, "create-pon-control-full.pcov");
    private static final File GT_TARGET_VAR_FILE = new File(TEST_DIR, "dummy_pon_target_variances_matlab.txt");
    private static final File TEST_FULL_PON_FILE = new File(TEST_DIR, "create-pon-all-targets.pon");

    @Test
    public void testRedoReductionWithSameEigensampleCount() {
        final int numEigensamples = 20;
        final File tempOutputPoN = IOUtils.createTempFile("redo-reduction-same-count-", ".pon");
        final File ponFile = PoNTestUtils.createDummyHDF5FilePoN(TEST_PCOV_FILE, numEigensamples);
        HDF5PCACoveragePoNCreationUtils.redoReduction(null, OptionalInt.of(numEigensamples), ponFile, tempOutputPoN, HDF5File.OpenMode.CREATE);
        PoNTestUtils.assertEquivalentPoN(ponFile, tempOutputPoN);
    }

    @Test
    public void testCalculateVariance() {
        /*
        This tests that the post-projection variance is correct.  Ground truth was acquired through a breakpoint (to get input values),
        some text parsing, and matlab.  This simply tests that the variances are correctly calculated from a matrix.
        */
        final int numEigensamples = 20;
        final File ponFile = PoNTestUtils.createDummyHDF5FilePoN(TEST_PCOV_FILE, numEigensamples);
        final double[] gtTargetVariances = ParamUtils.readValuesFromFile(GT_TARGET_VAR_FILE);

        try (final HDF5File ponHDF5File = new HDF5File(ponFile)) {
            final HDF5PCACoveragePoN pon = new HDF5PCACoveragePoN(ponHDF5File);
            final double[] targetVariances = pon.getTargetVariances();
            PoNTestUtils.assertEqualsDoubleArrays(targetVariances, gtTargetVariances, 1E-4);
        }
    }

    @Test
    public void testCalculateVarianceCalculation() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        try (final HDF5File ponHDF5File = new HDF5File(TEST_FULL_PON_FILE)) {
            final PCACoveragePoN pon = new HDF5PCACoveragePoN(ponHDF5File);
            final PCATangentNormalizationResult tnWithSpark = pon.normalizeNormalsInPoN(ctx);
            final double[] variances = MatrixSummaryUtils.getRowVariances(tnWithSpark.getTangentNormalized().counts());
            final double[] targetVariances = pon.getTargetVariances();
            PoNTestUtils.assertEqualsDoubleArrays(targetVariances, variances);
        }
    }

    @Test(dataProvider="readCountAndPercentileData")
    public void testSubsetTargetToUsableOnes(final ReadCountCollection readCount, final double percentile) {
        final Median median = new Median();
        final RealMatrix counts = readCount.counts();
        final double[] targetMedians = IntStream.range(0, counts.getRowDimension())
                .mapToDouble(i -> median.evaluate(counts.getRow(i))).toArray();
        final double threshold = new Percentile(percentile).evaluate(targetMedians);
        final Boolean[] toBeKept = DoubleStream.of(targetMedians)
                .mapToObj(d -> d >= threshold).toArray(Boolean[]::new);
        final int toBeKeptCount = (int) Stream.of(toBeKept).filter(b -> b).count();
        final Pair<ReadCountCollection, double[]> result = HDF5PCACoveragePoNCreationUtils.subsetReadCountsToUsableTargets(readCount, percentile, NULL_LOGGER);
        Assert.assertEquals(result.getLeft().targets().size(), toBeKeptCount);
        Assert.assertEquals(result.getRight().length, toBeKeptCount);
        int nextIndex = 0;
        for (int i = 0; i < toBeKept.length; i++) {
            if (toBeKept[i]) {
                int index = result.getLeft().targets().indexOf(readCount.targets().get(i));
                Assert.assertEquals(index, nextIndex++);
                Assert.assertEquals(counts.getRow(i), result.getLeft().counts().getRow(index));
                Assert.assertEquals(result.getRight()[index], targetMedians[i]);
            } else {
                Assert.assertEquals(result.getLeft().targets().indexOf(readCount.targets().get(i)), -1);
            }
        }
    }

    @Test(dataProvider = "readCountOnlyData")
    public void testNormalizeAndLogReadCounts(final ReadCountCollection readCounts) {
        final RealMatrix counts = readCounts.counts();
        final Median median = new Median();
        final double[] columnMedians = IntStream.range(0, counts.getColumnDimension())
                .mapToDouble(i -> median.evaluate(counts.getColumn(i))).toArray();
        final double epsilon = HDF5PCACoveragePoNCreationUtils.EPSILON;
        final double[][] expected = new double[counts.getRowDimension()][];
        for (int i = 0; i < expected.length; i++) {
            expected[i] = counts.getRow(i).clone();
            for (int j = 0; j < expected[i].length; j++) {
                expected[i][j] /= columnMedians[j];
                if (expected[i][j] < epsilon) {
                    expected[i][j] = epsilon;
                }
                expected[i][j] = Math.log(expected[i][j]) / Math.log(2);
            }
        }
        HDF5PCACoveragePoNCreationUtils.normalizeAndLogReadCounts(readCounts, NULL_LOGGER);
        final RealMatrix newCounts = readCounts.counts();
        Assert.assertEquals(newCounts.getColumnDimension(), expected[0].length);
        Assert.assertEquals(newCounts.getRowDimension(), expected.length);
        for (int i = 0; i < expected.length; i++) {
            for (int j = 0; j < expected[i].length; j++) {
                Assert.assertEquals(newCounts.getEntry(i, j), expected[i][j], 0.000001);
            }
        }
    }

    @Test(dataProvider = "readCountOnlyData")
    public void testSubtractMedianOfMedians(final ReadCountCollection readCounts) {
        final RealMatrix counts = readCounts.counts();
        final Median median = new Median();
        final double[] columnMedians = IntStream.range(0, counts.getColumnDimension())
                .mapToDouble(i -> median.evaluate(counts.getColumn(i))).toArray();
        final double center = median.evaluate(columnMedians);
        final double[][] expected = new double[counts.getRowDimension()][];
        for (int i = 0; i < expected.length; i++) {
            expected[i] = counts.getRow(i).clone();
            for (int j = 0; j < expected[i].length; j++) {
                expected[i][j] -= center;
            }
        }
        HDF5PCACoveragePoNCreationUtils.subtractMedianOfMedians(readCounts, NULL_LOGGER);
        final RealMatrix newCounts = readCounts.counts();
        Assert.assertEquals(newCounts.getColumnDimension(), expected[0].length);
        Assert.assertEquals(newCounts.getRowDimension(), expected.length);
        for (int i = 0; i < expected.length; i++) {
            for (int j = 0; j < expected[i].length; j++) {
                Assert.assertEquals(newCounts.getEntry(i, j), expected[i][j], 0.000001);
            }
        }
    }

    @Test(dataProvider = "readCountOnlyWithDiverseShapeData")
    public void testCalculateReducedPanelAndPInversesUsingJollifesRule(final ReadCountCollection readCounts) {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final ReductionResult result = HDF5PCACoveragePoNCreationUtils.calculateReducedPanelAndPInverses(readCounts, OptionalInt.empty(), NULL_LOGGER, ctx);
        final RealMatrix counts = readCounts.counts();
        Assert.assertNotNull(result);
        Assert.assertNotNull(result.getPseudoInverse());
        Assert.assertNotNull(result.getReducedCounts());
        Assert.assertNotNull(result.getReducedPseudoInverse());
        Assert.assertNotNull(result.getAllSingularValues());
        Assert.assertEquals(counts.getColumnDimension(), result.getAllSingularValues().length);
        Assert.assertEquals(result.getReducedCounts().getRowDimension(), counts.getRowDimension());
        final int eigensamples = result.getReducedCounts().getColumnDimension();
        final Mean mean = new Mean();
        final double meanSingularValue = mean.evaluate(result.getAllSingularValues());
        final double threshold = HDF5PCACoveragePoNCreationUtils.JOLLIFES_RULE_MEAN_FACTOR * meanSingularValue;
        final int expectedEigensamples = (int) DoubleStream.of(result.getAllSingularValues()).filter(d -> d >= threshold).count();
        Assert.assertTrue(eigensamples <= counts.getColumnDimension());
        Assert.assertEquals(eigensamples, expectedEigensamples);
        assertPseudoInverse(counts, result.getPseudoInverse());
        assertPseudoInverse(result.getReducedCounts(), result.getReducedPseudoInverse());
    }

    @Test(dataProvider = "readCountOnlyWithDiverseShapeData")
    public void testCalculateReducedPanelAndPInversesKeepingHalfOfAllColumns(final ReadCountCollection readCounts) {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final ReductionResult result = HDF5PCACoveragePoNCreationUtils.calculateReducedPanelAndPInverses(readCounts, OptionalInt.of(readCounts.columnNames().size() / 2), NULL_LOGGER, ctx);
        final RealMatrix counts = readCounts.counts();
        Assert.assertNotNull(result);
        Assert.assertNotNull(result.getPseudoInverse());
        Assert.assertNotNull(result.getReducedCounts());
        Assert.assertNotNull(result.getReducedPseudoInverse());
        Assert.assertNotNull(result.getAllSingularValues());
        Assert.assertEquals(counts.getColumnDimension(), result.getAllSingularValues().length);
        Assert.assertEquals(result.getReducedCounts().getRowDimension(), counts.getRowDimension());
        Assert.assertEquals(result.getReducedCounts().getColumnDimension(), readCounts.columnNames().size() / 2);
        final int eigensamples = result.getReducedCounts().getColumnDimension();
        Assert.assertEquals(eigensamples, readCounts.columnNames().size() / 2);
        assertPseudoInverse(counts, result.getPseudoInverse());
        assertPseudoInverse(result.getReducedCounts(), result.getReducedPseudoInverse());
    }

    @Test(dataProvider = "readCountOnlyWithDiverseShapeData")
    public void testCalculateReducedPanelAndPInversesKeepingAllColumns(final ReadCountCollection readCounts) {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final ReductionResult result = HDF5PCACoveragePoNCreationUtils.calculateReducedPanelAndPInverses(readCounts, OptionalInt.of(readCounts.columnNames().size()), NULL_LOGGER, ctx);
        final RealMatrix counts = readCounts.counts();
        Assert.assertNotNull(result);
        Assert.assertNotNull(result.getPseudoInverse());
        Assert.assertNotNull(result.getReducedCounts());
        Assert.assertNotNull(result.getReducedPseudoInverse());
        Assert.assertNotNull(result.getAllSingularValues());
        Assert.assertEquals(counts.getColumnDimension(), result.getAllSingularValues().length);
        Assert.assertEquals(result.getReducedCounts().getRowDimension(), counts.getRowDimension());
        Assert.assertEquals(result.getReducedCounts().getColumnDimension(), readCounts.columnNames().size());
        final int eigensamples = result.getReducedCounts().getColumnDimension();
        Assert.assertEquals(eigensamples, readCounts.columnNames().size());
        assertPseudoInverse(counts, result.getPseudoInverse());
        assertPseudoInverse(result.getReducedCounts(), result.getReducedPseudoInverse());
    }

    private static void assertPseudoInverse(final RealMatrix A, final RealMatrix pinvA) {
        Assert.assertEquals(A.getRowDimension(), pinvA.getColumnDimension());
        Assert.assertEquals(A.getColumnDimension(), pinvA.getRowDimension());

        // 1. condition: A * pinvA * A = A.
        final RealMatrix mustBeZeros1 = A.multiply(pinvA).multiply(A).subtract(A);
        for (int i = 0; i < mustBeZeros1.getRowDimension(); i++) {
            for (int j = 0; j < mustBeZeros1.getColumnDimension(); j++) {
                Assert.assertEquals(mustBeZeros1.getEntry(i, j), 0.0, 0.00001);
            }
        }
        // 2. condition: pinvA * A * pinvA = pinvA.
        final RealMatrix mustBeZeros2 = pinvA.multiply(A).multiply(pinvA).subtract(pinvA);
        for (int i = 0; i < mustBeZeros2.getRowDimension(); i++) {
            for (int j = 0; j < mustBeZeros2.getColumnDimension(); j++) {
                Assert.assertEquals(mustBeZeros2.getEntry(i, j), 0.0, 0.00001);
            }
        }

        // 3. condition: t( A * pinvA ) = A * pinvA
        final RealMatrix mustBeZeros3 = A.multiply(pinvA).transpose().subtract(A.multiply(pinvA));
        for (int i = 0; i < mustBeZeros3.getRowDimension(); i++) {
            for (int j = 0; j < mustBeZeros3.getColumnDimension(); j++) {
                Assert.assertEquals(mustBeZeros3.getEntry(i, j), 0.0, 0.00001);
            }
        }

        // 4. condition t( pinvA * A ) = pinvA * A
        final RealMatrix mustBeZeros4 = pinvA.multiply(A).transpose().subtract(pinvA.multiply(A));
        for (int i = 0; i < mustBeZeros4.getRowDimension(); i++) {
            for (int j = 0; j < mustBeZeros4.getColumnDimension(); j++) {
                Assert.assertEquals(mustBeZeros4.getEntry(i, j), 0.0, 0.00001);
            }
        }

    }

    @DataProvider(name="readCountOnlyWithDiverseShapeData")
    public Object[][] readCountOnlyWithDiverseShapeData() {
        final List<Object[]> result = new ArrayList<>(4);
        final Random rdn = new Random(31);
        final int[] columnCounts = new int[] { 10, 100, 100, 200};
        final int[] targetCounts = new int[] { 100, 100, 200, 200 };
        for (int k = 0; k < columnCounts.length; k++) {
            final List<String> columnNames = IntStream.range(0, columnCounts[k]).mapToObj(i -> "sample_" + (i + 1)).collect(Collectors.toList());
            final List<Target> targets = IntStream.range(0, targetCounts[k]).mapToObj(i -> new Target("target_" + (i+1))).collect(Collectors.toList());
            final double[][] counts = new double[targetCounts[k]][columnCounts[k]];
            for (int i = 0; i < counts.length; i++) {
                for (int j = 0; j < counts[0].length; j++) {
                    counts[i][j] = rdn.nextDouble();
                }
            }
            final ReadCountCollection readCounts = new ReadCountCollection(targets, columnNames,
                    new Array2DRowRealMatrix(counts, false));
            result.add(new Object[]{readCounts });
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="readCountOnlyData")
    public Object[][] readCountOnlyData() {
        final int repeats = 4;
        final List<Object[]> result = new ArrayList<>(repeats);
        final Random rdn = new Random(13);
        final int columnCount = 100;
        final int targetCount = 100;
        final List<String> columnNames = IntStream.range(0, columnCount).mapToObj(i -> "sample_" + (i + 1)).collect(Collectors.toList());
        final List<Target> targets = IntStream.range(0, targetCount).mapToObj(i -> new Target("target_" + (i+1))).collect(Collectors.toList());
        for (int k = 0; k < repeats; k++) {
            final double[][] counts = new double[columnCount][targetCount];
            for (int i = 0; i < counts.length; i++) {
                for (int j = 0; j < counts[0].length; j++) {
                    counts[i][j] = rdn.nextDouble();
                }
            }
            final ReadCountCollection readCounts = new ReadCountCollection(targets, columnNames,
                    new Array2DRowRealMatrix(counts, false));
            result.add(new Object[]{readCounts });
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="singleEigensample")
    public Object[][] simpleEigensampleData() {
        final List<Object[]> result = new ArrayList<>();
        final int NUM_TARGETS = 10;
        final int NUM_SAMPLES = 5;
        final List<Target> targets = IntStream.range(0, NUM_TARGETS).boxed()
                .map(i -> new Target("target_" + i, new SimpleInterval("1", 100*i + 1, 100*i + 5)))
                .collect(Collectors.toList());
        final List<String> columnNames = IntStream.range(0, NUM_SAMPLES).boxed()
                .map(i -> "sample_" + i)
                .collect(Collectors.toList());

        double [][] countsArray = new double[NUM_TARGETS][NUM_SAMPLES];
        final RealMatrix counts = new Array2DRowRealMatrix(countsArray);

        // All row data is the same (0,1,2,3,4...)
        final double [] rowData = IntStream.range(0, NUM_SAMPLES).boxed()
                .mapToDouble(i -> i).toArray();
        for (int i = 0; i < NUM_TARGETS; i++) {
            counts.setRow(i, rowData);
        }

        new ReadCountCollection(targets, columnNames, counts);

        result.add(new Object[] {new ReadCountCollection(targets, columnNames, counts)});
        return result.toArray(new Object[result.size()][]);
    }

    @Test(dataProvider = "singleEigensample" )
    public void testDetermineNumberOfEigensamplesNoSpark(final ReadCountCollection logNormals){

        final SVD logNormalsSVD = SVDFactory.createSVD(logNormals.counts());

        final int actualNumber = HDF5PCACoveragePoNCreationUtils.determineNumberOfEigensamples(OptionalInt.empty(), logNormals.columnNames().size(), logNormalsSVD, NULL_LOGGER);
        Assert.assertEquals(actualNumber, 1);
    }

    @Test(dataProvider = "singleEigensample" )
    public void testDetermineNumberOfEigensamplesSpark(final ReadCountCollection logNormals){
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final SVD logNormalsSVD = SVDFactory.createSVD(logNormals.counts(), ctx);

        final int actualNumber = HDF5PCACoveragePoNCreationUtils.determineNumberOfEigensamples(OptionalInt.empty(), logNormals.columnNames().size(), logNormalsSVD, NULL_LOGGER);
        Assert.assertEquals(actualNumber, 1);
    }

    @SuppressWarnings("serial")
    private static final Logger NULL_LOGGER = new AbstractLogger() {
        @Override
        public Level getLevel() {
            return Level.ALL;
        }

        @Override
        public boolean isEnabled(Level level, Marker marker, Message message, Throwable t) {
            return false;
        }

        @Override
        public boolean isEnabled(Level level, Marker marker, Object message, Throwable t) {
            return false;
        }

        @Override
        public boolean isEnabled(Level level, Marker marker, String message, Throwable t) {
            return false;
        }

        @Override
        public boolean isEnabled(Level level, Marker marker, String message) {
            return false;
        }

        @Override
        public boolean isEnabled(Level level, Marker marker, String message, Object... params) {
            return false;
        }

        @Override
        public void logMessage(String fqcn, Level level, Marker marker, Message message, Throwable t) {
            // Do nothing.
        }
    };
}