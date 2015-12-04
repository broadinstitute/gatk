package org.broadinstitute.hellbender.utils.hdf5;

import org.apache.commons.collections4.list.SetUniqueList;
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
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
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
 * Unit test for most relevant static routines in {@link HDF5PoNCreator}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class HDF5PoNCreatorUnitTest extends BaseTest {
    private final static String TEST_DIR = "src/test/resources/org/broadinstitute/hellbender/tools/exome/";
    private final static String TEST_PCOV_FILE = "create-pon-control-full.pcov";
    private final static String GT_TARGET_VAR_FILE = TEST_DIR + "dummy_pon_target_variances_matlab.txt";
    private final static String TEST_FULL_PON = TEST_DIR + "create-pon-all-targets.pon";

    @Test
    public void testCreateVariance() {
        /*
        This tests that the post-projection variance is correct.  Ground truth was acquired through a breakpoint (to get input values),
        some text parsing, and matlab.  This simply tests that the variances are correctly calculated from a matrix.
        */
        final int numEigenSamples = 20;
        final File ponFile = PoNTestUtils.createDummyHDF5FilePoN(new File(TEST_DIR + TEST_PCOV_FILE), numEigenSamples);
        final double[] gtTargetVariances = ParamUtils.readValuesFromFile(new File(GT_TARGET_VAR_FILE));

        try (final HDF5File ponHDF5File = new HDF5File(ponFile)) {
            final HDF5PoN pon = new HDF5PoN(ponHDF5File);
            final double[] targetVariances = pon.getTargetVariances();
            PoNTestUtils.assertEqualsDoubleArrays(targetVariances, gtTargetVariances);
        }
    }

    @Test
    public void testCreateVarianceCalculation() {

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final File ponFile = new File(TEST_FULL_PON);
        try (final HDF5File ponHDF5File = new HDF5File(ponFile)) {
            final PoN pon = new HDF5PoN(ponHDF5File);
            final TangentNormalizationResult tnWithSpark = TangentNormalizer.tangentNormalizeNormalsInPoN(pon, ctx);
            final double[] variances = HDF5PoNCreator.calculateRowVariances(tnWithSpark.getTangentNormalized().counts());
            final double[] targetVariances = pon.getTargetVariances();
            PoNTestUtils.assertEqualsDoubleArrays(targetVariances, variances);
        }
    }

    @Test(dataProvider="normalizeReadCountByTargetFactorsData")
    public void testNormalizeReadCountByTargetFactors(final ReadCountCollection readCount, final double[] targetFactors) {
        final RealMatrix before = readCount.counts().copy();
        HDF5PoNCreator.normalizeReadCountsByTargetFactors(readCount, targetFactors);
        final RealMatrix after = readCount.counts();
        Assert.assertEquals(before.getColumnDimension(), after.getColumnDimension());
        Assert.assertEquals(before.getRowDimension(), after.getRowDimension());
        for (int i = 0; i < after.getColumnDimension(); i++) {
            for (int j = 0; j < after.getRowDimension(); j++) {
                final double beforeValue = before.getEntry(j, i);
                final double afterValue = after.getEntry(j, i);
                Assert.assertEquals(afterValue, beforeValue / targetFactors[j], 0.001);
            }
        }
    }

    @Test(dataProvider = "tooManyZerosData")
    public void testRemoveColumnsWithTooManyZeros(final ReadCountCollection readCount) {
        final RealMatrix counts = readCount.counts();
        final int[] numberOfZeros = IntStream.range(0, counts.getColumnDimension())
                .map(i -> (int) DoubleStream.of(counts.getColumn(i)).filter(d -> d == 0.0).count()).toArray();
        final int maximumNumberOfZeros = IntStream.of(numberOfZeros).max().getAsInt();
        for (int maxZeros = 0; maxZeros < maximumNumberOfZeros; maxZeros++) {
            final int maxZerosThres = maxZeros;
            final int expectedRemainingCount = (int) IntStream.of(numberOfZeros).filter(i -> i <= maxZerosThres).count();
            if (expectedRemainingCount == 0) {
                try {
                    HDF5PoNCreator.removeColumnsWithTooManyZeros(readCount, maxZeros, NULL_LOGGER);
                } catch (final UserException.BadInput ex) {
                    // expected.
                    continue;
                }
                Assert.fail("expects an exception");
            }
            final ReadCountCollection rc = HDF5PoNCreator.removeColumnsWithTooManyZeros(readCount, maxZeros, NULL_LOGGER);
            Assert.assertEquals(rc.columnNames().size(), expectedRemainingCount);
            final int[] newIndices = new int[expectedRemainingCount];
            int nextIndex = 0;
            for (int i = 0; i < readCount.columnNames().size(); i++) {
                final String name = readCount.columnNames().get(i);
                final int newIndex = rc.columnNames().indexOf(name);
                if (numberOfZeros[i] <= maxZeros) {
                    Assert.assertTrue(newIndex >= 0);
                    newIndices[nextIndex++] = i;
                } else {
                    Assert.assertEquals(newIndex, -1);
                }
            }
            Assert.assertEquals(nextIndex, expectedRemainingCount);
            for (int i = 1; i < newIndices.length; i++) {
                Assert.assertTrue(newIndices[i - 1] < newIndices[i]);
            }
        }
    }

    @Test(dataProvider = "tooManyZerosData")
    public void testRemoveTargetsWithTooManyZeros(final ReadCountCollection readCount) {
        final RealMatrix counts = readCount.counts();
        final int[] numberOfZeros = IntStream.range(0, counts.getRowDimension())
                .map(i -> (int) DoubleStream.of(counts.getRow(i)).filter(d -> d == 0.0).count()).toArray();
        final int maximumNumberOfZeros = IntStream.of(numberOfZeros).max().getAsInt();
        for (int maxZeros = 0; maxZeros < maximumNumberOfZeros; maxZeros++) {
            final int maxZerosThres = maxZeros;
            final int expectedRemainingCount = (int) IntStream.of(numberOfZeros).filter(i -> i <= maxZerosThres).count();
            if (expectedRemainingCount == 0) {
                try {
                    HDF5PoNCreator.removeTargetsWithTooManyZeros(readCount, maxZeros, NULL_LOGGER);
                } catch (final UserException.BadInput ex) {
                    // expected.
                    continue;
                }
                Assert.fail("expects an exception");
            }
            final ReadCountCollection rc = HDF5PoNCreator.removeTargetsWithTooManyZeros(readCount, maxZeros, NULL_LOGGER);
            Assert.assertEquals(rc.targets().size(), expectedRemainingCount);
            int nextIndex = 0;
            for (int i = 0; i < readCount.targets().size(); i++) {
                final Target target = readCount.targets().get(i);
                final int newIndex = rc.targets().indexOf(target);
                if (numberOfZeros[i] <= maxZeros) {
                    Assert.assertTrue(newIndex >= 0, " " + numberOfZeros[i] + " " + maxZeros);
                    Assert.assertEquals(newIndex, nextIndex++);
                } else {
                    Assert.assertEquals(newIndex, -1);
                }
            }
            Assert.assertEquals(nextIndex, expectedRemainingCount);
        }
    }

    @Test(dataProvider="readCountAndPercentileData")
    public void testTruncateExtremeCounts(final ReadCountCollection readCount, final double percentile) {
        final RealMatrix counts = readCount.counts();
        final double[] allCounts = Stream.of(counts.getData())
                .flatMap(row -> DoubleStream.of(row).boxed())
                .mapToDouble(Double::doubleValue).toArray();
        final double bottom = new Percentile(percentile).evaluate(allCounts);
        final double top = new Percentile(100 - percentile).evaluate(allCounts);
        final double[][] expected = new double[counts.getRowDimension()][];
        for (int i = 0; i < expected.length; i++) {
            expected[i] = DoubleStream.of(counts.getRow(i)).map(d -> d < bottom ? bottom : (d > top) ? top : d).toArray();
        }
        HDF5PoNCreator.truncateExtremeCounts(readCount, percentile, NULL_LOGGER);
        final RealMatrix newCounts = readCount.counts();
        Assert.assertEquals(newCounts.getRowDimension(), newCounts.getRowDimension());
        Assert.assertEquals(newCounts.getColumnDimension(), newCounts.getColumnDimension());
        for (int i = 0; i < expected.length; i++) {
            for (int j = 0; j < expected[i].length; j++) {
                Assert.assertEquals(newCounts.getEntry(i, j), expected[i][j]);
            }
        }
    }

    @Test(dataProvider="readCountAndPercentileData")
    public void testExtremeMedianColumnsData(final ReadCountCollection readCount, final double percentile) {
        final Median median = new Median();
        final RealMatrix counts = readCount.counts();
        final double[] columnMedians = IntStream.range(0, counts.getColumnDimension())
                .mapToDouble(i -> median.evaluate(counts.getColumn(i))).toArray();
        final double top = new Percentile(100 - percentile).evaluate(columnMedians);
        final double bottom = new Percentile(percentile).evaluate(columnMedians);
        final Boolean[] toBeKept = DoubleStream.of(columnMedians)
                .mapToObj(d -> d <= top && d >= bottom).toArray(Boolean[]::new);
        final int toBeKeptCount = (int) Stream.of(toBeKept).filter(b -> b).count();
        final ReadCountCollection result = HDF5PoNCreator.removeColumnsWithExtremeMedianCounts(readCount, percentile, NULL_LOGGER);
        Assert.assertEquals(result.columnNames().size(), toBeKeptCount);
        int nextIndex = 0;
        for (int i = 0; i < toBeKept.length; i++) {
            if (toBeKept[i]) {
                int index = result.columnNames().indexOf(readCount.columnNames().get(i));
                Assert.assertEquals(index, nextIndex++);
                Assert.assertEquals(counts.getColumn(i), result.counts().getColumn(index));
            } else {
                Assert.assertEquals(result.columnNames().indexOf(readCount.columnNames().get(i)), -1);
            }
        }
    }

    @Test(dataProvider = "tooManyZerosData")
    public void testImputeZeroCounts(final ReadCountCollection readCounts) {
        final Median median = new Median();
        final RealMatrix counts = readCounts.counts();
        final double[] targetNonZeroMedians = IntStream.range(0, counts.getRowDimension())
                .mapToDouble(i -> median.evaluate(DoubleStream.of(counts.getRow(i)).filter(d -> d != 0.0).toArray())).toArray();
        final double[][] expected = new double[counts.getRowDimension()][];
        final double[][] original = counts.getData();
        for (int i = 0; i < expected.length; i++) {
            final double[] rowCounts = counts.getRow(i).clone();
            expected[i] = rowCounts;
            for (int j = 0; j < expected[i].length; j++) {
                if (expected[i][j] == 0.0) {
                    expected[i][j] = targetNonZeroMedians[i];
                }
            }
        }
        HDF5PoNCreator.imputeZerosCounts(readCounts, NULL_LOGGER);
        final RealMatrix newCounts = readCounts.counts();
        Assert.assertEquals(newCounts.getColumnDimension(), expected[0].length);
        Assert.assertEquals(newCounts.getRowDimension(), expected.length);
        for (int i = 0; i < expected.length; i++) {
            for (int j = 0; j < expected[i].length; j++) {
                Assert.assertEquals(newCounts.getEntry(i, j), expected[i][j], "i,j == " + i + "," + j + " " + original[i][j]);
            }
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
        final Pair<ReadCountCollection, double[]> result = HDF5PoNCreator.subsetReadCountsToUsableTargets(readCount, percentile, NULL_LOGGER);
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
        final double epsilon = HDF5PoNCreator.EPSILON;
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
        HDF5PoNCreator.normalizeAndLogReadCounts(readCounts, NULL_LOGGER);
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
    public void testSubtractBGSCenters(final ReadCountCollection readCounts) {
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
        HDF5PoNCreator.subtractBGSCenter(readCounts, NULL_LOGGER);
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
        final ReductionResult result = HDF5PoNCreator.calculateReducedPanelAndPInverses(readCounts, OptionalInt.empty(), NULL_LOGGER, ctx);
        final RealMatrix counts = readCounts.counts();
        Assert.assertNotNull(result);
        Assert.assertNotNull(result.getPseudoInverse());
        Assert.assertNotNull(result.getReducedCounts());
        Assert.assertNotNull(result.getReducedInverse());
        Assert.assertNotNull(result.getAllSingularValues());
        Assert.assertEquals(counts.getColumnDimension(), result.getAllSingularValues().length);
        Assert.assertEquals(result.getReducedCounts().getRowDimension(), counts.getRowDimension());
        final int eigenSamples = result.getReducedCounts().getColumnDimension();
        final Mean mean = new Mean();
        final double meanSingularValue = mean.evaluate(result.getAllSingularValues());
        final double threshold = HDF5PoNCreator.JOLLIFES_RULE_MEAN_FACTOR * meanSingularValue;
        final int expectedEigenSamples = (int) DoubleStream.of(result.getAllSingularValues()).filter(d -> d >= threshold).count();
        Assert.assertTrue(eigenSamples <= counts.getColumnDimension());
        Assert.assertEquals(eigenSamples, expectedEigenSamples);
        assertPseudoInverse(counts, result.getPseudoInverse());
        assertPseudoInverse(result.getReducedCounts(), result.getReducedInverse());
    }

    @Test(dataProvider = "readCountOnlyWithDiverseShapeData")
    public void testCalculateReducedPanelAndPInversesKeepingHalfOfAllColumns(final ReadCountCollection readCounts) {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final ReductionResult result = HDF5PoNCreator.calculateReducedPanelAndPInverses(readCounts, OptionalInt.of(readCounts.columnNames().size() / 2), NULL_LOGGER, ctx);
        final RealMatrix counts = readCounts.counts();
        Assert.assertNotNull(result);
        Assert.assertNotNull(result.getPseudoInverse());
        Assert.assertNotNull(result.getReducedCounts());
        Assert.assertNotNull(result.getReducedInverse());
        Assert.assertNotNull(result.getAllSingularValues());
        Assert.assertEquals(counts.getColumnDimension(), result.getAllSingularValues().length);
        Assert.assertEquals(result.getReducedCounts().getRowDimension(), counts.getRowDimension());
        Assert.assertEquals(result.getReducedCounts().getColumnDimension(), readCounts.columnNames().size() / 2);
        final int eigenSamples = result.getReducedCounts().getColumnDimension();
        Assert.assertEquals(eigenSamples, readCounts.columnNames().size() / 2);
        assertPseudoInverse(counts, result.getPseudoInverse());
        assertPseudoInverse(result.getReducedCounts(), result.getReducedInverse());
    }

    @Test
    public void testRedoReductionWithSameEigenSampleCount() {
        final int numEigenSamples = 20;
        final File tempOutputPoN = IOUtils.createTempFile("redo-reduction-same-count-", ".pon");
        final File ponFile = PoNTestUtils.createDummyHDF5FilePoN(new File(TEST_DIR + TEST_PCOV_FILE), numEigenSamples);
        HDF5PoNCreator.redoReduction(null, OptionalInt.of(numEigenSamples), ponFile, tempOutputPoN);
        PoNTestUtils.assertEquivalentPoN(ponFile, tempOutputPoN);
    }

    @Test(dataProvider = "readCountOnlyWithDiverseShapeData")
    public void testCalculateReducedPanelAndPInversesKeepingAllColumns(final ReadCountCollection readCounts) {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final ReductionResult result = HDF5PoNCreator.calculateReducedPanelAndPInverses(readCounts, OptionalInt.of(readCounts.columnNames().size()), NULL_LOGGER, ctx);
        final RealMatrix counts = readCounts.counts();
        Assert.assertNotNull(result);
        Assert.assertNotNull(result.getPseudoInverse());
        Assert.assertNotNull(result.getReducedCounts());
        Assert.assertNotNull(result.getReducedInverse());
        Assert.assertNotNull(result.getAllSingularValues());
        Assert.assertEquals(counts.getColumnDimension(), result.getAllSingularValues().length);
        Assert.assertEquals(result.getReducedCounts().getRowDimension(), counts.getRowDimension());
        Assert.assertEquals(result.getReducedCounts().getColumnDimension(), readCounts.columnNames().size());
        final int eigenSamples = result.getReducedCounts().getColumnDimension();
        Assert.assertEquals(eigenSamples, readCounts.columnNames().size());
        assertPseudoInverse(counts, result.getPseudoInverse());
        assertPseudoInverse(result.getReducedCounts(), result.getReducedInverse());
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
            final ReadCountCollection readCounts = new ReadCountCollection(
                    SetUniqueList.setUniqueList(new ArrayList<>(targets)),
                    SetUniqueList.setUniqueList(new ArrayList<>(columnNames)),
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
            final ReadCountCollection readCounts = new ReadCountCollection(
                    SetUniqueList.setUniqueList(new ArrayList<>(targets)),
                    SetUniqueList.setUniqueList(new ArrayList<>(columnNames)),
                    new Array2DRowRealMatrix(counts, false));
            result.add(new Object[]{readCounts });
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="readCountAndPercentileData")
    public Object[][] readCountAndPercentileData() {
        final double[] percentiles = new double[] { 1.0, 2.5, 5.0 , 10.0, 25.0 };
        final List<Object[]> result = new ArrayList<>();
        final Random rdn = new Random(13);
        final int columnCount = 100;
        final int targetCount = 100;
        final List<String> columnNames = IntStream.range(0, columnCount).mapToObj(i -> "sample_" + (i + 1)).collect(Collectors.toList());
        final List<Target> targets = IntStream.range(0, targetCount).mapToObj(i -> new Target("target_" + (i+1))).collect(Collectors.toList());
        for (final double percentile : percentiles) {
            final double[][] counts = new double[columnCount][targetCount];
            for (int i = 0; i < counts.length; i++) {
                for (int j = 0; j < counts[0].length; j++) {
                    counts[i][j] = rdn.nextDouble();
                }
            }
            final ReadCountCollection readCounts = new ReadCountCollection(
                    SetUniqueList.setUniqueList(new ArrayList<>(targets)),
                    SetUniqueList.setUniqueList(new ArrayList<>(columnNames)),
                    new Array2DRowRealMatrix(counts, false));
            result.add(new Object[]{readCounts, percentile});
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="tooManyZerosData")
    public Object[][] tooManyZerosData() {
        final double[] zeroProbabilities = new double[] { .001, .01, .02, 0.1 };
        final List<Object[]> result = new ArrayList<>();
        final Random rdn = new Random(13);
        final int columnCount = 100;
        final int targetCount = 100;
        final List<String> columnNames = IntStream.range(0, columnCount).mapToObj(i -> "sample_" + (i + 1)).collect(Collectors.toList());
        final List<Target> targets = IntStream.range(0, targetCount).mapToObj(i -> new Target("target_" + (i+1))).collect(Collectors.toList());
        for (final double zeroProbability : zeroProbabilities) {
            final double[][] counts = new double[columnCount][targetCount];
            for (int i = 0; i < counts.length; i++) {
                for (int j = 0; j < counts[0].length; j++) {
                    counts[i][j] = rdn.nextDouble() <= zeroProbability ? 0.0 : rdn.nextDouble();
                }
            }
            final ReadCountCollection readCounts = new ReadCountCollection(
                    SetUniqueList.setUniqueList(new ArrayList<>(targets)),
                    SetUniqueList.setUniqueList(new ArrayList<>(columnNames)),
                    new Array2DRowRealMatrix(counts, false));
            result.add(new Object[] { readCounts });
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="normalizeReadCountByTargetFactorsData")
    public Object[][] normalizeReadCountByTargetFactorsData() {
        final List<Object[]> result = new ArrayList<>(1);
        @SuppressWarnings("serial")
        final List<Target> targets = new ArrayList<Target>() {{
            add(new Target("A"));
            add(new Target("B"));
            add(new Target("C"));
        }};

        @SuppressWarnings("serial")
        final List<String> columnNames = new ArrayList<String>() {{
                add("1");
                add("2");
                add("3");
        }};
        result.add(new Object[] {
                new ReadCountCollection(
                   SetUniqueList.setUniqueList(targets),
                   SetUniqueList.setUniqueList(columnNames),
                   new Array2DRowRealMatrix(new double[][] {
                           new double[] { 1.1, 2.2, 3.3 },
                           new double[] { 0.1, 0.2, 0.3},
                           new double[] { 11.1, 22.2, 33.3 } }, false))
        , new double[] { 100.0, 200.0, 300.0 }});
        return result.toArray(new Object[1][]);
    }

    @DataProvider(name="singleEigenExample")
    public Object[][] simpleEigenExampleData() {
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

        new ReadCountCollection(
                SetUniqueList.setUniqueList(targets),
                SetUniqueList.setUniqueList(columnNames),
                counts);

        result.add(new Object[] {new ReadCountCollection(
                SetUniqueList.setUniqueList(targets),
                SetUniqueList.setUniqueList(columnNames),
                counts)});
        return result.toArray(new Object[result.size()][]);
    }

    @Test(dataProvider = "singleEigenExample" )
    public void testDetermineNumberOfEigenSamplesNoSpark(final ReadCountCollection logNormals){

        final SVD logNormalsSVD = SVDFactory.createSVD(logNormals.counts());

        final int actualNumber = HDF5PoNCreator.determineNumberOfEigenSamples(OptionalInt.empty(), logNormals.columnNames().size(), logNormalsSVD, NULL_LOGGER);
        Assert.assertEquals(actualNumber, 1);
    }

    @Test(dataProvider = "singleEigenExample" )
    public void testDetermineNumberOfEigenSamplesSpark(final ReadCountCollection logNormals){
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final SVD logNormalsSVD = SVDFactory.createSVD(logNormals.counts(), ctx);

        final int actualNumber = HDF5PoNCreator.determineNumberOfEigenSamples(OptionalInt.empty(), logNormals.columnNames().size(), logNormalsSVD, NULL_LOGGER);
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
