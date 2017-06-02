package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Unit tests for {@link ReadCountRecordUnitTest}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ReadCountRecordUnitTest {

    private static final int[] TEST_VALUE_COUNT = new int[] {0, 1, 10, 50};

    private static final Target TEST_TARGET = new Target("TEST");

    private static final Random rdn = new Random(13);

    private static final String LONG_ARRAY_CONSTRUCTOR = "long[]";
    private static final String DOUBLE_ARRAY_CONSTRUCTOR = "double[]";

    /**
     * Use to give a easy to understand name to tests involving the corresponding record constructor
     * in {@link #CONSTRUCTORS}.
     */
    private static final String[] CONSTRUCTOR_CLASSES =  {DOUBLE_ARRAY_CONSTRUCTOR, LONG_ARRAY_CONSTRUCTOR};

    private static final List<BiFunction<Target, double[], ? extends ReadCountRecord>> CONSTRUCTORS =
            Arrays.asList(ReadCountRecord::new,
                    ReadCountRecordUnitTest::readLongCountRecord);

    private static ReadCountRecord readLongCountRecord(final Target target, final double[] counts) {
        if (counts == null) {
            return new ReadCountRecord(target, (long[]) null);
        }
        final long[] longCounts = new long[counts.length];
        for (int i = 0; i < longCounts.length; i++) {
            longCounts[i] = Math.round(counts[i]);
        }
        return new ReadCountRecord(target, longCounts);
    }

    @Test(dataProvider="testData")
    public void testCreation(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        constructor.apply(TEST_TARGET, counts);
    }

    @Test(dataProvider="testData", dependsOnMethods = "testCreation", expectedExceptions = IllegalArgumentException.class)
    public void testNullTargetCreation(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        constructor.apply(null, counts);
    }

    @Test(dataProvider="testData", dependsOnMethods = "testCreation", expectedExceptions = IllegalArgumentException.class)
    public void testNullCountsCreation(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, @SuppressWarnings("unused") final int size) {
        constructor.apply(TEST_TARGET, null);
    }

    @Test(dataProvider="testData", dependsOnMethods = "testCreation")
    public void testGetDouble(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        final boolean round = testName.equals("long[]");
        final ReadCountRecord record = constructor.apply(TEST_TARGET, counts);
        if (!round) {
            for (int i = 0; i < counts.length; i++) {
                Assert.assertEquals(record.getDouble(i), counts[i], 0.0);
            }
        } else {
            for (int i = 0; i < counts.length; i++) {
                Assert.assertEquals(record.getDouble(i), Math.round(counts[i]), 0.00001);
            }
        }
    }

    @Test(dataProvider="testData", dependsOnMethods = "testGetDouble", expectedExceptions = IllegalArgumentException.class)
    public void testGetDoubleNegativeIndex(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        final ReadCountRecord record = constructor.apply(TEST_TARGET, counts);
        record.getDouble(-1);
    }

    @Test(dataProvider="testData", dependsOnMethods = "testGetDouble", expectedExceptions = IllegalArgumentException.class)
    public void testGetDoubleToolLargeIndex(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        final ReadCountRecord record = constructor.apply(TEST_TARGET, counts);
        record.getDouble(counts.length + 1);
    }

    @Test(dataProvider="testData", dependsOnMethods = "testCreation")
    public void testSize(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        final ReadCountRecord record = constructor.apply(TEST_TARGET, counts);
        Assert.assertEquals(record.size(), size);
    }

    @Test(dataProvider="testData", dependsOnMethods = "testCreation")
    public void testGetTarget(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        final ReadCountRecord record = constructor.apply(TEST_TARGET, counts);
        Assert.assertSame(record.getTarget(), TEST_TARGET);
    }

    @Test(dataProvider="testData", dependsOnMethods = "testCreation")
    public void testGetDoubleCounts(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        final boolean round = testName.equals("long[]");
        final ReadCountRecord record = constructor.apply(TEST_TARGET, counts);
        if (!round) {
            Assert.assertEquals(record.getDoubleCounts(), counts);
        } else {
            final double[] roundedCounts = new double[counts.length];
            for (int i = 0; i < roundedCounts.length; i++) {
                roundedCounts[i] = Math.round(counts[i]);
            }
            Assert.assertEquals(record.getDoubleCounts(), roundedCounts);
        }
    }

    @Test(dataProvider="testNonZeroCountsData", dependsOnMethods = "testCopyCountsTo", expectedExceptions = IllegalArgumentException.class)
    public void testCopyCountsToBeforeBeginning(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        final ReadCountRecord record = constructor.apply(TEST_TARGET, counts);
        final double[] copiedCounts = new double[counts.length];
        record.copyCountsTo(copiedCounts, -1);
    }

    @Test(dataProvider="testNonZeroCountsData", dependsOnMethods = "testCopyCountsTo", expectedExceptions = IllegalArgumentException.class)
    public void testCopyCountsToBeyondEnd(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        final ReadCountRecord record = constructor.apply(TEST_TARGET, counts);
        final double[] copiedCounts = new double[counts.length];
        record.copyCountsTo(copiedCounts, 10);
    }

    @Test(dataProvider="testData", dependsOnMethods = "testCreation")
    public void testCopyCountsTo(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        final boolean round = testName.equals("long[]");
        final ReadCountRecord record = constructor.apply(TEST_TARGET, counts);
        final double[] copiedCounts = new double[counts.length + 20];
        Arrays.fill(copiedCounts, -11);
        record.copyCountsTo(copiedCounts, 10);
        // Check the copied values.
        if (!round) {
            for (int i = 0; i < counts.length; i++) {
                Assert.assertEquals(copiedCounts[10 + i], counts[i], 0.0);
            }
        } else {
            for (int i = 0; i < counts.length; i++) {
                Assert.assertEquals(copiedCounts[10 + i], Math.round(counts[i]), 0.00001);
            }
        }
        // Check that the padding remains intact:
        for (int i = 0; i < 10; i++) {
            Assert.assertEquals(copiedCounts[i], -11, 0.0);
        }
        for (int i = counts.length + 10; i < copiedCounts.length; i++) {
            Assert.assertEquals(copiedCounts[i], -11, 0.0);
        }
    }

    @Test(dataProvider="testNonZeroCountsData", dependsOnMethods = "testAppendCountsTo", expectedExceptions = IllegalStateException.class)
    public void testAppendCountsToBeyondEnd(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);

        final ReadCountRecord record = constructor.apply(TEST_TARGET, counts);
        final List<String> columnNames = Stream.concat(Stream.concat(
                        IntStream.range(0, 10).mapToObj(i -> "pre-padding_" + i),
                        IntStream.range(0, counts.length).mapToObj(i -> "column_" + i)),
                IntStream.range(0, 10).mapToObj(i -> "post-padding_" + i)).collect(Collectors.toList());
        final TableColumnCollection columns = new TableColumnCollection(columnNames);
        final DataLine dataLine = new DataLine(columns, RuntimeException::new);
        final double[] copiedCounts = new double[counts.length + 20];
        Arrays.fill(copiedCounts, -11);
        dataLine.seek(columnNames.size());
        record.appendCountsTo(dataLine);
    }

    @Test(dataProvider="testData", dependsOnMethods = "testCreation")
    public void testAppendCountsTo(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        final double[] counts = generateCounts(size);
        final boolean round = testName.equals("long[]");
        final ReadCountRecord record = constructor.apply(TEST_TARGET, counts);
        final List<String> columnNames = Stream.concat(Stream.concat(
                IntStream.range(0, 10).mapToObj(i -> "pre-padding_" + i),
                IntStream.range(0, counts.length).mapToObj(i -> "column_" + i)),
                IntStream.range(0, 10).mapToObj(i -> "post-padding_" + i)).collect(Collectors.toList());
        final TableColumnCollection columns = new TableColumnCollection(columnNames);
        final DataLine dataLine = new DataLine(columns, RuntimeException::new);
        final double[] copiedCounts = new double[counts.length + 20];
        Arrays.fill(copiedCounts, -11);
        for (int i = 0; i < 10 + 10 + counts.length; i++) {
            dataLine.append("-11");
        }
        dataLine.seek(10);
        record.appendCountsTo(dataLine);
        // Check the copied values.
        if (!round) {
            for (int i = 0; i < counts.length; i++) {
                Assert.assertEquals(dataLine.getDouble(10 + i), counts[i], 0.0);
            }
        } else {
            for (int i = 0; i < counts.length; i++) {
                Assert.assertEquals(dataLine.getDouble(10 + i), Math.round(counts[i]), 0.00001);
            }
        }
        // Check that the padding remains intact:
        for (int i = 0; i < 10; i++) {
            Assert.assertEquals(dataLine.get(i), "-11");
        }
        for (int i = counts.length + 10; i < copiedCounts.length; i++) {
            Assert.assertEquals(dataLine.get(i), "-11");
        }
    }

    @Test(dataProvider="testData", dependsOnMethods = "testCreation", expectedExceptions = IllegalArgumentException.class)
    public void testNullTarget(@SuppressWarnings("unused") final String testName, final BiFunction<Target, double[], ? extends ReadCountRecord> constructor, final int size) {
        constructor.apply(null, generateCounts(size));
    }

    @DataProvider(name="testData")
    public Object[][] testData() {
        final List<Object[]> result = new ArrayList<>();
        for (final int size : TEST_VALUE_COUNT) {
            for (int i = 0; i < CONSTRUCTORS.size(); i++) {
                final BiFunction<Target, double[], ? extends ReadCountRecord> constructor = CONSTRUCTORS.get(i);
                result.add(new Object[] { CONSTRUCTOR_CLASSES[i], constructor, size });
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="testNonZeroCountsData")
    public Object[][] testNonZeroCountsData() {
        final List<Object[]> result = new ArrayList<>();
        for (final int size : TEST_VALUE_COUNT) {
            if (size == 0) {
                continue;
            }
            for (int i = 0; i < CONSTRUCTORS.size(); i++) {
                final BiFunction<Target, double[], ? extends ReadCountRecord> constructor = CONSTRUCTORS.get(i);
                result.add(new Object[] { CONSTRUCTOR_CLASSES[i], constructor, size });
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    private static double[] generateCounts(int countSize) {
        final double[] result = new double[countSize];
        for (int i = 0; i < result.length; i++) {
            result[i] = rdn.nextDouble() * 50;
        }
        return result;
    }
}
