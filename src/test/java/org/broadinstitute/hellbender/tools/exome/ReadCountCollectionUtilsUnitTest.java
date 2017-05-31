package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.Marker;
import org.apache.logging.log4j.message.Message;
import org.apache.logging.log4j.spi.AbstractLogger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Unit tests for {@link ReadCountCollectionUtils}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ReadCountCollectionUtilsUnitTest {

    private static final File TEST_FILE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private static final File FULL_CORRECT_FILE = new File(TEST_FILE_DIR, "rcc-test-full-counts.txt");

    private static final String CONTIG_START_END = TargetTableColumn.CONTIG.toString()
            + "\t" + TargetTableColumn.START.toString() + "\t" + TargetTableColumn.END.toString();

    private static final String CONTIG_START_END_NAME = CONTIG_START_END + "\t" + TargetTableColumn.NAME.toString();

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullFile() throws IOException {
        ReadCountCollectionUtils.parse(null);
    }

    @Test
    public void testReadFullCorrectFile() throws IOException {
        ReadCountCollectionUtils.parse(FULL_CORRECT_FILE);
    }

    @Test(expectedExceptions = IOException.class)
    public void testReadNonExistingFile() throws IOException {
        ReadCountCollectionUtils.parse(new File(TEST_FILE_DIR, "false-" + FULL_CORRECT_FILE.getName()));
    }

    @Test(expectedExceptions = IOException.class)
    public void testReadFromDirectory() throws IOException {
        ReadCountCollectionUtils.parse(TEST_FILE_DIR);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testEmptyFile() throws IOException {
        final File emptyFile = createTempFile();
        ReadCountCollectionUtils.parse(emptyFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testCommentHeaderOnly() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoCountColumns() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println(CONTIG_START_END_NAME);
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoSpecialColumns() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println("SAMPLE1\tSAMPLE2");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testRepeatedSampleNames() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println(CONTIG_START_END_NAME + "\tSAMPLE1\tSAMPLE1");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testRepeatedSpecialColumnNames() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println("CONTIG\tCONTIG\tEND\tNAME\tSAMPLE1\tSAMPLE1");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadEmptyCounts() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println(CONTIG_START_END_NAME + "\tSAMPLE1\tSAMPLE2");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testRepeatedTargets() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println(CONTIG_START_END_NAME + "\tSAMPLE1\tSAMPLE2");
        writer.println("1\t100\t200\ttgt_0\t1.1\t2.2");
        writer.println("1\t201\t300\ttgt_1\t1.1\t2.2");
        writer.println("2\t400\t500\ttgt_0\t-1.1E-7\t-2.2E-8");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test
    public void testReadFullFormattedFile() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println(CONTIG_START_END_NAME + "\tSAMPLE1\tSAMPLE2");
        writer.println("1\t100\t200\ttgt_0\t1.1\t2.2");
        writer.println("2\t200\t300\ttgt_1\t-1.1E-7\t-2.2E-8");
        writer.close();
        final ReadCountCollection subject = ReadCountCollectionUtils.parse(testFile);
        Assert.assertNotNull(subject);
        Assert.assertEquals(subject.columnNames(), Arrays.asList("SAMPLE1", "SAMPLE2"));
        Assert.assertEquals(subject.targets().stream().map(Target::getName).collect(Collectors.toList()), Arrays.asList("tgt_0", "tgt_1"));
        Assert.assertEquals(subject.targets().stream().map(Target::getInterval).collect(Collectors.toList()), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("2", 200, 300)));
        Assert.assertEquals(subject.targets().size(), 2);
        final RealMatrix counts = subject.counts();
        Assert.assertEquals(counts.getRowDimension(), 2);
        Assert.assertEquals(counts.getColumnDimension(), 2);
        Assert.assertEquals(counts.getEntry(0, 0), 1.1, 0.0001);
        Assert.assertEquals(counts.getEntry(0, 1), 2.2, 0.0001);
        Assert.assertEquals(counts.getEntry(1, 0), -1.1E-7, 0.000000001);
        Assert.assertEquals(counts.getEntry(1, 1), -2.2E-8, 0.000000001);
    }

    @Test
    public void testReadTargetNameOnlyFormattedFile() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println("SAMPLE2\tSAMPLE1\t" + TargetTableColumn.NAME.toString());
        writer.println("1.1\t2.2\ttgt_0");
        writer.println("-1.1E-7\t-2.2E-8\ttgt_1");
        writer.close();
        final ReadCountCollection subject = ReadCountCollectionUtils.parse(testFile);
        Assert.assertNotNull(subject);
        Assert.assertEquals(subject.columnNames(), Arrays.asList("SAMPLE2", "SAMPLE1"));
        Assert.assertEquals(subject.targets().stream().map(Target::getName).collect(Collectors.toList()), Arrays.asList("tgt_0", "tgt_1"));
        Assert.assertEquals(subject.targets().stream().map(Target::getInterval).collect(Collectors.toList()), Arrays.asList(null, null));
        Assert.assertEquals(subject.targets().size(), 2);
        final RealMatrix counts = subject.counts();
        Assert.assertEquals(counts.getRowDimension(), 2);
        Assert.assertEquals(counts.getColumnDimension(), 2);
        Assert.assertEquals(counts.getEntry(0, 0), 1.1, 0.0001);
        Assert.assertEquals(counts.getEntry(0, 1), 2.2, 0.0001);
        Assert.assertEquals(counts.getEntry(1, 0), -1.1E-7, 0.000000001);
        Assert.assertEquals(counts.getEntry(1, 1), -2.2E-8, 0.000000001);
    }

    @Test
    public void testReadIntervalsOnlyFile() throws IOException {
        final File targetFile = createTempFile();
        final PrintWriter targetWriter = new PrintWriter(targetFile);
        targetWriter.println(CONTIG_START_END_NAME);
        targetWriter.println("1\t100\t200\tTGT_0");
        targetWriter.println("2\t200\t300\tTGT_1");
        targetWriter.close();
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("contig\tSAMPLE3\tstop\tSAMPLE2\tstart");
        writer.println("1\t1.1\t200\t2.2\t100");
        writer.println("2\t-1.1E-7\t300\t-2.2E-8\t200");
        writer.close();

        final ReadCountCollection subject = ReadCountCollectionUtils.parse(testFile,
                new HashedListTargetCollection<>(TargetTableReader.readTargetFile(targetFile)), false);
        Assert.assertNotNull(subject);
        Assert.assertEquals(subject.columnNames(), Arrays.asList("SAMPLE3", "SAMPLE2"));
        Assert.assertEquals(subject.targets().stream().map(Target::getInterval).collect(Collectors.toList()), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("2", 200, 300)));
        Assert.assertEquals(subject.targets().stream().map(Target::getName).collect(Collectors.toList()), Arrays.asList("TGT_0", "TGT_1"));
        Assert.assertEquals(subject.targets().size(), 2);
        final RealMatrix counts = subject.counts();
        Assert.assertEquals(counts.getRowDimension(), 2);
        Assert.assertEquals(counts.getColumnDimension(), 2);
        Assert.assertEquals(counts.getEntry(0, 0), 1.1, 0.0001);
        Assert.assertEquals(counts.getEntry(0, 1), 2.2, 0.0001);
        Assert.assertEquals(counts.getEntry(1, 0), -1.1E-7, 0.000000001);
        Assert.assertEquals(counts.getEntry(1, 1), -2.2E-8, 0.000000001);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testOneTooFewValueInLine() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("contig\tSAMPLE3\tstop\tSAMPLE2\tstart");
        writer.println("1.1\t200\t2.2\t100");
        writer.println("2\t-1.1E-7\t300\t-2.2E-8\t200");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testOneTooManyValueInLine() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("contig\tSAMPLE3\tstop\tSAMPLE2\tstart");
        writer.println("1\t1.1\t200\t2.2\t100\t21.0");
        writer.println("2\t-1.1E-7\t300\t-2.2E-8\t200");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
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
                    ReadCountCollectionUtils.removeColumnsWithTooManyZeros(readCount, maxZeros, false, NULL_LOGGER);
                } catch (final UserException.BadInput ex) {
                    // expected.
                    continue;
                }
                Assert.fail("expects an exception");
            }
            final ReadCountCollection rc = ReadCountCollectionUtils.removeColumnsWithTooManyZeros(readCount, maxZeros, false, NULL_LOGGER);
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
                    ReadCountCollectionUtils.removeTargetsWithTooManyZeros(readCount, maxZeros, false, NULL_LOGGER);
                } catch (final UserException.BadInput ex) {
                    // expected.
                    continue;
                }
                Assert.fail("expects an exception");
            }
            final ReadCountCollection rc = ReadCountCollectionUtils.removeTargetsWithTooManyZeros(readCount, maxZeros, false, NULL_LOGGER);
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
        ReadCountCollectionUtils.truncateExtremeCounts(readCount, percentile, NULL_LOGGER);
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
        final ReadCountCollection result = ReadCountCollectionUtils.removeColumnsWithExtremeMedianCounts(readCount, percentile, NULL_LOGGER);
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
        ReadCountCollectionUtils.imputeZeroCountsAsTargetMedians(readCounts, NULL_LOGGER);
        final RealMatrix newCounts = readCounts.counts();
        Assert.assertEquals(newCounts.getColumnDimension(), expected[0].length);
        Assert.assertEquals(newCounts.getRowDimension(), expected.length);
        for (int i = 0; i < expected.length; i++) {
            for (int j = 0; j < expected[i].length; j++) {
                Assert.assertEquals(newCounts.getEntry(i, j), expected[i][j], "i,j == " + i + "," + j + " " + original[i][j]);
            }
        }
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
            final ReadCountCollection readCounts = new ReadCountCollection(targets, columnNames,
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
            final ReadCountCollection readCounts = new ReadCountCollection(targets, columnNames, new Array2DRowRealMatrix(counts, false));
            result.add(new Object[] { readCounts });
        }
        return result.toArray(new Object[result.size()][]);
    }

    private File createTempFile() throws IOException {
        final File result = File.createTempFile("file", ".test");
        result.deleteOnExit();
        return result;
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
