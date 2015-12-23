package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link CalculateCoverageStats}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class CalculateCoverageStatsIntegrationTest extends CommandLineProgramTest {

    public static final double EPSILON = 0.00001;

    private static final int[] TEST_SAMPLE_COUNTS = new int[] { 1, 2, 5, 10, 20, 100 };

    private static final int[] TEST_TARGET_COUNTS = new int[] { 1, 2, 3, 11, 31 };

    private static final Random TEST_RANDOM = new Random(13);

    @Override
    public String getTestedClassName() {
        return CalculateCoverageStats.class.getSimpleName();
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoData()
            throws IOException
    {
        final Object[] arguments = composeTestDataArguments(0, 0);
        @SuppressWarnings("unchecked")
        final List<String> sampleNames = (List<String>) arguments[0];
        @SuppressWarnings("unchecked")
        final List<Target> targets = (List<Target>) arguments[1];
        @SuppressWarnings("unchecked")
        final double[][] values = (double[][]) arguments[2];
        final File inputFile = createReadCountsFile(sampleNames, targets, values, true);
        runCommandLine(inputFile, true, true);
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testNoExistingInputFile() {
        runCommandLine(new File("bogus-file.tab"), true, true);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoSamples()
            throws IOException
    {
        final Object[] arguments = composeTestDataArguments(0, 10);
        @SuppressWarnings("unchecked")
        final List<String> sampleNames = (List<String>) arguments[0];
        @SuppressWarnings("unchecked")
        final List<Target> targets = (List<Target>) arguments[1];
        @SuppressWarnings("unchecked")
        final double[][] values = (double[][]) arguments[2];
        final File inputFile = createReadCountsFile(sampleNames, targets, values, true);
        runCommandLine(inputFile, true, true);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoTargets()
            throws IOException
    {
        final Object[] arguments = composeTestDataArguments(10, 0);
        @SuppressWarnings("unchecked")
        final List<String> sampleNames = (List<String>) arguments[0];
        @SuppressWarnings("unchecked")
        final List<Target> targets = (List<Target>) arguments[1];
        @SuppressWarnings("unchecked")
        final double[][] values = (double[][]) arguments[2];
        final File inputFile = createReadCountsFile(sampleNames, targets, values, true);
        runCommandLine(inputFile, true, true);
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testNoOutputFiles()
        throws IOException
    {
        final Object[] arguments = composeTestDataArguments(10, 10);
        @SuppressWarnings("unchecked")
        final List<String> sampleNames = (List<String>) arguments[0];
        @SuppressWarnings("unchecked")
        final List<Target> targets = (List<Target>) arguments[1];
        @SuppressWarnings("unchecked")
        final double[][] values = (double[][]) arguments[2];
        final File inputFile = createReadCountsFile(sampleNames, targets, values, true);
        runCommandLine(inputFile, false, false);
    }

    @Test(dataProvider = "testData")
    public void testWithIntervals(final List<String> sampleNames, final List<Target> targets, final double[][] values)
            throws IOException {
        final File inputFile = createReadCountsFile(sampleNames, targets, values, true);
        final Pair<File, File> outputFiles = runCommandLine(inputFile, true, true);
        final File targetOutputFile = outputFiles.getFirst();
        final File sampleOutputFile = outputFiles.getSecond();
        Assert.assertTrue(targetOutputFile.isFile());
        Assert.assertTrue(targetOutputFile.canRead());
        Assert.assertTrue(sampleOutputFile.isFile());
        Assert.assertTrue(sampleOutputFile.canRead());
        assertTargetFileContent(targetOutputFile, targets, values, true);
        assertSampleFileContent(sampleOutputFile, sampleNames, values);
        Assert.assertTrue(inputFile.delete());
        Assert.assertTrue(targetOutputFile.delete());
        Assert.assertTrue(sampleOutputFile.delete());
    }

    @Test(dataProvider = "testData")
    public void testOnlyTargetOutput(final List<String> sampleNames, final List<Target> targets, final double[][] values)
            throws IOException {
        final File inputFile = createReadCountsFile(sampleNames, targets, values, true);
        final Pair<File, File> outputFiles = runCommandLine(inputFile, true, false);
        final File targetOutputFile = outputFiles.getFirst();
        Assert.assertNull(outputFiles.getSecond());
        Assert.assertTrue(targetOutputFile.isFile());
        Assert.assertTrue(targetOutputFile.canRead());
        assertTargetFileContent(targetOutputFile, targets, values, true);
        Assert.assertTrue(inputFile.delete());
        Assert.assertTrue(targetOutputFile.delete());
    }

    @Test(dataProvider = "testData")
    public void testOnlySampleOutput(final List<String> sampleNames, final List<Target> targets, final double[][] values)
            throws IOException {
        final File inputFile = createReadCountsFile(sampleNames, targets, values, true);
        final Pair<File, File> outputFiles = runCommandLine(inputFile, false, true);
        Assert.assertNull(outputFiles.getFirst());
        final File sampleOutputFile = outputFiles.getSecond();
        Assert.assertTrue(sampleOutputFile.isFile());
        Assert.assertTrue(sampleOutputFile.canRead());
        assertSampleFileContent(sampleOutputFile, sampleNames, values);
        Assert.assertTrue(inputFile.delete());
        Assert.assertTrue(sampleOutputFile.delete());
    }

    @Test(dataProvider = "testData")
    public void testWithoutIntervals(final List<String> sampleNames, final List<Target> targets, final double[][] values)
            throws IOException {
        final File inputFile = createReadCountsFile(sampleNames, targets, values, false);
        final Pair<File, File> outputFiles = runCommandLine(inputFile, true, true);
        final File targetOutputFile = outputFiles.getFirst();
        final File sampleOutputFile = outputFiles.getSecond();
        Assert.assertTrue(targetOutputFile.isFile());
        Assert.assertTrue(targetOutputFile.canRead());
        Assert.assertTrue(sampleOutputFile.isFile());
        Assert.assertTrue(sampleOutputFile.canRead());
        assertTargetFileContent(targetOutputFile, targets, values, false);
        assertSampleFileContent(sampleOutputFile, sampleNames, values);
        Assert.assertTrue(inputFile.delete());
        Assert.assertTrue(targetOutputFile.delete());
        Assert.assertTrue(sampleOutputFile.delete());
    }

    private void assertSampleFileContent(final File sampleOutputFile, final List<String> sampleNames,
                                         final double[][] values) throws IOException {
        Assert.assertNotNull(sampleOutputFile);
        try (final SampleCoverageStatsReader reader = new SampleCoverageStatsReader(sampleOutputFile)) {
            final List<SampleCoverageStats> stats = reader.stream().collect(Collectors.toList());
            Assert.assertEquals(stats.size(), sampleNames.size());
            for (int i = 0; i < stats.size(); i++) {
                Assert.assertEquals(stats.get(i).sample, sampleNames.get(i));
                final double observedMean = stats.get(i).mean;
                final double observedVariance = stats.get(i).variance;
                final int index = i;
                final double[] sampleValues = IntStream.range(0, values.length)
                        .mapToDouble(j -> values[j][index])
                        .toArray();
                final double expectedMean = MathUtils.mean(sampleValues, 0, sampleValues.length);
                Assert.assertEquals(observedMean, expectedMean, EPSILON);
                final double expectedVariance = sampleValues.length == 1 ? 0 : DoubleStream.of(sampleValues)
                        .map(d -> Math.pow(d - expectedMean , 2))
                        .sum() / (sampleValues.length - 1);
                Assert.assertTrue(Double.isNaN(observedVariance) == Double.isNaN(expectedVariance) || Math.abs(observedVariance - expectedVariance) <= EPSILON);
            }
        }
    }

    private void assertTargetFileContent(final File targetOutputFile, final List<Target> targets, final double[][] values, final boolean withIntervals) throws IOException {
        Assert.assertNotNull(targetOutputFile);
        try (final TargetCoverageStatsReader reader = new TargetCoverageStatsReader(targetOutputFile)) {
            if (withIntervals) {
                Assert.assertTrue(reader.columns().containsAll(TargetTableColumns.MANDATORY_COLUMN_NAME_ARRAY));
            } else {
                Assert.assertTrue(reader.columns().contains(TargetTableColumns.NAME.toString()));
            }
            Assert.assertTrue(reader.columns().containsAll(TargetCoverageStats.MEAN_COLUMN_NAME, TargetCoverageStats.VARIANCE_COLUMN_NAME));
            final List<TargetCoverageStats> stats = reader.stream().collect(Collectors.toList());
            Assert.assertEquals(stats.size(), targets.size());
            for (int i = 0; i < stats.size(); i++) {
                Assert.assertNotNull(stats.get(i).target);
                Assert.assertEquals(stats.get(i).target, targets.get(i));
                if (withIntervals) {
                    Assert.assertNotNull(stats.get(i).target.getInterval());
                    Assert.assertEquals(stats.get(i).target.getInterval(), targets.get(i).getInterval());
                }
                final double observedMean = stats.get(i).mean;
                final double expectedMean = MathUtils.mean(values[i], 0, values[i].length);
                Assert.assertEquals(observedMean, expectedMean, 0.00001);
                final double observedVariance = stats.get(i).variance;
                final double expectedVariance = values[i].length == 1 ? 0 : DoubleStream.of(values[i])
                        .map(d -> Math.pow(d - expectedMean , 2))
                        .sum() / (values[i].length - 1);
                Assert.assertTrue(Double.isNaN(observedVariance) == Double.isNaN(expectedVariance) ||  Math.abs(observedVariance - expectedVariance) <= EPSILON);
            }
        }
    }

    private Pair<File,File> runCommandLine(final File inputFile, final boolean generateTargetOutput, final boolean generateSampleOutput) {
        final File sampleOutput = generateSampleOutput ? createTempFile("ccsi-sample-out", ".tsv") : null;
        final File targetOutput = generateTargetOutput ? createTempFile("ccsi-target-out", ".tsv") : null;
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(inputFile.getPath());
        if (sampleOutput != null) {
            arguments.add("-" + CalculateCoverageStats.SAMPLE_OUTPUT_FILE_SHORT_NAME);
            arguments.add(sampleOutput.getPath());
        }
        if (targetOutput != null) {
            arguments.add("-" + CalculateCoverageStats.TARGET_OUTPUT_FILE_SHORT_NAME);
            arguments.add(targetOutput.getPath());
        }
        runCommandLine(arguments);
        return new Pair<>(targetOutput, sampleOutput);
    }

    private File createReadCountsFile(final List<String> sampleNames, final List<Target> targets, final double[][] values, final boolean withIntervals) throws IOException {
        final File result = createTempFile("ccsi-input", ".tab");
        if (values.length == 0) {
            ReadCountCollectionUtils.writerWithIntervals(new FileWriter(result), sampleNames).close();
        } else if (sampleNames.size() == 0) {
            try (final TargetTableWriter writer = new TargetTableWriter(result)) {
                writer.writeAllRecords(targets);
            }
        } else {
            final ReadCountCollection coverage = new ReadCountCollection(
                    SetUniqueList.setUniqueList(withIntervals ? new ArrayList<>(targets) : targets.stream()
                            .map(t -> new Target(t.getName())).collect(Collectors.toList())),
                    SetUniqueList.setUniqueList(new ArrayList<>(sampleNames)),
                    new Array2DRowRealMatrix(values));
            ReadCountCollectionUtils.write(result, coverage);
        }
        return result;
    }

    @DataProvider(name = "testData")
    public Object[][] testData() {
        final List<Object[]> result = new ArrayList<>(TEST_SAMPLE_COUNTS.length * TEST_TARGET_COUNTS.length);
        for (final int sampleCount : TEST_SAMPLE_COUNTS) {
            for (final int targetCount : TEST_TARGET_COUNTS)  {
                final Object[] arguments = composeTestDataArguments(sampleCount, targetCount);
                result.add(arguments);
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    private Object[] composeTestDataArguments(final int sampleCount, final int targetCount) {
        final List<String> sampleNames = createSampleNames(sampleCount);
        final List<Target> targets = createTargets(targetCount);
        final double[][] values = new double[targetCount][sampleCount];
        for (int i = 0; i < targetCount; i++) {
            for (int j = 0; j < sampleCount; j++) {
                values[i][j] = TEST_RANDOM.nextGaussian() * 100 + 200;
            }
        }
        return new Object[] {sampleNames, targets, values};
    }

    private List<Target> createTargets(final int targetCount) {
        final List<Target> result = new ArrayList<>(targetCount);
        for (int i = 0; i < targetCount; i++) {
            result.add(new Target("TARGET_" + (i + 1), new SimpleInterval("chr1", (100 * i + 1), (100 * i + 51))));
        }
        return result;
    }

    private List<String> createSampleNames(final int sampleCount) {
        final List<String> result = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++) {
            result.add("SAMPLE_" + (i + 1));
        }
        return result;
    }
}
