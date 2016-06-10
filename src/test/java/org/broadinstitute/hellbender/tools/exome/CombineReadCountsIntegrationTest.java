package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.fakedata.SimulatedTargets;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link CombineReadCounts}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class CombineReadCountsIntegrationTest extends CommandLineProgramTest {

    private final int[] TEST_SAMPLE_COUNT = new int[] { 1, 2, 10, 100, 1000 };

    private final int TEST_TOTAL_COUNT = 10000;

    private final String[] TEST_CHR = new String[] { "1", "2", "10", "X"};

    private final int[] TEST_CHR_LEN = new int[] {10000000, 10000000, 10000000};
    private final int TEST_GENOME_LEN = IntStream.of(TEST_CHR_LEN).sum();

    private final int TEST_TARGET_LEN = 111;

    @Test(expectedExceptions = UserException.class)
    public void testEmptyInputFileListProvided() throws Exception {
        @SuppressWarnings("serial")
        final List<Target> phonyTargets = SimulatedTargets.phonyTargets(3);
        final List<File> inputFiles = Collections.emptyList();
        final File targetFile = createTargetFile(phonyTargets);
        final File inputListFile = createInputListFile(inputFiles);
        try {
            runTool(targetFile, inputFiles, inputListFile);
        } catch (final Exception ex) {
            throw ex;
        } finally {
            targetFile.delete();
            inputListFile.delete();
        }
    }

    @Test(expectedExceptions = UserException.class)
    public void testNoInputFileProvided() throws Exception {
        @SuppressWarnings("serial")
        final List<Target> phonyTargets = SimulatedTargets.phonyTargets(3);
        final List<File> inputFiles = Collections.emptyList();
        final File targetFile = createTargetFile(phonyTargets);
        try {
            runTool(targetFile, inputFiles, null);
        } catch (final Exception ex) {
            throw ex;
        } finally {
            targetFile.delete();
        }
    }

    @Test(dataProvider = "testData")
    public void testArbitraryTargetOrder(final List<Target> targets, final List<String> sampleNames, final double[][] counts) throws IOException {
        final List<Target> finalTargets = new ArrayList<>(targets);
        Collections.shuffle(finalTargets, new Random(13));
        final List<File> inputFiles = createInputCountFiles(finalTargets, sampleNames, counts, true, true);
        final File targetFile = createTargetFile(targets);
        try {
            runTool(targetFile, inputFiles, null);
        } catch (final Exception ex) {
            throw ex;
        } finally {
            targetFile.delete();
        }
    }

    @Test(dataProvider = "testData", expectedExceptions = UserException.BadInput.class)
    public void testRepeatedInputSamples(final List<Target> targets, final List<String> sampleNames, final double[][] counts) throws IOException {
        final List<File> inputFiles = createInputCountFiles(targets, sampleNames, counts, true, true);
        inputFiles.add(inputFiles.get(inputFiles.size() >> 1));
        final File targetFile = createTargetFile(targets);
        try {
            runTool(targetFile, inputFiles, null);
        } catch (final Exception ex) {
            throw ex;
        } finally {
            targetFile.delete();
        }
    }

    @Test(dataProvider = "testData", expectedExceptions = UserException.BadInput.class)
    public void testMissingTargets(final List<Target> targets, final List<String> sampleNames, final double[][] counts) throws IOException {
        final List<Target> finalTargets = new ArrayList<>(targets);
        finalTargets.remove(finalTargets.size() >> 1);
        final List<File> inputFiles = createInputCountFiles(finalTargets, sampleNames, counts, true, true);
        final File targetFile = createTargetFile(targets);
        try {
            runTool(targetFile, inputFiles, null);
        } catch (final Exception ex) {
            throw ex;
        } finally {
            targetFile.delete();
        }
    }

    @Test(dataProvider="testData")
    public void testMixedSampleInputApproachFullTargetInfo(final List<Target> targets, final List<String> sampleNames, final double[][] counts) throws IOException {
        final List<File> inputFiles = createInputCountFiles(targets, sampleNames, counts, true, true);
        final List<File> listedInputFiles = inputFiles.subList(0, (inputFiles.size() >> 1) + 1);
        final List<File> directInputFiles = inputFiles.subList((inputFiles.size() >> 1) + 1, inputFiles.size());
        final File targetFile = createTargetFile(targets);
        final File inputFileList = createInputListFile(listedInputFiles);
        final File output = runTool(targetFile, directInputFiles, inputFileList);
        // Make sure that the input files are not removed by the tool.
        Assert.assertFalse(inputFiles.stream().anyMatch(f -> !f.canRead()));
        inputFiles.forEach(File::delete);
        targetFile.delete();
        Assert.assertTrue(output.canRead());
        assertOutputContents(output, targets, sampleNames, counts);
        // do the testing here:
        output.delete();
    }

    @Test(dataProvider="testData")
    public void testOneSampleOneFileFullTargetInfo(final List<Target> targets, final List<String> sampleNames, final double[][] counts) throws IOException {
        final List<File> inputFiles = createInputCountFiles(targets, sampleNames, counts, true, true);
        final File targetFile = createTargetFile(targets);
        final File output = runTool(targetFile, inputFiles, null);
        // Make sure that the input files are not removed by the tool.
        Assert.assertFalse(inputFiles.stream().anyMatch(f -> !f.canRead()));
        inputFiles.forEach(File::delete);
        targetFile.delete();
        Assert.assertTrue(output.canRead());
        assertOutputContents(output, targets, sampleNames, counts);
        // do the testing here:
        output.delete();
    }

    @Test(dataProvider="testData")
    public void testOneSampleOneFileOnlyNames(final List<Target> targets, final List<String> sampleNames, final double[][] counts) throws IOException {
        final List<File> inputFiles = createInputCountFiles(targets, sampleNames, counts, true, false);
        final File targetFile = createTargetFile(targets);
        final File output = runTool(targetFile, inputFiles, null);
        // Make sure that the input files are not removed by the tool.
        Assert.assertFalse(inputFiles.stream().anyMatch(f -> !f.canRead()));
        inputFiles.forEach(File::delete);
        targetFile.delete();
        Assert.assertTrue(output.canRead());
        assertOutputContents(output, targets, sampleNames, counts);
        // do the testing here:
        output.delete();
    }

    @Test(dataProvider="testData")
    public void testOneSampleOneFileOnlyCoordinates(final List<Target> targets, final List<String> sampleNames, final double[][] counts) throws IOException {
        final List<File> inputFiles = createInputCountFiles(targets, sampleNames, counts, false, true);
        final File targetFile = createTargetFile(targets);
        final File output = runTool(targetFile, inputFiles, null);
        // Make sure that the input files are not removed by the tool.
        Assert.assertFalse(inputFiles.stream().anyMatch(f -> !f.canRead()));
        inputFiles.forEach(File::delete);
        targetFile.delete();
        Assert.assertTrue(output.canRead());
        assertOutputContents(output, targets, sampleNames, counts);
        // do the testing here:
        output.delete();
    }

    @Test(dataProvider="testData")
    public void testAllSamplesInOneFileFullTargetInfoNoTargetsFile(final List<Target> targets, final List<String> sampleNames, final double[][] counts) throws IOException {
        final List<File> inputFiles = createInputCountFiles(targets, sampleNames, counts, true, true);
        final File inputListFile = createInputListFile(inputFiles);
        final File output = runTool(null, Collections.emptyList(), inputListFile);
        // Make sure that the input files are not removed by the tool.
        Assert.assertFalse(inputFiles.stream().anyMatch(f -> !f.canRead()));
        inputFiles.forEach(File::delete);
        Assert.assertTrue(output.canRead());
        assertOutputContents(output, targets, sampleNames, counts);
        // do the testing here:
        output.delete();
    }

    @Test(dataProvider="testData")
    public void testAllSamplesInOneFileFullTargetInfo(final List<Target> targets, final List<String> sampleNames, final double[][] counts) throws IOException {
        final List<File> inputFiles = createInputCountFiles(targets, sampleNames, counts, true, true);
        final File targetFile = createTargetFile(targets);
        final File inputListFile = createInputListFile(inputFiles);
        final File output = runTool(targetFile, Collections.emptyList(), inputListFile);
        // Make sure that the input files are not removed by the tool.
        Assert.assertFalse(inputFiles.stream().anyMatch(f -> !f.canRead()));
        inputFiles.forEach(File::delete);
        targetFile.delete();
        Assert.assertTrue(output.canRead());
        assertOutputContents(output, targets, sampleNames, counts);
        // do the testing here:
        output.delete();
    }

    private File createInputListFile(final List<File> inputFiles) throws IOException {
        final File result = createTempFile("inputs", ".list");
        final PrintWriter writer = new PrintWriter(new FileWriter(result));
        inputFiles.forEach(f -> writer.println(f.getAbsoluteFile()));
        writer.close();
        return result;
    }

    private void assertOutputContents(final File output, final List<Target> targets, final List<String> sampleNames, final double[][] counts) throws IOException {
        Assert.assertTrue(output.length() > 10);
        final ReadCountCollection readCounts = ReadCountCollectionUtils.parse(output);
        final List<String> sortedSampleNames = new ArrayList<>(sampleNames);
        Collections.sort(sortedSampleNames);
        Assert.assertEquals(readCounts.targets(), targets);
        Assert.assertEquals(readCounts.columnNames(), sortedSampleNames);
        final RealMatrix actualCounts = readCounts.counts();
        for (int i = 0; i < actualCounts.getColumnDimension(); i++) {
            final int actualColumn = readCounts.columnNames().indexOf(sampleNames.get(i));
            for (int j = 0; j < actualCounts.getRowDimension(); j++) {
                Assert.assertEquals(actualCounts.getEntry(j, actualColumn), counts[i][j], 0.00001);
            }
        }
    }

    private File runTool(final File targetFile, final List<File> inputFiles, final File inputFileList) {
        final List<String> args = new ArrayList<>();
        if (targetFile != null) {
            args.add("-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME);
            args.add(targetFile.getAbsolutePath());
        }
        if (inputFileList != null) {
            args.add("-" + CombineReadCounts.READ_COUNT_FILE_LIST_SHORT_NAME);
            args.add(inputFileList.getAbsolutePath());
        }
        for (final File inputFile : inputFiles) {
            args.add("-" + CombineReadCounts.READ_COUNT_FILES_SHORT_NAME);
            args.add(inputFile.getAbsolutePath());
        }
        final File outputFile = createTempFile("output", ".tab");
        args.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        args.add(outputFile.getAbsolutePath());
        args.add("-" + CombineReadCounts.MAX_GROUP_SIZE_SHORT_NAME);
        args.add("7");
        runCommandLine(args);
        return outputFile;
    }

    private List<File> createInputCountFiles(final List<Target> targets, final List<String> sampleNames, final double[][] counts, final boolean withTargetName, final boolean withTargetCoordinates) throws IOException {
        final List<File> inputFiles = sampleNames.stream()
                .map(s -> BaseTest.createTempFile(s, ".tab"))
                .collect(Collectors.toList());
        for (int i = 0; i < sampleNames.size(); i++) {
            createCountFile(inputFiles.get(i), targets, sampleNames.get(i), counts[i], withTargetName, withTargetCoordinates);
        }
        return inputFiles;
    }

    private File createTargetFile(final List<Target> targets) throws IOException {
        final File result = BaseTest.createTempFile("targets",".tab");
        final TableWriter<Target> writer = TableUtils.writer(result, new TableColumnCollection(
                        TargetTableColumn.CONTIG,
                        TargetTableColumn.START,
                        TargetTableColumn.END,
                        TargetTableColumn.NAME),
                (t,dataLine) -> {
                    final SimpleInterval interval = t.getInterval();
                    dataLine.append(interval.getContig())
                            .append(interval.getStart()).append(interval.getEnd())
                            .append(t.getName());
                }
        );
        for (final Target target : targets)
            writer.writeRecord(target);
        writer.close();
        return result;
    }

    private void createCountFile(final File output, List<Target> targets, final String sample, final double[] count, final boolean useName, final boolean useCoordinates) throws IOException {
        final List<String> columnNames = new ArrayList<>();
        if (useCoordinates) {
            columnNames.add(TargetTableColumn.CONTIG.toString());
            columnNames.add(TargetTableColumn.START.toString());
            columnNames.add(TargetTableColumn.END.toString());
        }
        if (useName) {
            columnNames.add(TargetTableColumn.NAME.toString());
        }
        columnNames.add(sample);
        final TableColumnCollection columns = new TableColumnCollection(columnNames);
        final TableWriter<Integer> writer = TableUtils.writer(output, columns, (i,d) -> {
            final Target target = targets.get(i);
            final SimpleInterval interval = target.getInterval();
            if (useCoordinates)
                d.append(interval.getContig()).append(interval.getStart()).append(interval.getEnd());
            if (useName)
                d.append(target.getName());
            d.append(count[i]);
        });
        for (int i = 0; i < targets.size(); i++)
            writer.writeRecord(i);
        writer.close();
    }


    @DataProvider(name="testData")
    public Object[][] testData() {
        final List<Object[]> result = new ArrayList<>();
        for (final int sampleCount : TEST_SAMPLE_COUNT) {
            final int targetCount = TEST_TOTAL_COUNT / sampleCount;
            final List<Target> targets = createTargetList(targetCount);
            final Random rdn = new Random(11 * sampleCount + targets.size() * 31);
            final List<String> sampleNames = createSampleNames(sampleCount);
            final double[][] counts = createCounts(sampleCount, targets.size(), rdn);
            result.add(new Object[] { targets, sampleNames, counts });
        }
        return result.toArray(new Object[result.size()][]);
    }

    private double[][] createCounts(final int sampleCount, final int targetCount, final Random rdn) {
        final double[][] result = new double[sampleCount][targetCount];
        for (int i = 0; i < sampleCount; i++)
            for (int j = 0; j < targetCount; j++)
                result[i][j] = rdn.nextDouble();
        return result;
    }

    private List<String> createSampleNames(final int sampleCount) {
        final int digits = (int) Math.ceil(Math.log10(sampleCount));
        return IntStream.range(0, sampleCount)
                .mapToObj(i -> String.format("sample_%d" +  digits + "d",i))
                .collect(Collectors.toCollection(() -> {
                    @SuppressWarnings("serial")
                    final List<String> result = new ArrayList<String>(sampleCount) {
                        @Override
                        public String toString() {
                            final String superResult = super.toString();
                            return (superResult.length() > 80) ? superResult.substring(0, 76) + " ...]" : superResult;
                        }
                    };
                    return result;
                }));
    }

    private List<Target> createTargetList(final int targetCount) {
        final double geneLength = TEST_GENOME_LEN / targetCount;
        final int targetLength = (int) Math.min(geneLength / 2.0, TEST_TARGET_LEN);
        @SuppressWarnings("serial")
        final List<Target> result = new ArrayList<Target>(targetCount) {
            @Override
            public String toString() {
                final String superResult = super.toString();
                return (superResult.length() > 80) ? superResult.substring(0, 76) + " ...]" : superResult;
            }
        };
        for (int i = 0; i < targetCount; i++) {
            final String name = "target_" + i;
            final double geneStart = geneLength * i;
            final int absoluteStart = (int) Math.round(geneStart + (geneLength - targetLength) / 2.0);
            int chr;
            int relativeStart;
            for (chr = 0, relativeStart = absoluteStart; chr < TEST_CHR.length; relativeStart -= TEST_CHR_LEN[chr++]) {
                if (relativeStart < TEST_CHR_LEN[chr]) {
                    break;
                }
            }
            final int end = Math.min(relativeStart + targetLength, TEST_CHR_LEN[chr] - 1);
            final SimpleInterval interval = new SimpleInterval(TEST_CHR[chr], relativeStart + 1, end);
            final Target target = new Target(name, interval);
            result.add(target);
        }

        return result;
    }
}
