package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.util.IntervalList.IntervalListScatterMode;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;


/**
 * Created by David Benjamin on 4/25/17.
 */
public class SplitIntervalsIntegrationTest extends CommandLineProgramTest {

    private static final Path REFERENCE = Paths.get(b37_reference_20_21);
    private static final GenomeLocParser GLP = new GenomeLocParser(ReferenceDataSource.of(REFERENCE).getSequenceDictionary());
    private static final DecimalFormat formatter = new DecimalFormat("0000");

    @Test
    public void testOneInterval() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");
        final String[] args = {
                "-L", "20:1000000-2000000",
                "-R", REFERENCE.toAbsolutePath().toString(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount),
                "-O", outputDir.getAbsolutePath()
        };
        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_EXTENSION);
        checkIntervalSizes(scatterCount, outputDir, 1000000, SplitIntervals.DEFAULT_EXTENSION);
    }

    @Test
    public void testOneIntervalAlternateExtension() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");
        final String extension = "-scattered.with.a.wierd.extension";
        final String[] args = {
                "-L", "20:1000000-2000000",
                "-R", REFERENCE.toAbsolutePath().toString(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount),
                "-O", outputDir.getAbsolutePath(),
                "--extension", extension
        };
        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, extension);
        checkIntervalSizes(scatterCount, outputDir, 1000000, extension);
    }

    @Test
    public void testSingleScatter() {
        final int scatterCount = 1;
        final File outputDir = createTempDir("output");
        final String[] args = {
                "-L", "20:1000000-2000000",
                "-R", REFERENCE.toAbsolutePath().toString(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount),
                "-O", outputDir.getAbsolutePath()
        };
        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_EXTENSION);
        checkIntervalSizes(scatterCount, outputDir, 1000000, SplitIntervals.DEFAULT_EXTENSION);

    }

    @Test
    public void testTwoIntervals() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");
        final String[] args = {
                "-L", "20:1000000-2000000",
                "-L", "20:3000000-4000000",
                "-R", REFERENCE.toAbsolutePath().toString(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount),
                "-O", outputDir.getAbsolutePath()
        };
        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_EXTENSION);
        checkIntervalSizes(scatterCount, outputDir, 2000000, SplitIntervals.DEFAULT_EXTENSION);

    }

    @Test
    public void testDontMixContigs() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval("20:1000000-2000000")
                .addInterval("21:1000000-2000000")
                .addReference(b37Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(SplitIntervals.DONT_MIX_CONTIGS_LONG_NAME, true)
                .addOutput(outputDir);
        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount + 1, outputDir, SplitIntervals.DEFAULT_EXTENSION);
        checkTotalSize(scatterCount + 1, outputDir, 2000002, SplitIntervals.DEFAULT_EXTENSION);
    }

    @Test
    public void testNoIntervals() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");
        final String[] args = {
                "-R", REFERENCE.toAbsolutePath().toString(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount),
                "-O", outputDir.getAbsolutePath()
        };
        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_EXTENSION);
        final int totalLengthInRef = GLP.getSequenceDictionary().getSequences().stream().mapToInt(SAMSequenceRecord::getSequenceLength).sum();
        checkIntervalSizes(scatterCount, outputDir, totalLengthInRef, SplitIntervals.DEFAULT_EXTENSION);

    }

    @Test
    public void testAsInWGSJointCalling() {
        //these intervals are small and adjacent -- we want to distribute them among multiple interval lists, but the default behavior is to merge
        final File wgsCallingIntervals = new File(publicTestDir + "hg38.handcurated.noCentromeres.noTelomeres.interval_list");
        final int scatterCount = 100;
        final File outputDir = createTempDir("output");
        final String[] args = {
                "-R", hg38Reference,
                "-L", wgsCallingIntervals.getAbsolutePath(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount),
                "--" + SplitIntervals.SUBDIVISION_MODE_lONG_NAME, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW.toString(),
                "-imr", IntervalMergingRule.OVERLAPPING_ONLY.toString(),
                "-O", outputDir.getAbsolutePath()
        };
        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_EXTENSION);

        final File outputDir2 = createTempDir("output2");
        final String[] args2 = {
                "-R", hg38Reference,
                "-L", wgsCallingIntervals.getAbsolutePath(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount),
                "--" + SplitIntervals.SUBDIVISION_MODE_lONG_NAME, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW.toString(),
                "-O", outputDir2.getAbsolutePath()
        };
        runCommandLine(args2);
        Assert.assertTrue(outputDir2.listFiles().length < scatterCount);
    }

    //generates the files to look for given a scatter count, directory and extension
    private static Stream<File> getExpectedScatteredFiles(final int scatterCount, final File outputDir, String extension) {
        return IntStream.range(0, scatterCount).mapToObj(n -> new File(outputDir, formatter.format(n) + extension));
    }

    private static void verifyScatteredFilesExist(final int scatterCount, final File outputDir, String extension) {
        getExpectedScatteredFiles(scatterCount, outputDir, extension).forEach(f -> Assert.assertTrue(f.exists()));
        Assert.assertFalse(new File(outputDir, formatter.format(scatterCount) + extension).exists());
    }

    private static List<SimpleInterval> readIntervals(final File intervalsFile) {
        return Utils.stream(IntervalList.fromFile(intervalsFile)).map(SimpleInterval::new).collect(Collectors.toList());
    }

    private static void checkIntervalSizes(final int scatterCount, final File outputDir, final int expectedTotalLength, String extension) {
        final int splitLength = expectedTotalLength / scatterCount;
        getExpectedScatteredFiles(scatterCount, outputDir, extension).forEach(f -> Assert.assertEquals(readIntervals(f).stream().mapToInt(SimpleInterval::size).sum(), splitLength, 100));
    }

    private static void checkTotalSize(final int scatterCount, final File outputDir, final int expectedTotalLength, String extension) {
        final int totalLength = getExpectedScatteredFiles(scatterCount, outputDir, extension)
                .mapToInt(f -> readIntervals(f).stream().mapToInt(SimpleInterval::size).sum())
                .sum();
        Assert.assertEquals(totalLength, expectedTotalLength);
    }

}