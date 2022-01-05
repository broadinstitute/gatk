package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.GATKIntervalListScatterMode;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.util.IntervalList.IntervalListScatterMode;

import java.io.File;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;


/**
 * Created by David Benjamin on 4/25/17.
 */
public class    SplitIntervalsIntegrationTest extends CommandLineProgramTest {

    private static final GenomeLocParser GLP = new GenomeLocParser(ReferenceDataSource.of(IOUtils.getPath(b37Reference)).getSequenceDictionary());
    private static final DecimalFormat formatter = new DecimalFormat("0000");

    @Test
    public void testOneInterval() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval("20:1000000-2000000")
                .addReference(b37Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .addOutput(outputDir);

        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);
        checkIntervalSizes(scatterCount, outputDir, 1000000, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);
    }

    @Test
    public void testOneIntervalWithWeights() {
        final int scatterCount = 5;

        final File outputDir = new File("out");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval("1:1-100")
                .addReference(b37Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(SplitIntervals.WEIGHTS_BED_FILE_FULL_NAME, "intervals.bed")
                .add(SplitIntervals.SUBDIVISION_MODE_lONG_NAME, GATKIntervalListScatterMode.INTERVAL_SUBDIVISION_BY_WEIGHT)
                .addOutput(outputDir);

        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);
        checkIntervalSizes(scatterCount, outputDir, 100, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);
    }


    @Test
    public void testOneIntervalFilenamePrefix() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");
        final String prefix = "scattered-with-weird-prefix-";
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval("20:1000000-2000000")
                .addReference(b37Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(SplitIntervals.INTERVAL_FILE_PREFIX_FULL_NAME, prefix)
                .addOutput(outputDir);

        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, prefix, SplitIntervals.DEFAULT_EXTENSION);
        checkIntervalSizes(scatterCount, outputDir, 1000000,  SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, prefix, SplitIntervals.DEFAULT_EXTENSION);
    }


    @Test
    public void testOneIntervalAlternateExtension() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");
        final String extension = "-scattered.with.a.weird.extension";
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval("20:1000000-2000000")
                .addReference(b37Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(SplitIntervals.INTERVAL_FILE_EXTENSION_FULL_NAME, extension)
                .addOutput(outputDir);

        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, extension);
        checkIntervalSizes(scatterCount, outputDir, 1000000, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, extension);
    }

    @Test
    public void testSingleScatter() {
        final int scatterCount = 1;
        final File outputDir = createTempDir("output");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval("20:1000000-2000000")
                .addReference(b37Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .addOutput(outputDir);

        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);
        checkIntervalSizes(scatterCount, outputDir, 1000000, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);

    }

    @Test
    public void testTwoIntervals() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval("20:1000000-2000000")
                .addInterval("20:3000000-4000000")
                .addReference(b37Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .addOutput(outputDir);

        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);
        checkIntervalSizes(scatterCount, outputDir, 2000000, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);

    }

    @Test
    public void testTwoIntervalsWithDigits() {
        final int scatterCount = 5;
        final int numDigits = 10;
        final File outputDir = createTempDir("output");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval("20:1000000-2000000")
                .addInterval("20:3000000-4000000")
                .addReference(b37Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(SplitIntervals.INTERVAL_NUMBER_OF_DIGITS_FULL_NAME, numDigits)
                .addOutput(outputDir);

        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, numDigits, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);
        checkIntervalSizes(scatterCount, outputDir, 2000000, numDigits, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);

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
        verifyScatteredFilesExist(scatterCount + 1, outputDir, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);
        checkTotalSize(scatterCount + 1, outputDir, 2000002, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);
    }

    @Test
    public void testNoIntervals() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(b37Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .addOutput(outputDir);

        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);
        final long totalLengthInRef = GLP.getSequenceDictionary().getSequences().stream().mapToLong(SAMSequenceRecord::getSequenceLength).sum();
        checkIntervalSizes(scatterCount, outputDir, totalLengthInRef, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);
    }

    @Test
    public void testAsInWGSJointCalling() {
        //these intervals are small and adjacent -- we want to distribute them among multiple interval lists, but the default behavior is to merge
        final File wgsCallingIntervals = new File(publicTestDir + "hg38.handcurated.noCentromeres.noTelomeres.interval_list");
        final int scatterCount = 100;
        final File outputDir = createTempDir("output");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addIntervals(wgsCallingIntervals)
                .addReference(hg38Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(SplitIntervals.SUBDIVISION_MODE_lONG_NAME, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW.toString())
                .add(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputDir);

        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir, SplitIntervals.DEFAULT_NUMBER_OF_DIGITS, SplitIntervals.DEFAULT_PREFIX, SplitIntervals.DEFAULT_EXTENSION);

        final File outputDir2 = createTempDir("output2");

        final ArgumentsBuilder args2 = new ArgumentsBuilder()
                .addIntervals(wgsCallingIntervals)
                .addReference(hg38Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(SplitIntervals.SUBDIVISION_MODE_lONG_NAME, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW.toString())
                .addOutput(outputDir2);

        runCommandLine(args2);
        Assert.assertTrue(outputDir2.listFiles().length < scatterCount);
    }

    //generates the files to look for given a scatter count, directory and extension
    private static Stream<File> getExpectedScatteredFiles(final int scatterCount, final File outputDir, final int numDigits, String prefix, String extension) {
        final int maxNumberOfPlaces = (int)Math.max(Math.floor(Math.log10(scatterCount-1))+1, numDigits);
        final String fileIndexFormatString = "%0" + maxNumberOfPlaces + "d";
        return IntStream.range(0, scatterCount).mapToObj(n -> new File(outputDir, prefix + String.format(fileIndexFormatString, n) + extension));
    }

    private static void verifyScatteredFilesExist(final int scatterCount, final File outputDir, final int numDigits, String prefix, String extension) {
        getExpectedScatteredFiles(scatterCount, outputDir, numDigits, prefix, extension).forEach(f -> Assert.assertTrue(f.exists()));
        Assert.assertFalse(new File(outputDir, prefix + formatter.format(scatterCount) + extension).exists());
    }

    private static List<SimpleInterval> readIntervals(final File intervalsFile) {
        return Utils.stream(IntervalList.fromFile(intervalsFile))
                .map(SimpleInterval::new)
                .collect(Collectors.toList());
    }

    private static void checkIntervalSizes(final int scatterCount, final File outputDir, final long expectedTotalLength, final int numDigits, String prefix, String extension) {
        final long splitLength = expectedTotalLength / scatterCount;
        getExpectedScatteredFiles(scatterCount, outputDir, numDigits, prefix, extension)
                .forEach(f -> Assert.assertEquals(readIntervals(f).stream()
                        .mapToLong(SimpleInterval::size)
                        .sum(), splitLength, 100));
    }

    private static void checkTotalSize(final int scatterCount, final File outputDir, final long expectedTotalLength, final int numDigits, String prefix, String extension) {
        final long totalLength = getExpectedScatteredFiles(scatterCount, outputDir, numDigits, prefix, extension)
                .mapToLong(f -> readIntervals(f).stream().mapToLong(SimpleInterval::size).sum())
                .sum();
        Assert.assertEquals(totalLength, expectedTotalLength);
    }

}



