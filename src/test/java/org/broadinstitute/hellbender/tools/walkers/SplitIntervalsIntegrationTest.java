package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceRecord;
import junit.framework.Assert;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.testng.Assert.*;

/**
 * Created by David Benjamin on 4/25/17.
 */
public class SplitIntervalsIntegrationTest extends CommandLineProgramTest {

    private static final File REFERENCE = new File(b37_reference_20_21);
    private static final GenomeLocParser GLP = new GenomeLocParser(ReferenceDataSource.of(REFERENCE).getSequenceDictionary());

    @Test
    public void testOneInterval() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");
        final String[] args = {
                "-L", "20:1000000-2000000",
                "-R", REFERENCE.getAbsolutePath(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount),
                "-O", outputDir.getAbsolutePath()
        };
        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir);
        checkIntervalSizes(scatterCount, outputDir, 1000000);
    }

    @Test
    public void testSingleScatter() {
        final int scatterCount = 1;
        final File outputDir = createTempDir("output");
        final String[] args = {
                "-L", "20:1000000-2000000",
                "-R", REFERENCE.getAbsolutePath(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount),
                "-O", outputDir.getAbsolutePath()
        };
        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir);
        checkIntervalSizes(scatterCount, outputDir, 1000000);

    }

    @Test
    public void testTwoIntervals() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");
        final String[] args = {
                "-L", "20:1000000-2000000",
                "-L", "20:3000000-4000000",
                "-R", REFERENCE.getAbsolutePath(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount),
                "-O", outputDir.getAbsolutePath()
        };
        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir);
        checkIntervalSizes(scatterCount, outputDir, 2000000);

    }

    @Test
    public void testNoIntervals() {
        final int scatterCount = 5;
        final File outputDir = createTempDir("output");
        final String[] args = {
                "-R", REFERENCE.getAbsolutePath(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount),
                "-O", outputDir.getAbsolutePath()
        };
        runCommandLine(args);
        verifyScatteredFilesExist(scatterCount, outputDir);
        final int totalLengthInRef = GLP.getSequenceDictionary().getSequences().stream().mapToInt(SAMSequenceRecord::getSequenceLength).sum();
        checkIntervalSizes(scatterCount, outputDir, totalLengthInRef);

    }

    private static Stream<File> getScatteredFiles(final int scatterCount, final File outputDir) {
        return IntStream.range(0, scatterCount).mapToObj(n -> new File(outputDir, "000" + n + "-scattered.intervals"));
    }

    private static void verifyScatteredFilesExist(final int scatterCount, final File outputDir) {
        getScatteredFiles(scatterCount, outputDir).forEach(f -> Assert.assertTrue(f.exists()));
        Assert.assertFalse(new File(outputDir, "000" + scatterCount + "-scattered.intervals").exists());
    }

    private static List<SimpleInterval> readIntervals(final File intervalsFile) {
        return IntervalUtils.intervalFileToList(GLP, intervalsFile.getAbsolutePath()).stream().map(SimpleInterval::new).collect(Collectors.toList());
    }

    private static void checkIntervalSizes(final int scatterCount, final File outputDir, final int expectedTotalLength) {
        final int splitLength = expectedTotalLength / scatterCount;
        getScatteredFiles(scatterCount, outputDir).forEach(f -> Assert.assertEquals(readIntervals(f).stream().mapToInt(SimpleInterval::size).sum(), splitLength, 100));
    }

}