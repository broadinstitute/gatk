package org.broadinstitute.hellbender.tools.gvs.common;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.SplitIntervals;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.stream.Collectors;


// Suppresses the "expression might be null" warnings which would be genuine concerns in production code but which
// are fine in tests (unexpected and dereferenced nulls should crash the test which will appropriately signal failure).
@SuppressWarnings("ConstantConditions")
public class WeightedSplitIntervalsIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testNoLossSimple() {
        final File weights = new File(publicTestDir + "example_weights_chr20_chr21.bed.gz");
        final int scatterCount = 100;
        final File outputDir = createTempDir("output");
        final Interval interval = new Interval("chr20", 1000000, 2000000);

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval(interval)
                .addReference(b38_reference_20_21)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(WeightedSplitIntervals.WEIGHTS_BED_FILE_FULL_NAME, weights)
                .addOutput(outputDir);

        runCommandLine(args);

        // verify we haven't lost any bases
        IntervalList outList = IntervalList.fromFiles(Arrays.stream(outputDir.listFiles()).collect(Collectors.toList())).uniqued().sorted();

        // total size is the same
        Assert.assertEquals(outList.getBaseCount(),1000001);

        // we should have just one interval
        Assert.assertEquals(outList.getIntervals().size(), 1);

        // and it should be the one we started with!
        Assert.assertEquals(outList.getIntervals().get(0), interval);
    }

    @Test
    public void testNoLossRealisticWgs() {
        final File weights = new File(publicTestDir + "example_weights_chr20_chr21.bed.gz");
        final File wgsIntervalList = new File(publicTestDir + "hg38.handcurated.noCentromeres.noTelomeres.chr20_chr21.interval_list");

        final int scatterCount = 100;
        final File outputDir = createTempDir("output");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval(wgsIntervalList.getAbsolutePath())
                .addReference(b38_reference_20_21)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(WeightedSplitIntervals.WEIGHTS_BED_FILE_FULL_NAME, weights)
                .addOutput(outputDir);

        runCommandLine(args);

        // verify we haven't lost any bases
        IntervalList original = IntervalList.fromFile(wgsIntervalList).uniqued().sorted();
        IntervalList outList = IntervalList.fromFiles(Arrays.stream(outputDir.listFiles()).collect(Collectors.toList())).uniqued().sorted();

        // total size is the same
        Assert.assertEquals(outList.getBaseCount(),original.getBaseCount());

        // and the number of intervals are the same
        Assert.assertEquals(outList.getIntervals().size(), original.getIntervals().size());

        // ensure they are the same
        Assert.assertEquals(0, IntervalList.difference(original, outList).size());

    }

    @Test
    public void testHandleZeroBasesToTake() {
        final File weights = new File(publicTestDir + "example_weights_chr20_chr21.bed.gz");
        final File wgsIntervalList = new File(publicTestDir + "hg38.handcurated.noCentromeres.noTelomeres.chr20_chr21.interval_list");

        // This value of scatter width exercises the `WeightedSplitIntervals` zero `basesToTake` condition of VS-384.
        // While the zero `basesToTake` condition is known to be exercised by other scatter widths (18, 40, 59, and 77),
        // this test uses 94 as that's the closest to the "natural" choice of scatter width 100 in the
        // `testNoLossRealisticWgs` test from which this test descended.
        final int scatterCount = 94;
        final File outputDir = createTempDir("output");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval(wgsIntervalList.getAbsolutePath())
                .addReference(b38_reference_20_21)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(WeightedSplitIntervals.WEIGHTS_BED_FILE_FULL_NAME, weights)
                .addOutput(outputDir);

        runCommandLine(args);

        // verify we haven't lost any bases
        IntervalList original = IntervalList.fromFile(wgsIntervalList).uniqued().sorted();
        IntervalList outList = IntervalList.fromFiles(Arrays.stream(outputDir.listFiles()).collect(Collectors.toList())).uniqued().sorted();

        // total size is the same
        Assert.assertEquals(outList.getBaseCount(),original.getBaseCount());

        // and the number of intervals are the same
        Assert.assertEquals(outList.getIntervals().size(), original.getIntervals().size());

        // ensure they are the same
        Assert.assertEquals(0, IntervalList.difference(original, outList).size());
    }

    @Test
    public void testHandleNotOneMoreBase() {
        // Sometimes, the remaining space in a scatter interval is non-zero, yet can't take
        // even a single base of the next interval.
        // e.g.
        // 1 bases, 200 weight each
        // 10 bases 1000 weight each
        // target: sum 10200 / 10 -> 1020 per shard
        // chr20	60000	60001	.	200
        // chr20	60001	60011	.	10000
        final File weights = new File(publicTestDir + "example_weights_chr20_test_zero.bed");
        final int scatterCount = 10;
        final File outputDir = createTempDir("output");
        final Interval inputInterval = new Interval("chr20",60000,60109);

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval(inputInterval.getContig() + ":" + inputInterval.getStart() + "-" + inputInterval.getEnd())
                .addReference(hg38Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(WeightedSplitIntervals.WEIGHTS_BED_FILE_FULL_NAME, weights)
                .addOutput(outputDir);
        runCommandLine(args);

        // even though we asked for 10 scatter, we should get 11 because the target per shard is 1020 weight
        // and due to the size of our weights each of the 11 bases must be in their own shard
        Assert.assertEquals(outputDir.listFiles().length, 11);

        // verify we have exactly the input intervals
        IntervalList outList = IntervalList.fromFiles(Arrays.stream(outputDir.listFiles()).collect(Collectors.toList())).uniqued().sorted();

        // assert it's a single interval
        Assert.assertEquals(outList.getIntervals().size(), 1);

        // and the interval itself is the same
        Assert.assertEquals(outList.getIntervals().get(0).compareTo(inputInterval), 0);

    }

    @Test
    public void testDontMixContigs() {
        final File weights = new File(publicTestDir + "example_weights_chr20_chr21.bed.gz");
        final int scatterCount = 1;
        final File outputDir = createTempDir("output");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInterval("chr20:1000000-2000000")
                .addInterval("chr21:1000000-2000000")
                .addReference(hg38Reference)
                .add(SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount)
                .add(WeightedSplitIntervals.WEIGHTS_BED_FILE_FULL_NAME, weights)
                .add(SplitIntervals.DONT_MIX_CONTIGS_LONG_NAME, true)
                .addOutput(outputDir);
        runCommandLine(args);

        // even though we asked for one scatter, we should get two because of not mixing contigs
        Assert.assertEquals(outputDir.listFiles().length, 2);
    }

}


