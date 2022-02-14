package org.broadinstitute.hellbender.tools.gvs.common;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.SplitIntervals;
import org.broadinstitute.hellbender.tools.walkers.SplitIntervalsIntegrationTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;

import static java.util.Arrays.asList;


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
        final Interval interval = new Interval("chr20", 1000000, 2000000);

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


