package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import com.google.common.collect.Iterables;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
import org.broadinstitute.hellbender.utils.read.markduplicates.PairedEnds;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class MarkDuplicatesSparkUtilsUnitTest extends BaseTest {

    @DataProvider(name = "fragments")
    public Object[][] loadFragments() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();

        GATKRead lowScoreRead = ArtificialReadUtils.createArtificialRead(header, "lowScoreRead", 1, 200, 10);
        lowScoreRead.setBaseQualities(Utils.dupBytes((byte) 15, 10));

        GATKRead lowScoreReadDup = ArtificialReadUtils.createArtificialRead(header, "lowScoreRead", 1, 200, 10);
        lowScoreReadDup.setBaseQualities(Utils.dupBytes((byte) 15, 10));
        lowScoreReadDup.setIsDuplicate(true);

        GATKRead highScoreRead = ArtificialReadUtils.createArtificialRead(header, "highScoreRead", 1, 200, 10);
        highScoreRead.setBaseQualities(Utils.dupBytes((byte) 30, 10));
        highScoreRead.setFragmentLength(0);

        GATKRead highScoreReadDup = ArtificialReadUtils.createArtificialRead(header, "highScoreRead", 1, 200, 10);
        highScoreReadDup.setBaseQualities(Utils.dupBytes((byte) 30, 10));
        highScoreReadDup.setFragmentLength(0);
        highScoreReadDup.setIsDuplicate(true);

        GATKRead highScoreReadTieWinner = ArtificialReadUtils.createArtificialRead(header, "highScoreRead", 1, 200, 10);
        highScoreRead.setBaseQualities(Utils.dupBytes((byte) 30, 10));
        highScoreReadTieWinner.setFragmentLength(1);

        GATKRead pairedRead = ArtificialReadUtils.createArtificialRead(header, "pairedRead", 1, 200, 10);
        pairedRead.setBaseQualities(Utils.dupBytes((byte) 18, 10));
        pairedRead.setMatePosition("1", 200);

        // Iterable<PairedEnds> pairedEnds, final SAMFileHeader header, List<GATKRead> expected, String message.
        return new Object[][]{
                {Collections.emptyList(), header, Collections.emptyList(),
                        "There should be no reads emitted if none are passed in."},
                {Collections.singletonList(PairedEnds.of(lowScoreRead)), header, Collections.singletonList(lowScoreRead),
                        "A single read should not be parked as a duplicate."},
                {Arrays.asList(PairedEnds.of(lowScoreRead), PairedEnds.of(highScoreRead)), header, Arrays.asList(highScoreRead, lowScoreReadDup),
                        "The low quality read should be marked as a duplicate."},
                {Arrays.asList(PairedEnds.of(lowScoreRead), PairedEnds.of(highScoreRead), PairedEnds.of(highScoreReadTieWinner)), header,
                        Arrays.asList(highScoreReadTieWinner, highScoreReadDup, lowScoreReadDup),
                        "The scores of highScoreRead/2 are the same, the fragment length should be used as the tie-breaker."},
                {Arrays.asList(PairedEnds.of(lowScoreRead), PairedEnds.of(highScoreRead), PairedEnds.of(pairedRead)), header,
                        Arrays.asList(lowScoreReadDup, highScoreReadDup),
                        "If a read with a mapped mate is present, all fragments should be marked as duplicates (the paired read is not emitted)."},
        };
    }

    @DataProvider(name = "pairedEnds")
    public Object[][] loadPairedEnds() {
        JavaSparkContext context = SparkContextFactory.getTestSparkContext();
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();

        List<Iterable<PairedEnds>> noPairs = Collections.emptyList();
        List<GATKRead> noReads = Collections.emptyList();
        JavaRDD<Iterable<PairedEnds>> noPairsRDD = context.parallelize(noPairs);
        JavaPairRDD<String, Iterable<PairedEnds>> noPairsKeyed = noPairsRDD.keyBy(v1 -> Iterables.getFirst(v1, null).key(header));

        GATKRead lowScoreRead = ArtificialReadUtils.createArtificialRead(header, "lowScoreRead", 1, 200, 10);
        List<Iterable<PairedEnds>> singlePair = Collections.singletonList(Collections.singletonList(PairedEnds.of(lowScoreRead)));

        JavaPairRDD<String, Iterable<PairedEnds>> singleKeyed = context.parallelize(singlePair).keyBy(v1 -> Iterables.getFirst(v1, null).key(header));


        // final JavaPairRDD<String, Iterable<PairedEnds>> keyedPairs, final OpticalDuplicateFinder finder,
        // final SAMFileHeader header, JavaRDD<GATKRead> expected, String message.
        return new Object[][]{
                {noPairsKeyed, null, header, noReads,
                        "There should be no reads emitted if none are passed in."},
                {singleKeyed, null, header, Collections.singletonList(lowScoreRead),
                        "No-op for fragments."},
        };
    }

    @Test(dataProvider = "fragments")
    void handleFragmentsTest(Iterable<PairedEnds> pairedEnds, final SAMFileHeader header, List<GATKRead> expected, String message) {
        List<GATKRead> reads = MarkDuplicatesSparkUtils.handleFragments(pairedEnds, header);
        Assert.assertEquals(reads, expected, message);
    }

    @Test(dataProvider = "pairedEnds")
    void markPairedEndsTest(final JavaPairRDD<String, Iterable<PairedEnds>> keyedPairs,
                            final OpticalDuplicateFinder finder, final SAMFileHeader header, List<GATKRead> expected, String message) {
        List<GATKRead> reads = MarkDuplicatesSparkUtils.markPairedEnds(keyedPairs, finder, header).collect();
        Assert.assertEquals(reads, expected, message);
    }
}