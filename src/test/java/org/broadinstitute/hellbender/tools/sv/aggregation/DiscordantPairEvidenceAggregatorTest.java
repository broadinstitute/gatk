package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class DiscordantPairEvidenceAggregatorTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICTIONARY = SVTestUtils.hg38Dict;
    private static final String TEST_EVIDENCE = toolsTestDir + "/walkers/sv/printevidence/test_hg38.pe.txt.gz";

    @DataProvider(name = "getDiscordantPairIntervalTestData")
    public Object[][] getDiscordantPairIntervalTestData() {
        return new Object[][]{
                {1100, true, 200, 500, 600, 1300},   // positive strand
                {1100, false, 200, 500, 900, 1600},  // negative strand
        };
    }

    @Test(dataProvider= "getDiscordantPairIntervalTestData")
    public void getDiscordantPairIntervalTest(final int pos, final boolean strand,
                                              final int innerWindow, final int outerWindow,
                                              final int expectedStart, final int expectedEnd) {
        final FeatureDataSource<DiscordantPairEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final DiscordantPairEvidenceAggregator aggregator = new DiscordantPairEvidenceAggregator(source, DICTIONARY,
                innerWindow, outerWindow);
        final String contig = "chr1";

        final SVCallRecord startRecord = SVTestUtils.newBndCallRecordWithPositionAndStrands(contig, pos, strand, contig, pos + 100, true);
        final SimpleInterval startWindow = aggregator.getDiscordantPairStartInterval(startRecord);
        Assert.assertEquals(startWindow.getContig(), contig);
        Assert.assertEquals(startWindow.getStart(), expectedStart);
        Assert.assertEquals(startWindow.getEnd(), expectedEnd);

        final SVCallRecord endRecord = SVTestUtils.newBndCallRecordWithPositionAndStrands(contig, pos - 100, true, contig, pos, strand);
        final SimpleInterval endWindow = aggregator.getDiscordantPairEndInterval(endRecord);
        Assert.assertEquals(endWindow.getContig(), contig);
        Assert.assertEquals(endWindow.getStart(), expectedStart);
        Assert.assertEquals(endWindow.getEnd(), expectedEnd);
    }

    @Test
    public void testGetters() {
        final FeatureDataSource<DiscordantPairEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final int innerWindow = 100;
        final int outerWindow = 200;
        final DiscordantPairEvidenceAggregator aggregator = new DiscordantPairEvidenceAggregator(source, DICTIONARY, innerWindow, outerWindow);
        Assert.assertEquals(aggregator.getInnerWindow(), innerWindow);
        Assert.assertEquals(aggregator.getOuterWindow(), outerWindow);
    }

    @DataProvider(name = "testCollectEvidenceData")
    public Object[][] testCollectEvidenceData() {
        return new Object[][] {
                // here in 1-based coordinates (raw data is 0-based)

                // Single pair, +- strands
                {25047743, true, 25048449, false, 0, 0, 1},
                {25047743, true, 25048449, true, 0, 0, 0},  // wrong strands
                {25047743, false, 25048449, true, 0, 0, 0},  // wrong strands
                {25047743, false, 25048449, false, 0, 0, 0},  // wrong strands
                {25047744, true, 25048448, false, 0, 1, 1},  // outer window +1
                {25047742, true, 25048450, false, 0, 1, 0},
                {25047742, true, 25048450, false, 1, 0, 1},  // inner window +1
                {25047744, true, 25048448, false, 1, 0, 0},

                // Single pair, -- strands
                {25052410, false, 25052645, false, 0, 0, 1},
                {25052410, true, 25052645, false, 0, 0, 0},  // wrong strands

                // 2 pairs
                {25125000, true, 31228000, false, 500, 500, 2},
        };
    }

    @Test(dataProvider= "testCollectEvidenceData")
    public void testCollectEvidence(final int posA, final boolean strandA, final int posB, final boolean strandB,
                                    final int innerWindow, final int outerWindow, final int expectedCount) {
        final FeatureDataSource<DiscordantPairEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);

        final SVCallRecord record = SVTestUtils.newBndCallRecordWithPositionAndStrands("chr21", posA, strandA,
                "chr21", posB, strandB);

        final Collection<SimpleInterval> cacheIntervals = new ArrayList<>();
        cacheIntervals.add(new SimpleInterval("chr21", 25004918, 25005918));
        if (strandA) {
            cacheIntervals.add(new SimpleInterval("chr21", posA - outerWindow, posA + innerWindow));
        } else {
            cacheIntervals.add(new SimpleInterval("chr21", posA - innerWindow, posA + outerWindow));
        }

        // No caching
        final List<DiscordantPairEvidence> testStart = new DiscordantPairEvidenceAggregator(source, DICTIONARY, innerWindow, outerWindow)
                .collectEvidence(record);
        Assert.assertEquals(testStart.size(), expectedCount);

        // With caching
        final DiscordantPairEvidenceAggregator cachedAggregatorStart = new DiscordantPairEvidenceAggregator(source, DICTIONARY, innerWindow, outerWindow);
        cachedAggregatorStart.setCacheIntervals(cacheIntervals);
        final List<DiscordantPairEvidence> testStartCached = cachedAggregatorStart.collectEvidence(record);
        Assert.assertEquals(testStartCached.size(), expectedCount);
    }

    @DataProvider(name = "testGetEvidenceQueryIntervalData")
    public Object[][] testGetEvidenceQueryIntervalData() {
        return new Object[][] {
                {1000, true, 0, 0, 1000, 1000}, // 0 window, positive strand
                {1000, false, 0, 0, 1000, 1000}, // 0 window, negative strand
                {1000, true, 1, 0, 1000, 1001}, // 1 inner window, positive strand
                {1000, true, 0, 1, 999, 1000}, // 1 outer window, positive strand
                {1000, false, 1, 0, 999, 1000}, // 1 inner window, negative strand
                {1000, false, 0, 1, 1000, 1001}, // 1 outer window, negative strand
                {1000, true, 100, 1000, 1, 1100}, // left clip window
                {1000, true, 1000000000, 1000000000, 1, 46709983}, // left and right clip window
        };
    }

    @Test(dataProvider= "testGetEvidenceQueryIntervalData")
    public void testGetEvidenceQueryInterval(final int pos, final boolean strand, final int innerWindow,
                                             final int outerWindow, final int expectedStart, final int expectedEnd) {
        final FeatureDataSource<DiscordantPairEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final SVCallRecord record = SVTestUtils.newBndCallRecordWithPositionAndStrands("chr21", pos, strand,
                "chr21", pos + 1000, true);
        final SimpleInterval expected = new SimpleInterval("chr21", expectedStart, expectedEnd);

        // Start position aggregation
        final DiscordantPairEvidenceAggregator aggregator = new DiscordantPairEvidenceAggregator(source, DICTIONARY, innerWindow, outerWindow);
        final SimpleInterval resultStart = aggregator.getEvidenceQueryInterval(record);
        Assert.assertEquals(resultStart, expected);
    }

    @Test
    public void testEvidenceFilter() {
        final String contig = "chr21";
        final int pos = 10000;
        final FeatureDataSource<DiscordantPairEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final SVCallRecord recordTrue = SVTestUtils.newBndCallRecordWithPositionAndStrands(contig, pos, true,
                contig, pos + 1000, true);
        final SVCallRecord recordFalse = SVTestUtils.newBndCallRecordWithPositionAndStrands(contig, pos, false,
                contig, pos + 1000, true);

        final DiscordantPairEvidence evidenceTrue = new DiscordantPairEvidence("", contig, pos,  true, contig, pos + 1000, true);
        final DiscordantPairEvidence evidenceFalse = new DiscordantPairEvidence("", contig, pos, false, contig, pos + 1000, true);
        final DiscordantPairEvidence evidenceTrueNonOverlapping = new DiscordantPairEvidence("", contig, pos + 10000, true, contig, pos + 1001, true);

        // Start position aggregation
        final DiscordantPairEvidenceAggregator aggregator = new DiscordantPairEvidenceAggregator(source, DICTIONARY, 0, 0);
        Assert.assertTrue(aggregator.evidenceFilter(recordTrue, evidenceTrue));
        Assert.assertFalse(aggregator.evidenceFilter(recordTrue, evidenceFalse));
        Assert.assertFalse(aggregator.evidenceFilter(recordTrue, evidenceTrueNonOverlapping));
        Assert.assertTrue(aggregator.evidenceFilter(recordFalse, evidenceFalse));
        Assert.assertFalse(aggregator.evidenceFilter(recordFalse, evidenceTrue));
    }

}