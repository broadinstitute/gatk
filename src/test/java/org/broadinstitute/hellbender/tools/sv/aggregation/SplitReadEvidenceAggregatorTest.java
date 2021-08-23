package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class SplitReadEvidenceAggregatorTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICTIONARY = SVTestUtils.hg38Dict;
    private static final String TEST_EVIDENCE = toolsTestDir + "/walkers/sv/printevidence/test_hg38.sr.txt.gz";

    @DataProvider(name = "testCollectEvidenceData")
    public Object[][] testCollectEvidenceData() {
        return new Object[][] {
                // 1-based coordinates, positive (true) strand <=> "left", negative (false) strand <=> "right"
                {25006096, true, 0, 3},
                {25006096, false, 0, 0},  // wrong strand
                {25006095, true, 0, 1},   // position off by 1
                {25006097, true, 0, 0},   // position off by 1
                {25006096, true, 1, 4},   // window +1
                {25006096, true, 2, 4},   // window +2
                {25006096, true, 3, 5},   // window +3
                {25005791, false, 0, 1},  // negative strand evidence
                {25006981, true, 40, 2},  // filter out one negative strand record within window
        };
    }

    @Test(dataProvider= "testCollectEvidenceData")
    public void testCollectEvidence(final int pos, final boolean strand, final int window, final int expectedCount) {
        final FeatureDataSource<SplitReadEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);

        final SVCallRecord startRecord = SVTestUtils.newBndCallRecordWithPositionAndStrands("chr21", pos, strand,
                "chr21", pos + 1000, true);
        final SVCallRecord endRecord = SVTestUtils.newBndCallRecordWithPositionAndStrands("chr21", pos - 1000, true,
                "chr21", pos, strand);
        final Collection<SimpleInterval> cacheIntervals = new ArrayList<>();
        cacheIntervals.add(new SimpleInterval("chr21", 25247871, 25247871));
        cacheIntervals.add(new SimpleInterval("chr21", 25249566, 25249566));
        cacheIntervals.add(new SimpleInterval("chr21", 25249466, 25249666));
        cacheIntervals.add(new SimpleInterval("chr22", 25247871, 25247871));
        cacheIntervals.add(new SimpleInterval("chr21", pos, pos));
        cacheIntervals.add(new SimpleInterval("chr21", pos - window, pos + window));

        // Start SR, no caching
        final List<SplitReadEvidence> testStart = new SplitReadEvidenceAggregator(source, DICTIONARY, window, true)
                .collectEvidence(startRecord);
        Assert.assertEquals(testStart.size(), expectedCount);

        // End SR, no caching
        final List<SplitReadEvidence> testEnd = new SplitReadEvidenceAggregator(source, DICTIONARY, window, false)
                .collectEvidence(endRecord);
        Assert.assertEquals(testEnd.size(), expectedCount);

        // Start SR, with caching
        final SplitReadEvidenceAggregator cachedAggregatorStart = new SplitReadEvidenceAggregator(source, DICTIONARY, window, true);
        cachedAggregatorStart.setCacheIntervals(cacheIntervals);
        final List<SplitReadEvidence> testStartCached = cachedAggregatorStart.collectEvidence(startRecord);
        Assert.assertEquals(testStartCached.size(), expectedCount);

        // End SR, with caching
        final SplitReadEvidenceAggregator cachedAggregatorEnd = new SplitReadEvidenceAggregator(source, DICTIONARY, window, false);
        cachedAggregatorEnd.setCacheIntervals(cacheIntervals);
        final List<SplitReadEvidence> testEndCached = cachedAggregatorEnd.collectEvidence(endRecord);
        Assert.assertEquals(testEndCached.size(), expectedCount);
    }

    @Test
    public void testGetters() {
        final FeatureDataSource<SplitReadEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final int window = 100;
        final SplitReadEvidenceAggregator aggregator = new SplitReadEvidenceAggregator(source, DICTIONARY, window, true);
        Assert.assertEquals(aggregator.getWindow(), window);
    }

    @DataProvider(name = "testGetEvidenceQueryIntervalData")
    public Object[][] testGetEvidenceQueryIntervalData() {
        return new Object[][] {
                {1000, 0, 1000, 1000}, // 0 window
                {1000, 1, 999, 1001}, // 1 window
                {1000, 1000, 1, 2000}, // left clip window
                {1000, 1000000000, 1, 46709983}, // left and right clip window
        };
    }

    @Test(dataProvider= "testGetEvidenceQueryIntervalData")
    public void testGetEvidenceQueryInterval(final int pos, final int window, final int expectedStart, final int expectedEnd) {
        final FeatureDataSource<SplitReadEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final SVCallRecord startRecord = SVTestUtils.newBndCallRecordWithPositionAndStrands("chr21", pos, true,
                "chr21", pos + 1000, true);
        final SVCallRecord endRecord = SVTestUtils.newBndCallRecordWithPositionAndStrands("chr21", pos - 1000, true,
                "chr21", pos, true);
        final SimpleInterval expected = new SimpleInterval("chr21", expectedStart, expectedEnd);

        // Start position aggregation
        final SplitReadEvidenceAggregator startAggregator = new SplitReadEvidenceAggregator(source, DICTIONARY, window, true);
        final SimpleInterval resultStart = startAggregator.getEvidenceQueryInterval(startRecord);
        Assert.assertEquals(resultStart, expected);

        // End position aggregation
        final SplitReadEvidenceAggregator endAggregator = new SplitReadEvidenceAggregator(source, DICTIONARY, window, false);
        final SimpleInterval resultEnd = endAggregator.getEvidenceQueryInterval(endRecord);
        Assert.assertEquals(resultEnd, expected);
    }

    @Test
    public void testEvidenceFilter() {
        final String contig = "chr21";
        final int pos = 10000;
        final FeatureDataSource<SplitReadEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final SVCallRecord startRecordTrue = SVTestUtils.newBndCallRecordWithPositionAndStrands(contig, pos, true,
                contig, pos + 1000, true);
        final SVCallRecord startRecordFalse = SVTestUtils.newBndCallRecordWithPositionAndStrands(contig, pos, false,
                contig, pos + 1000, true);
        final SVCallRecord endRecordTrue = SVTestUtils.newBndCallRecordWithPositionAndStrands(contig, pos - 1000, true,
                contig, 1000, true);
        final SVCallRecord endRecordFalse = SVTestUtils.newBndCallRecordWithPositionAndStrands(contig, pos - 1000, true,
                contig, 1000, false);

        final SplitReadEvidence evidenceTrue = new SplitReadEvidence("", contig, pos, 1, true);
        final SplitReadEvidence evidenceFalse = new SplitReadEvidence("", contig, pos, 1, false);

        // Evidence filter should not check for overlap
        final SplitReadEvidence evidenceTrueNonOverlapping = new SplitReadEvidence("", contig, pos + 10000, 1, true);
        final SplitReadEvidence evidenceFalseNonOverlapping = new SplitReadEvidence("", contig, pos + 10000, 1, false);

        // Start position aggregation
        final SplitReadEvidenceAggregator startAggregator = new SplitReadEvidenceAggregator(source, DICTIONARY, 0, true);
        Assert.assertTrue(startAggregator.evidenceFilter(startRecordTrue, evidenceTrue));
        Assert.assertFalse(startAggregator.evidenceFilter(startRecordTrue, evidenceFalse));
        Assert.assertTrue(startAggregator.evidenceFilter(startRecordTrue, evidenceTrueNonOverlapping));
        Assert.assertFalse(startAggregator.evidenceFilter(startRecordTrue, evidenceFalseNonOverlapping));
        Assert.assertTrue(startAggregator.evidenceFilter(startRecordFalse, evidenceFalse));
        Assert.assertFalse(startAggregator.evidenceFilter(startRecordFalse, evidenceTrue));

        // End position aggregation
        final SplitReadEvidenceAggregator endAggregator = new SplitReadEvidenceAggregator(source, DICTIONARY, 0, false);
        Assert.assertTrue(endAggregator.evidenceFilter(endRecordTrue, evidenceTrue));
        Assert.assertFalse(endAggregator.evidenceFilter(endRecordTrue, evidenceFalse));
        Assert.assertTrue(startAggregator.evidenceFilter(endRecordTrue, evidenceTrueNonOverlapping));
        Assert.assertFalse(startAggregator.evidenceFilter(endRecordTrue, evidenceFalseNonOverlapping));
        Assert.assertTrue(endAggregator.evidenceFilter(endRecordFalse, evidenceFalse));
        Assert.assertFalse(endAggregator.evidenceFilter(endRecordFalse, evidenceTrue));
    }
}