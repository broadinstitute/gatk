package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class BafEvidenceAggregatorTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICTIONARY = SVTestUtils.hg38Dict;
    private static final String TEST_EVIDENCE = toolsTestDir + "/walkers/sv/printevidence/test_hg38.baf.txt.gz";

    @DataProvider(name = "testCollectBafEvidenceData")
    public Object[][] testCollectBafEvidenceData() {
        return new Object[][] {
                {"chr21", 25794845, 1, 0, 0},
                {"chr21", 25794845, 1, 1, 0},
                {"chr21", 25794388, 1000, 0, 11},
                {"chr21", 25794389, 1000, 0, 11},
                {"chr21", 25794390, 1000, 0, 10},
                {"chr21", 25794845, 1000, 0.5, 18},
                {"chr21", 25984167, 1000, 0.5, 6},
                {"chr22", 30468112, 1000, 0.5, 14},
        };
    }

    @Test(dataProvider= "testCollectBafEvidenceData")
    public void testCollectBafEvidence(final String contig, final int pos, final int length, final double padding, final int expectedCount) {
        final FeatureDataSource<BafEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final SVCallRecord record = SVTestUtils.newDepthCallRecordWithIntervalAndType(contig, pos, pos + length - 1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);

        final List<BafEvidence> test = new BafEvidenceAggregator(source, DICTIONARY, padding)
                .collectEvidence(record);
        Assert.assertEquals(test.size(), expectedCount);
    }

    @Test
    public void testGetters() {
        final FeatureDataSource<BafEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final double padding = 0.1;
        final BafEvidenceAggregator aggregator = new BafEvidenceAggregator(source, DICTIONARY, padding);
        Assert.assertEquals(aggregator.getPaddingFraction(), padding);
    }

    @DataProvider(name = "testGetEvidenceQueryIntervalData")
    public Object[][] testGetEvidenceQueryIntervalData() {
        return new Object[][] {
                {"chr21", 1000, 1, 0, 1000, 1000},
                {"chr21", 1000, 1, 0.5, 999, 1001},
                {"chr21", 1000, 1, 1.1, 998, 1002},
                {"chr21", 1000, 1000, 0, 1000, 1999},
                {"chr21", 1000, 500, 0, 1000, 1499},
                {"chr21", 1000, 500, 0.5, 750, 1749},
                {"chr21", 1000, 501, 0.5, 749, 1751},
                {"chr21", 1000, 46708983, 0.5, 1, 46709983},  // fit to whole contig
        };
    }

    @Test(dataProvider= "testGetEvidenceQueryIntervalData")
    public void testGetEvidenceQueryInterval(final String contig, final int pos, final int length, final double padding, final int expectedStart, final int expectedEnd) {
        final FeatureDataSource<BafEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final SVCallRecord record = SVTestUtils.newDepthCallRecordWithIntervalAndType(contig, pos, pos + length - 1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final SimpleInterval expected = new SimpleInterval(contig, expectedStart, expectedEnd);
        final BafEvidenceAggregator startAggregator = new BafEvidenceAggregator(source, DICTIONARY, padding);
        final SimpleInterval resultStart = startAggregator.getEvidenceQueryInterval(record);
        Assert.assertEquals(resultStart, expected);
    }

    @Test
    public void testEvidenceFilter() {
        final String contig = "chr21";
        final int pos = 10000;
        final int length = 1000;
        final FeatureDataSource<BafEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final SVCallRecord record = SVTestUtils.newDepthCallRecordWithIntervalAndType(contig, pos, pos + length - 1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final BafEvidence evidence = new BafEvidence("", contig, pos, 1);
        final BafEvidenceAggregator aggregator = new BafEvidenceAggregator(source, DICTIONARY, 0);
        Assert.assertTrue(aggregator.evidenceFilter(record, evidence));
    }

}