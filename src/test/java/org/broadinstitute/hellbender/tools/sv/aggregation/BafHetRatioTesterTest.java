package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.mockito.internal.util.collections.Sets;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;
import java.util.Set;

public class BafHetRatioTesterTest {

    private static final String SAMPLE_1 = "sample1";
    private static final String SAMPLE_2 = "sample2";
    private static final String SAMPLE_3 = "sample3";
    private static final String TEST_CONTIG = "chr21";

    @Test
    public void test() {
        final double pSnp = 0.002;
        final double pMaxHomozygous = 0.06;
        final double tolerance = 1e-4;
        final SVCallRecord record = SVTestUtils.newDepthCallRecordWithIntervalAndType("chr21", 10000, 20000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);

        final List<BafEvidence> evidence = Lists.newArrayList(
            new BafEvidence(SAMPLE_1, TEST_CONTIG, 9500, 0.6),
            new BafEvidence(SAMPLE_2, TEST_CONTIG, 9500, 0.4),
            new BafEvidence(SAMPLE_3, TEST_CONTIG, 9500, 0.9),
            new BafEvidence(SAMPLE_3, TEST_CONTIG, 15000, 0.9),
            new BafEvidence(SAMPLE_3, TEST_CONTIG, 16000, 0.1),
            new BafEvidence(SAMPLE_1, TEST_CONTIG, 16000, 0.5),
            new BafEvidence(SAMPLE_1, TEST_CONTIG, 20500, 0.5),
            new BafEvidence(SAMPLE_3, TEST_CONTIG, 20500, 0.4)
        );

        final Set<String> allSamples = Sets.newSet(SAMPLE_1, SAMPLE_2, SAMPLE_3);
        final Set<String> carrierSamples = Sets.newSet(SAMPLE_3);
        final int flankSize = 1000;

        final BafHetRatioTester tester = new BafHetRatioTester(pSnp, pMaxHomozygous);
        final Double result = tester.test(record, evidence, allSamples, carrierSamples, flankSize);
        SVTestUtils.assertFloatWithinTolerance(result, -0.7520386983881369, tolerance);

        // Null test cases
        final SVCallRecord recordNotCnv = SVTestUtils.newCallRecordInsertionWithLengthAndCoordinates("chr21", 10000, 10000);
        Assert.assertNull(tester.test(recordNotCnv, evidence, allSamples, carrierSamples, flankSize));
        Assert.assertNull(tester.test(record, null, allSamples, carrierSamples, flankSize));
        Assert.assertNull(tester.test(record, Collections.emptyList(), allSamples, carrierSamples, flankSize));
    }

    @Test
    public void testApplyToRecord() {
        final BafHetRatioTester tester = new BafHetRatioTester(0.1, 0.2);
        final SVCallRecord record = SVTestUtils.newDepthCallRecordWithIntervalAndType("chr21", 10000, 20000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final SVCallRecord resultNulll = tester.applyToRecord(record, null);
        SVTestUtils.assertEquals(resultNulll, record);

        final SVCallRecord expected = SVCallRecordUtils.copyCallWithNewAttributes(record, Collections.singletonMap(GATKSVVCFConstants.BAF_HET_RATIO_ATTRIBUTE, -0.3));
        final SVCallRecord result = tester.applyToRecord(record, -0.3);
        SVTestUtils.assertEquals(result, expected);
    }
}