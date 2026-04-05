package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.mockito.internal.util.collections.Sets;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class BafHetRatioTesterTest {

    private static final String SAMPLE_1 = "sample1";
    private static final String SAMPLE_2 = "sample2";
    private static final String SAMPLE_3 = "sample3";
    private static final String TEST_CONTIG = "chr21";

    @Test
    public void test() {
        final int seed = 0;
        final double tolerance = 1e-4;
        // Use a DEL with enough het SNPs in flanks and inside to avoid ROH classification
        final SVCallRecord record = SVTestUtils.newDepthCallRecordWithIntervalAndType("chr21", 10000, 20000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);

        // Create evidence with enough het SNPs so samples don't get classified as ROH:
        // epsilon = min(50/10000, 0.0005) = 0.0005; threshold = 0.0005
        // Need F/L >= 0.0005 or M/L >= 0.0005 for at least one to be true
        // L = 10000, so need at least 5 hets in flanks or inside
        final List<BafEvidence> evidence = Lists.newArrayList();

        // sample1 (control): 6 before, 6 inside, 6 after
        for (int i = 0; i < 6; i++) {
            evidence.add(new BafEvidence(SAMPLE_1, TEST_CONTIG, 9000 + i * 100, 0.5));
            evidence.add(new BafEvidence(SAMPLE_1, TEST_CONTIG, 12000 + i * 1000, 0.5));
            evidence.add(new BafEvidence(SAMPLE_1, TEST_CONTIG, 21000 + i * 100, 0.5));
        }

        // sample2 (control): 6 before, 6 inside, 6 after
        for (int i = 0; i < 6; i++) {
            evidence.add(new BafEvidence(SAMPLE_2, TEST_CONTIG, 9000 + i * 100, 0.5));
            evidence.add(new BafEvidence(SAMPLE_2, TEST_CONTIG, 12000 + i * 1000, 0.5));
            evidence.add(new BafEvidence(SAMPLE_2, TEST_CONTIG, 21000 + i * 100, 0.5));
        }

        // sample3 (carrier): 6 before, 2 inside (deletion!), 6 after
        for (int i = 0; i < 6; i++) {
            evidence.add(new BafEvidence(SAMPLE_3, TEST_CONTIG, 9000 + i * 100, 0.5));
            evidence.add(new BafEvidence(SAMPLE_3, TEST_CONTIG, 21000 + i * 100, 0.5));
        }
        evidence.add(new BafEvidence(SAMPLE_3, TEST_CONTIG, 12000, 0.5));
        evidence.add(new BafEvidence(SAMPLE_3, TEST_CONTIG, 15000, 0.5));

        final Set<String> allSamples = Sets.newSet(SAMPLE_1, SAMPLE_2, SAMPLE_3);
        final Set<String> carrierSamples = Sets.newSet(SAMPLE_3);
        final int svLength = 10000;

        final BafHetRatioTester tester = new BafHetRatioTester(seed);
        final BafHetRatioTester.BafDelResult result = tester.test(record, evidence, allSamples, carrierSamples, svLength);
        Assert.assertNotNull(result);

        // Verify the carrier (sample3) has elevated het ratio (10^(-Ratio)):
        // sample3: F=6, M=2, E=6, epsilon=0.0005, L=10000
        // flank = min(6, 6) = 6, Ratio = log10((2 + 5)/(6 + 5)) = log10(7/11) = -0.1961
        // hetRatio = 10^(0.1961) = 1.571
        SVTestUtils.assertFloatWithinTolerance(result.hetRatio(), 1.5714, tolerance);

        // Null test cases
        final SVCallRecord recordNotCnv = SVTestUtils.newCallRecordInsertionWithLengthAndCoordinates("chr21", 10000, 10000);
        Assert.assertNull(tester.test(recordNotCnv, evidence, allSamples, carrierSamples, svLength));
        Assert.assertNull(tester.test(record, null, allSamples, carrierSamples, svLength));
        Assert.assertNull(tester.test(record, Collections.emptyList(), allSamples, carrierSamples, svLength));
    }

    @Test
    public void testApplyToRecord() {
        final BafHetRatioTester tester = new BafHetRatioTester(0);
        final SVCallRecord record = SVTestUtils.newDepthCallRecordWithIntervalAndType("chr21", 10000, 20000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final SVCallRecord resultNull = tester.applyToRecord(record, null);
        SVTestUtils.assertEquals(resultNull, record);

        final BafHetRatioTester.BafDelResult delResult = new BafHetRatioTester.BafDelResult(1.5, 0.8);
        final SVCallRecord resultRecord = tester.applyToRecord(record, delResult);
        Assert.assertEquals(resultRecord.getAttributes().get(GATKSVVCFConstants.BAF_HET_RATIO_ATTRIBUTE), 1.5);
        Assert.assertEquals(resultRecord.getAttributes().get(GATKSVVCFConstants.BAF_DEL_LOGLIK_ATTRIBUTE), 0.8);

        final BafHetRatioTester.BafDelResult nanLoglik = new BafHetRatioTester.BafDelResult(1.5, Double.NaN);
        final SVCallRecord resultWithoutLoglik = tester.applyToRecord(record, nanLoglik);
        Assert.assertEquals(resultWithoutLoglik.getAttributes().get(GATKSVVCFConstants.BAF_HET_RATIO_ATTRIBUTE), 1.5);
        Assert.assertFalse(resultWithoutLoglik.getAttributes().containsKey(GATKSVVCFConstants.BAF_DEL_LOGLIK_ATTRIBUTE));
    }

    @Test
    public void testCalculateHandlesUntestableGmmCases() {
        final BafHetRatioTester tester = new BafHetRatioTester(0);

        final Map<String, int[]> manyCarriersFewControls = new LinkedHashMap<>();
        manyCarriersFewControls.put("control", new int[]{6, 6, 6});
        manyCarriersFewControls.put("carrier1", new int[]{6, 2, 6});
        manyCarriersFewControls.put("carrier2", new int[]{6, 2, 6});
        final BafHetRatioTester.BafDelResult tooManyCarriers = tester.calculate(manyCarriersFewControls, Sets.newSet("carrier1", "carrier2"), 10000);
        Assert.assertNotNull(tooManyCarriers);
        SVTestUtils.assertFloatWithinTolerance(tooManyCarriers.hetRatio(), 11.0 / 7.0, 1e-6);
        Assert.assertTrue(Double.isNaN(tooManyCarriers.delLoglik()));

        final List<Double> tenControls = Arrays.asList(-0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15);
        final List<Double> oneCarrier = Collections.singletonList(-0.40);
        Assert.assertTrue(Double.isNaN(tester.computeGmmLoglik(tenControls, oneCarrier)));

        final List<Double> elevenFlatControls = Arrays.asList(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1);
        Assert.assertTrue(Double.isNaN(tester.computeGmmLoglik(elevenFlatControls, oneCarrier)));

        final List<Double> elevenVariableControls = Arrays.asList(-0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20);
        final double loglik = tester.computeGmmLoglik(elevenVariableControls, oneCarrier);
        Assert.assertFalse(Double.isNaN(loglik));
    }
}