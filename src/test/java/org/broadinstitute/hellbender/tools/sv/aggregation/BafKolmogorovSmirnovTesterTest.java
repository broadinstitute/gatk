package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public class BafKolmogorovSmirnovTesterTest {

    private static final String SAMPLE_PREFIX = "sample";
    private static final String CARRIER_PREFIX = "car";
    private static final String TEST_CONTIG = "chr21";

    private List<BafEvidence> generateTestEvidence(final int numCarriers, final int numBackground) {
        final List<BafEvidence> evidences = Lists.newArrayList();
        for (int i = 9000; i < 21000; i += 2000) {  // spatial coordinates
            for (int j = 0; j < numCarriers; j++) {  // carrier samples
                evidences.add(new BafEvidence(CARRIER_PREFIX + j, TEST_CONTIG, i, 0.2 + (0.7 * j / (double) numCarriers)));
            }
            for (int j = 0; j < numBackground; j++) {  // background samples
                evidences.add(new BafEvidence(SAMPLE_PREFIX + j, TEST_CONTIG, i, j / (double) numBackground));
            }
        }
        return evidences;
    }

    private Set<String> generatTestCarriers(final int numCarriers) {
        final Set<String> carrierSamples = new HashSet<>();
        for (int j = 0; j < numCarriers; j++) {
            carrierSamples.add(CARRIER_PREFIX + j);
        }
        return carrierSamples;
    }

    @Test
    public void test() {
        final int minBafCount = 1;
        final long seed = 42;
        final SVCallRecord record = SVTestUtils.newDepthCallRecordWithIntervalAndType("chr21", 10000, 20000, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP);

        // Small - exact method
        final int numCarriers = 10;
        final List<BafEvidence> evidence = generateTestEvidence(numCarriers, 30);
        final Set<String> carrierSamples = generatTestCarriers(numCarriers);
        final BafKolmogorovSmirnovTester tester = new BafKolmogorovSmirnovTester(minBafCount, seed);
        final BafKolmogorovSmirnovTester.KSTestResult result = tester.test(record, evidence, carrierSamples);
        SVTestUtils.assertFloatWithinTolerance(result.getStat(), 0.2, 1e-4);
        SVTestUtils.assertFloatWithinTolerance(result.getP(), 0.09290160749765941, 1e-4);

        // Large - approximate method
        final int numCarriersLarge = 100;
        final List<BafEvidence> evidenceLarge = generateTestEvidence(numCarriersLarge, 300);
        final Set<String> carrierSamplesLarge = generatTestCarriers(numCarriersLarge);
        final BafKolmogorovSmirnovTester.KSTestResult resultLarge = tester.test(record, evidenceLarge, carrierSamplesLarge);
        SVTestUtils.assertFloatWithinTolerance(resultLarge.getStat(), 0.2, 1e-4);
        SVTestUtils.assertFloatWithinTolerance(resultLarge.getP(), 1.871836019518014E-13, 1e-16);

        // Null test cases
        final SVCallRecord recordNotCnv = SVTestUtils.newCallRecordInsertionWithLengthAndCoordinates("chr21", 10000, 10000);
        Assert.assertNull(tester.test(recordNotCnv, evidence, carrierSamples));
        Assert.assertNull(tester.test(record, null, carrierSamples));
        Assert.assertNull(tester.test(record, Collections.emptyList(), carrierSamples));
    }

    /**
     * Tests that the singleton KS fallback matches scipy.stats.ks_2samp for cases where one array has length 1.
     * Reference values from: scipy.stats.ks_2samp(carrier, null)
     */
    @Test
    public void testSingletonKs() {
        // Case 1: 1 carrier vs 5 null -> scipy: stat=0.8, pval=0.666667
        double[] result1 = BafKolmogorovSmirnovTester.singletonKsTest(
                new double[]{0.45}, new double[]{0.30, 0.35, 0.40, 0.42, 0.50});
        SVTestUtils.assertFloatWithinTolerance(result1[0], 0.8, 1e-10);
        SVTestUtils.assertFloatWithinTolerance(result1[1], 2.0 / 3.0, 1e-10);

        // Case 2: 1 carrier vs 1 null (different) -> scipy: stat=1.0, pval=1.0
        double[] result2 = BafKolmogorovSmirnovTester.singletonKsTest(
                new double[]{0.45}, new double[]{0.30});
        SVTestUtils.assertFloatWithinTolerance(result2[0], 1.0, 1e-10);
        SVTestUtils.assertFloatWithinTolerance(result2[1], 1.0, 1e-10);

        // Case 3: 5 carrier vs 1 null (symmetric) -> scipy: stat=0.8, pval=0.666667
        double[] result3 = BafKolmogorovSmirnovTester.singletonKsTest(
                new double[]{0.30, 0.35, 0.40, 0.42, 0.50}, new double[]{0.45});
        SVTestUtils.assertFloatWithinTolerance(result3[0], 0.8, 1e-10);
        SVTestUtils.assertFloatWithinTolerance(result3[1], 2.0 / 3.0, 1e-10);

        // Case 4: identical singletons -> scipy: stat=0.0, pval=1.0
        double[] result4 = BafKolmogorovSmirnovTester.singletonKsTest(
                new double[]{0.45}, new double[]{0.45});
        SVTestUtils.assertFloatWithinTolerance(result4[0], 0.0, 1e-10);
        SVTestUtils.assertFloatWithinTolerance(result4[1], 1.0, 1e-10);

        // Case 5: 1 carrier vs 2 null -> scipy: stat=0.5, pval=1.0
        double[] result5 = BafKolmogorovSmirnovTester.singletonKsTest(
                new double[]{0.45}, new double[]{0.30, 0.50});
        SVTestUtils.assertFloatWithinTolerance(result5[0], 0.5, 1e-10);
        SVTestUtils.assertFloatWithinTolerance(result5[1], 1.0, 1e-10);
    }

    /**
     * Integration test: singleton carrier evidence should produce a non-null result (no longer blocked by the n>=2 guard).
     */
    @Test
    public void testSingletonCarrierEvidence() {
        final int minBafCount = 1;
        final long seed = 42;
        final SVCallRecord record = SVTestUtils.newDepthCallRecordWithIntervalAndType(
                "chr21", 10000, 20000, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP);

        // 1 carrier sample with 1 BAF observation, 5 background samples each with 1 BAF observation
        final List<BafEvidence> evidence = new ArrayList<>();
        evidence.add(new BafEvidence(CARRIER_PREFIX + "0", TEST_CONTIG, 15000, 0.45));
        for (int j = 0; j < 5; j++) {
            evidence.add(new BafEvidence(SAMPLE_PREFIX + j, TEST_CONTIG, 15000, 0.30 + j * 0.05));
        }
        final Set<String> carrierSamples = generatTestCarriers(1);

        final BafKolmogorovSmirnovTester tester = new BafKolmogorovSmirnovTester(minBafCount, seed);
        final BafKolmogorovSmirnovTester.KSTestResult result = tester.test(record, evidence, carrierSamples);
        Assert.assertNotNull(result, "Singleton carrier should produce a non-null KS result");
        Assert.assertTrue(result.getStat() > 0, "KS stat should be positive");
        Assert.assertTrue(result.getP() > 0 && result.getP() <= 1, "p-value should be in (0, 1]");
    }

    @Test
    public void testApplyToRecord() {
        final BafKolmogorovSmirnovTester tester = new BafKolmogorovSmirnovTester(2, 42L);
        final SVCallRecord record = SVTestUtils.newDepthCallRecordWithIntervalAndType("chr21", 10000, 20000, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP);
        // Null test
        final SVCallRecord resultNulll = tester.applyToRecord(record, null);
        SVTestUtils.assertEquals(resultNulll, record);

        // Non-null test
        final Map<String, Object> attr = new HashMap<>();
        attr.put(GATKSVVCFConstants.BAF_KS_STAT_ATTRIBUTE, 1.5);
        attr.put(GATKSVVCFConstants.BAF_KS_Q_ATTRIBUTE, 20.0);
        final SVCallRecord expected = SVCallRecordUtils.copyCallWithNewAttributes(record, attr);
        final SVCallRecord result = tester.applyToRecord(record, new BafKolmogorovSmirnovTester.KSTestResult(1.5, 0.013));
        SVTestUtils.assertEqualsExceptExcludedAttributes(result, expected, Lists.newArrayList(GATKSVVCFConstants.BAF_KS_STAT_ATTRIBUTE, GATKSVVCFConstants.BAF_KS_Q_ATTRIBUTE));
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.BAF_KS_Q_ATTRIBUTE));
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.BAF_KS_STAT_ATTRIBUTE));
        SVTestUtils.assertFloatWithinTolerance((Double) result.getAttributes().get(GATKSVVCFConstants.BAF_KS_Q_ATTRIBUTE), 1.8860566476931634, 1e-6);
        SVTestUtils.assertFloatWithinTolerance((Double) result.getAttributes().get(GATKSVVCFConstants.BAF_KS_STAT_ATTRIBUTE), 1.5, 1e-6);

        final SVCallRecord capped = tester.applyToRecord(record, new BafKolmogorovSmirnovTester.KSTestResult(0.5, 0.0));
        SVTestUtils.assertFloatWithinTolerance((Double) capped.getAttributes().get(GATKSVVCFConstants.BAF_KS_Q_ATTRIBUTE), 50.0, 1e-6);
    }
}