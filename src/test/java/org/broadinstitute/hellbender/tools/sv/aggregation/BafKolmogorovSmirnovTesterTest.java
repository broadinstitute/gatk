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