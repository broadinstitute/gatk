package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.mockito.internal.util.collections.Sets;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public class PESREvidenceTesterTest {

    private static final String TEST_CONTIG = "chr21";
    private static final String SAMPLE_1 = "sample1";
    private static final String SAMPLE_2 = "sample2";
    private static final String SAMPLE_3 = "sample3";

    @Test
    public void test() {
        final Map<String, Double> sampleCoverageMap = new HashMap<>();
        sampleCoverageMap.put(SAMPLE_1, 30.);
        sampleCoverageMap.put(SAMPLE_2, 25.);
        sampleCoverageMap.put(SAMPLE_3, 35.);
        final Map<String, Integer> splitReadCountStart = new HashMap<>();
        splitReadCountStart.put(SAMPLE_1, 4);
        splitReadCountStart.put(SAMPLE_2, 2);
        final Map<String, Integer> splitReadCountEnd = new HashMap<>();
        splitReadCountStart.put(SAMPLE_1, 2);
        splitReadCountStart.put(SAMPLE_3, 1);
        final EvidenceStatUtils.PoissonTestResult splitReadResultStart = new EvidenceStatUtils.PoissonTestResult(0.01, 0.2, 0.3);
        final EvidenceStatUtils.PoissonTestResult splitReadResultEnd = new EvidenceStatUtils.PoissonTestResult(0.02, 0.1, 0.15);
        final EvidenceStatUtils.PoissonTestResult bothsideResult = new EvidenceStatUtils.PoissonTestResult(0.05, 0.6, 0.7);
        final SplitReadSite firstSite = new SplitReadSite(TEST_CONTIG, 1000, true, splitReadCountStart, splitReadResultStart);
        final SplitReadSite secondSite = new SplitReadSite(TEST_CONTIG, 2000, false, splitReadCountEnd, splitReadResultEnd);
        final SplitReadEvidenceTester.SplitReadTestResult splitReadTestResult = new SplitReadEvidenceTester.SplitReadTestResult(firstSite, secondSite, bothsideResult);

        final Map<String, Integer> discordantPairCounts = new HashMap<>();
        discordantPairCounts.put(SAMPLE_1, 8);
        discordantPairCounts.put(SAMPLE_3, 3);
        final EvidenceStatUtils.PoissonTestResult discordantPairResult = new EvidenceStatUtils.PoissonTestResult(0.04, 1.4, 0.98);
        final List<DiscordantPairEvidence> discordantPairEvidenceList = new ArrayList<>();
        final DiscordantPairEvidenceTester.DiscordantPairTestResult discordantPairTestResult = new DiscordantPairEvidenceTester.DiscordantPairTestResult(discordantPairResult, discordantPairCounts, discordantPairEvidenceList);

        final Set<String> carrierSamples = Collections.singleton(SAMPLE_1);
        final Set<String> backgroundSamples = Sets.newSet(SAMPLE_2, SAMPLE_3);

        final PESREvidenceTester tester = new PESREvidenceTester(sampleCoverageMap);
        final PESREvidenceTester.PESRTestResult result = tester.test(splitReadTestResult, discordantPairTestResult, carrierSamples, backgroundSamples);

        final double tolerance = 1e-6;
        SVTestUtils.assertFloatWithinTolerance(result.getPesrResult().getP(), 2.5512249585630067E-4, tolerance);
        SVTestUtils.assertFloatWithinTolerance(result.getPesrResult().getCarrierSignal(), 20., tolerance);
        SVTestUtils.assertFloatWithinTolerance(result.getPesrResult().getBackgroundSignal(), 6., tolerance);
        SVTestUtils.assertFloatWithinTolerance(result.getCombinedCarrierSignal(), 54.34782608695652, tolerance);
    }

    @Test
    public void testEmpty() {
        final Map<String, Double> sampleCoverageMap = new HashMap<>();
        final Map<String, Integer> splitReadCountStart = new HashMap<>();
        final Map<String, Integer> splitReadCountEnd = new HashMap<>();
        final EvidenceStatUtils.PoissonTestResult splitReadResultStart = new EvidenceStatUtils.PoissonTestResult(Double.NaN, 0, 0);
        final EvidenceStatUtils.PoissonTestResult splitReadResultEnd = new EvidenceStatUtils.PoissonTestResult(Double.NaN, 0, 0);
        final EvidenceStatUtils.PoissonTestResult bothsideResult = new EvidenceStatUtils.PoissonTestResult(Double.NaN, 0, 0);
        final SplitReadSite firstSite = new SplitReadSite(TEST_CONTIG, 1000, true, splitReadCountStart, splitReadResultStart);
        final SplitReadSite secondSite = new SplitReadSite(TEST_CONTIG, 2000, false, splitReadCountEnd, splitReadResultEnd);
        final SplitReadEvidenceTester.SplitReadTestResult splitReadTestResult = new SplitReadEvidenceTester.SplitReadTestResult(firstSite, secondSite, bothsideResult);

        final Map<String, Integer> discordantPairCounts = new HashMap<>();
        final EvidenceStatUtils.PoissonTestResult discordantPairResult = new EvidenceStatUtils.PoissonTestResult(Double.NaN, 0, 0);
        final List<DiscordantPairEvidence> discordantPairEvidenceList = new ArrayList<>();
        final DiscordantPairEvidenceTester.DiscordantPairTestResult discordantPairTestResult = new DiscordantPairEvidenceTester.DiscordantPairTestResult(discordantPairResult, discordantPairCounts, discordantPairEvidenceList);

        final Set<String> carrierSamples = Collections.emptySet();
        final Set<String> backgroundSamples = Collections.emptySet();

        final PESREvidenceTester tester = new PESREvidenceTester(sampleCoverageMap);
        final PESREvidenceTester.PESRTestResult result = tester.test(splitReadTestResult, discordantPairTestResult, carrierSamples, backgroundSamples);

        Assert.assertEquals(result.getPesrResult().getP(), 1.);
        Assert.assertEquals(result.getPesrResult().getCarrierSignal(), 0.);
        Assert.assertEquals(result.getPesrResult().getBackgroundSignal(), 0.);
        Assert.assertEquals(result.getCombinedCarrierSignal(), 0.);
    }

    @Test
    public void testApplyToRecord() {
        final Map<String, Double> sampleCoverageMap = new HashMap<>();
        final PESREvidenceTester tester = new PESREvidenceTester(sampleCoverageMap);
        final SVCallRecord record = SVTestUtils.newDepthCallRecordWithIntervalAndType("chr21", 10000, 20000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final Map<String, Object> attr = new HashMap<>();
        attr.put(GATKSVVCFConstants.PESR_QUALITY_ATTRIBUTE, 20.0);
        attr.put(GATKSVVCFConstants.PESR_CARRIER_SIGNAL_ATTRIBUTE, 0.1);
        final PESREvidenceTester.PESRTestResult tesetResult = new PESREvidenceTester.PESRTestResult(new EvidenceStatUtils.PoissonTestResult(0.11, 0.4, 0.5), 0.8);
        final SVCallRecord result = tester.applyToRecord(record, tesetResult);
        SVTestUtils.assertEqualsExceptExcludedAttributes(result, record, Lists.newArrayList(GATKSVVCFConstants.PESR_QUALITY_ATTRIBUTE, GATKSVVCFConstants.PESR_CARRIER_SIGNAL_ATTRIBUTE));
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.PESR_QUALITY_ATTRIBUTE));
        Assert.assertTrue(result.getAttributes().containsKey(GATKSVVCFConstants.PESR_CARRIER_SIGNAL_ATTRIBUTE));
        SVTestUtils.assertFloatWithinTolerance((Double) result.getAttributes().get(GATKSVVCFConstants.PESR_QUALITY_ATTRIBUTE),9.58607314841775, 1e-6);
        SVTestUtils.assertFloatWithinTolerance((Double) result.getAttributes().get(GATKSVVCFConstants.PESR_CARRIER_SIGNAL_ATTRIBUTE), 0.8, 1e-6);
    }

}