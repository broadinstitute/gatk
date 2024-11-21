package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class PESREvidenceTester {

    public static final int MAX_QUAL = 99;

    private final Map<String,Double> sampleCoverageMap;
    private final SAMSequenceDictionary dictionary;
    protected int representativeDepth;

    /**
     * @param sampleCoverageMap map with (sample id, per-base sample coverage) entries
     * @param dictionary reference dictionary
     */
    public PESREvidenceTester(final Map<String, Double> sampleCoverageMap, final SAMSequenceDictionary dictionary) {
        this.sampleCoverageMap = Utils.nonNull(sampleCoverageMap);
        this.dictionary = Utils.nonNull(dictionary);
        this.representativeDepth = EvidenceStatUtils.computeRepresentativeDepth(sampleCoverageMap.values());
    }


    /**
     * Test combined PE/SR evidence
     */
    public PESRTestResult test(final SplitReadEvidenceTester.SplitReadTestResult splitReadTestResult,
                               final DiscordantPairEvidenceTester.DiscordantPairTestResult discordantPairTestResult,
                               final Set<String> carrierSamples,
                               final Set<String> backgroundSamples) {
        final Map<String, Integer> sampleCountSums = new HashMap<>(SVUtils.hashMapCapacity(carrierSamples.size() + backgroundSamples.size()));
        final Map<String, Integer> discordantPairCounts = discordantPairTestResult.getSampleCounts();
        final SplitReadSite startSite = splitReadTestResult.getFirst();
        final SplitReadSite endSite = splitReadTestResult.getSecond();
        final EvidenceStatUtils.PoissonTestResult bothsideResult = splitReadTestResult.getBothsidesResult();
        for (final String sample : Sets.union(carrierSamples, backgroundSamples)) {
            sampleCountSums.put(sample, startSite.getCount(sample) + endSite.getCount(sample) + discordantPairCounts.getOrDefault(sample, 0));
        }
        final EvidenceStatUtils.PoissonTestResult result =  EvidenceStatUtils.calculateOneSamplePoissonTest(sampleCountSums,
                carrierSamples, backgroundSamples, sampleCoverageMap, representativeDepth);

        final Integer combinedCarrierSignal = EvidenceStatUtils.carrierSignalFraction(
                discordantPairTestResult.getTest().getCarrierSignal() + bothsideResult.getCarrierSignal(),
                discordantPairTestResult.getTest().getBackgroundSignal() + bothsideResult.getBackgroundSignal());

        return new PESRTestResult(result, combinedCarrierSignal);
    }

    /**
     * Annotates record with PESR test results
     */
    public SVCallRecord applyToRecord(final SVCallRecord record,
                                      final PESRTestResult result) {
        Utils.nonNull(record);
        Utils.nonNull(result);
        final EvidenceStatUtils.PoissonTestResult discordantPairTest = result.getPesrResult();
        final Map<String, Object> refinedAttr = new HashMap<>(record.getAttributes());
        refinedAttr.put(GATKSVVCFConstants.PESR_CARRIER_SIGNAL_ATTRIBUTE, result.getCombinedCarrierSignal());
        final Integer pesrQuality = Double.isNaN(discordantPairTest.getP()) ?
                null : EvidenceStatUtils.probToQual(discordantPairTest.getP(), (byte) MAX_QUAL);
        refinedAttr.put(GATKSVVCFConstants.PESR_QUALITY_ATTRIBUTE, pesrQuality);
        return SVCallRecordUtils.copyCallWithNewAttributes(record, refinedAttr);
    }

    public final class PESRTestResult {
        private final EvidenceStatUtils.PoissonTestResult pesrResult;
        private final Integer combinedCarrierSignal;

        public PESRTestResult(final EvidenceStatUtils.PoissonTestResult pesrResult,
                              final Integer combinedCarrierSignal) {
            this.pesrResult = pesrResult;
            this.combinedCarrierSignal = combinedCarrierSignal;
        }

        public EvidenceStatUtils.PoissonTestResult getPesrResult() {
            return pesrResult;
        }

        public Integer getCombinedCarrierSignal() {
            return combinedCarrierSignal;
        }
    }
}
