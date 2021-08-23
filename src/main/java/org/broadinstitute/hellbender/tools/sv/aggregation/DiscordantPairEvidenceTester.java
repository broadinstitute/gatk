package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class DiscordantPairEvidenceTester {

    private final Map<String,Double> sampleCoverageMap;
    private final SAMSequenceDictionary dictionary;

    public DiscordantPairEvidenceTester(final Map<String, Double> sampleCoverageMap,
                                        final SAMSequenceDictionary dictionary) {
        this.sampleCoverageMap = Utils.nonNull(sampleCoverageMap);
        this.dictionary = Utils.nonNull(dictionary);
    }

    public DiscordantPairTestResult poissonTestRecord(final SVCallRecord record,
                                                      final List<DiscordantPairEvidence> evidence,
                                                      final Set<String> carrierSamples,
                                                      final Set<String> backgroundSamples) {
        Utils.nonNull(record);
        SVCallRecordUtils.validateCoordinatesWithDictionary(record, dictionary);

        final int representativeDepth = EvidenceStatUtils.computeRepresentativeDepth(sampleCoverageMap.values());
        return poissonTest(evidence, carrierSamples, backgroundSamples, representativeDepth);
    }

    public SVCallRecord applyToRecord(final SVCallRecord record, final DiscordantPairTestResult discordantPairResult) {
        Utils.nonNull(record);
        if (discordantPairResult == null) {
            return record;
        }
        final EvidenceStatUtils.PoissonTestResult test = discordantPairResult.getTest();
        final double p = test == null ? 1. : test.getP();
        final Integer q = EvidenceStatUtils.probToQual(p, (byte) 99);
        final Integer carrierSignal = test == null ? 0 :
                EvidenceStatUtils.carrierSignalFraction(test.getCarrierSignal(), test.getBackgroundSignal());
        final Map<String, Object> attributes = new HashMap<>();
        attributes.put(GATKSVVCFConstants.DISCORDANT_PAIR_QUALITY_ATTRIBUTE, q);
        attributes.put(GATKSVVCFConstants.DISCORDANT_PAIR_CARRIER_SIGNAL_ATTRIBUTE, carrierSignal);
        final SVCallRecord newRecord = SVCallRecordUtils.copyCallWithNewAttributes(record, attributes);
        return SVCallRecordUtils.assignDiscordantPairCountsToGenotypes(newRecord, discordantPairResult.getDiscordantPairEvidence());
    }

    public DiscordantPairTestResult poissonTest(final List<DiscordantPairEvidence> evidence,
                                                final Collection<String> carrierSamples,
                                                final Collection<String> backgroundSamples,
                                                final int representativeDepth) {
        Utils.validateArg(sampleCoverageMap.keySet().containsAll(carrierSamples),
                "One or more carrier samples not found in sample coverage map");
        Utils.validateArg(sampleCoverageMap.keySet().containsAll(backgroundSamples),
                "One or more non-carrier samples not found in sample coverage map");

        // Default case
        if (evidence.isEmpty() || carrierSamples.isEmpty() || backgroundSamples.isEmpty()) {
            return null;
        }
        final Map<String, Integer> sampleCounts = evidence.stream()
                .collect(Collectors.groupingBy(DiscordantPairEvidence::getSample,
                        Collectors.collectingAndThen(Collectors.toList(), List::size)));
        final EvidenceStatUtils.PoissonTestResult test = EvidenceStatUtils.calculateOneSamplePoissonTest(sampleCounts, carrierSamples, backgroundSamples,
                sampleCoverageMap, representativeDepth);
        return new DiscordantPairTestResult(test, sampleCounts, evidence);
    }

    public final class DiscordantPairTestResult {
        private final EvidenceStatUtils.PoissonTestResult test;
        private final Map<String, Integer> sampleCounts;
        private final List<DiscordantPairEvidence> discordantPairEvidence;
        public DiscordantPairTestResult(final EvidenceStatUtils.PoissonTestResult test,
                                        final Map<String, Integer> sampleCounts,
                                        final List<DiscordantPairEvidence> discordantPairEvidence) {
            this.test = test;
            this.sampleCounts = sampleCounts;
            this.discordantPairEvidence = discordantPairEvidence;
        }

        public EvidenceStatUtils.PoissonTestResult getTest() {
            return test;
        }

        public Map<String, Integer> getSampleCounts() {
            return sampleCounts;
        }

        public List<DiscordantPairEvidence> getDiscordantPairEvidence() {
            return discordantPairEvidence;
        }
    }
}
