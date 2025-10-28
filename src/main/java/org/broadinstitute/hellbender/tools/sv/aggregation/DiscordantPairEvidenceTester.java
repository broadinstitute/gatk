package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Analyzes and refines SV breakpoints using discordant pair (PE) evidence.
 * See {@link SplitReadEvidenceTester} for description of statistical methods.
 */

public class DiscordantPairEvidenceTester {

    private final Map<String,Double> sampleCoverageMap;
    private final SAMSequenceDictionary dictionary;

    public static final int MAX_QUAL = 99;

    /**
     * @param sampleCoverageMap map with (sample id, per-base sample coverage) entries
     * @param dictionary reference dictionary
     */
    public DiscordantPairEvidenceTester(final Map<String, Double> sampleCoverageMap,
                                        final SAMSequenceDictionary dictionary) {
        this.sampleCoverageMap = Utils.nonNull(sampleCoverageMap);
        this.dictionary = Utils.nonNull(dictionary);
    }

    /**
     * Applies statistical model to PE evidence around variant breakends.
     *
     * @param record  variant to test
     * @param evidence  PE evidence for the variant
     * @param carrierSamples  non-ref samples for this variant
     * @param backgroundSamples  hom-ref samples for this variant
     */
    public DiscordantPairTestResult test(final SVCallRecord record,
                                         final List<DiscordantPairEvidence> evidence,
                                         final Set<String> carrierSamples,
                                         final Set<String> backgroundSamples) {
        Utils.nonNull(record);
        SVCallRecordUtils.validateCoordinatesWithDictionary(record, dictionary);
        return poissonTest(evidence, carrierSamples, backgroundSamples);
    }

    /**
     * Annotates record with PE test results
     */
    public SVCallRecord applyToRecord(final SVCallRecord record, final DiscordantPairTestResult discordantPairResult) {
        Utils.nonNull(record);
        Utils.nonNull(discordantPairResult);
        final EvidenceStatUtils.PoissonTestResult test = discordantPairResult.getTest();
        final double p = test == null ? 1. : test.getP();
        final Double q = EvidenceStatUtils.probToQual(p, (byte) MAX_QUAL);
        final Double carrierSignal = test == null ? 0 :
                EvidenceStatUtils.carrierSignalFraction(test.getCarrierSignal(), test.getBackgroundSignal());
        final Map<String, Object> attributes = new HashMap<>(record.getAttributes());
        attributes.put(GATKSVVCFConstants.DISCORDANT_PAIR_QUALITY_ATTRIBUTE, q);
        attributes.put(GATKSVVCFConstants.DISCORDANT_PAIR_CARRIER_SIGNAL_ATTRIBUTE, carrierSignal);
        final SVCallRecord newRecord = SVCallRecordUtils.copyCallWithNewAttributes(record, attributes);
        return SVCallRecordUtils.assignDiscordantPairCountsToGenotypes(newRecord, discordantPairResult.getDiscordantPairEvidence());
    }

    public static Map<String, Integer> countEvidence(final List<DiscordantPairEvidence> evidence) {
        return evidence.stream()
                .collect(Collectors.groupingBy(DiscordantPairEvidence::getSample,
                        Collectors.collectingAndThen(Collectors.toList(), List::size)));
    }

    /**
     * Runs statistical tests
     */
    protected DiscordantPairTestResult poissonTest(final List<DiscordantPairEvidence> evidence,
                                                   final Collection<String> carrierSamples,
                                                   final Collection<String> backgroundSamples) {
        Utils.validateArg(sampleCoverageMap.keySet().containsAll(carrierSamples),
                "One or more carrier samples not found in sample coverage map");
        Utils.validateArg(sampleCoverageMap.keySet().containsAll(backgroundSamples),
                "One or more non-carrier samples not found in sample coverage map");
        final Map<String, Integer> sampleCounts = countEvidence(evidence);
        final EvidenceStatUtils.PoissonTestResult test = EvidenceStatUtils.calculateOneSamplePoissonTest(sampleCounts, carrierSamples, backgroundSamples,
                sampleCoverageMap, PESREvidenceTester.DEPTH_BASIS);
        return new DiscordantPairTestResult(test, sampleCounts, evidence);
    }

    public static final class DiscordantPairTestResult {
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

        /**
         * Gets stats
         */
        public EvidenceStatUtils.PoissonTestResult getTest() {
            return test;
        }

        /**
         * Gets map from sample ID to PE evidence count
         */
        public Map<String, Integer> getSampleCounts() {
            return sampleCounts;
        }

        /**
         * Gets all PE evidence for the test
         */
        public List<DiscordantPairEvidence> getDiscordantPairEvidence() {
            return discordantPairEvidence;
        }
    }
}
