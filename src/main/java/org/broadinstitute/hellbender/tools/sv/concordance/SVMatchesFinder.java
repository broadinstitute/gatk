package org.broadinstitute.hellbender.tools.sv.concordance;

import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class SVMatchesFinder implements SVMatcher {
    protected final Map<Long, SVCallRecord> truthIdToItemMap;
    protected final Map<Long, ActiveCluster> idToClusterMap;
    private final SVConcordanceLinkage linkage;

    private Integer lastItemStart;
    private String lastItemContig;

    public SVMatchesFinder(final SVConcordanceLinkage linkage) {
        this.linkage = Utils.nonNull(linkage);
        truthIdToItemMap = new HashMap<>();
        idToClusterMap = new HashMap<>();
        lastItemStart = null;
        lastItemContig = null;
    }

    private SVCallRecord annotate(final ActiveCluster cluster) {
        final Map<String, Object> attributes = new HashMap<>(cluster.getItem().getAttributes());
        final ConcordanceState variantStatus = cluster.getMatchVids().isEmpty() ? ConcordanceState.FALSE_POSITIVE : ConcordanceState.TRUE_POSITIVE;
        final List<String> matchVids = cluster.getMatchVids().isEmpty() ? null : cluster.getMatchVids();
        attributes.put(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE, variantStatus.getAbbreviation());
        attributes.put(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, matchVids);
        if (variantStatus.equals(ConcordanceState.TRUE_POSITIVE)) {
            final List<CanonicalSVLinkage.CanonicalLinkageResult> linkageResults = cluster.getLinkageResults();
            attributes.put(GATKSVVCFConstants.TRUTH_RECIPROCAL_OVERLAP_INFO, linkageResults.stream()
                    .map(CanonicalSVLinkage.CanonicalLinkageResult::getReciprocalOverlap)
                    .collect(Collectors.toList()));
            attributes.put(GATKSVVCFConstants.TRUTH_SIZE_SIMILARITY_INFO, linkageResults.stream()
                    .map(CanonicalSVLinkage.CanonicalLinkageResult::getSizeSimilarity)
                    .collect(Collectors.toList()));
            attributes.put(GATKSVVCFConstants.TRUTH_DISTANCE_START_INFO, linkageResults.stream()
                    .map(CanonicalSVLinkage.CanonicalLinkageResult::getBreakpointDistance1)
                    .collect(Collectors.toList()));
            attributes.put(GATKSVVCFConstants.TRUTH_DISTANCE_END_INFO, linkageResults.stream()
                    .map(CanonicalSVLinkage.CanonicalLinkageResult::getBreakpointDistance2)
                    .collect(Collectors.toList()));
            // TODO: add header line for log AF diff and annotate?
            attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO, cluster.getMatchAlleleCounts());
            attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO, cluster.getMatchAlleleFrequencies());
            attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO, cluster.getMatchAlleleNumbers());

        } else {
            attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO, null);
            attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO, null);
            attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO, null);
        }
        return SVCallRecordUtils.copyCallWithNewAttributes(cluster.getItem(), attributes);
    }

    public List<ClosestSVFinder.LinkageConcordanceRecord> flush(final boolean force) {
        // output LinkageConcordanceRecord to keep Long id with record
        // use first linkage as placeholder - not used in StratifiedConcordanceEngine
        return flushClusters(force).stream()
                .map(c -> new ClosestSVFinder.LinkageConcordanceRecord(c.getItemId(), annotate(c), c.getOneLinkageResult()))
                .collect(Collectors.toList());
    }


    /**
     * Flushes active clusters
     */
    private List<ActiveCluster> flushClusters(final boolean force) {
        if (force) {
            final List<ActiveCluster> output = new ArrayList<>(idToClusterMap.values());
            truthIdToItemMap.clear();
            idToClusterMap.clear();
            lastItemStart = null;
            lastItemContig = null;
            return output;
        } else {
            // Remove finalized ref items
            truthIdToItemMap.values().removeIf(v -> linkage.getMaxClusterableStartingPosition(v) < lastItemStart);
            // Find and remove finalized clusters
            final List<Map.Entry<Long, ActiveCluster>> finalizedClusters = idToClusterMap.entrySet().stream()
                    .filter(e -> e.getValue().getMaxClusterableStartingPosition() < lastItemStart)
                    .toList();
            finalizedClusters.forEach(e -> idToClusterMap.remove(e.getKey()));
            return finalizedClusters.stream().map(Map.Entry::getValue).collect(Collectors.toList());
        }
    }

    public String getLastItemContig() {
        return lastItemContig;
    }

    public void add(final SVCallRecord item, final Long id, final boolean isTruthVariant) {
        Utils.validateArg(lastItemContig == null || lastItemContig.equals(item.getContigA()), "Attempted to add item on a new contig; please run a force flush beforehand");
        Utils.validateArg(lastItemStart == null || lastItemStart <= item.getPositionA(), "Items must be added in dictionary-sorted order");
        Utils.validateArg(!idToClusterMap.containsKey(id), "ID already in use: " + id);
        lastItemContig = item.getContigA();
        lastItemStart = item.getPositionA();
        if (isTruthVariant) {
            Utils.validateArg(!truthIdToItemMap.containsKey(id), "ID already in use: " + id);
            truthIdToItemMap.put(id, item);
            for (final ActiveCluster cluster : idToClusterMap.values()) {
                final CanonicalSVLinkage.CanonicalLinkageResult result = linkage.areClusterable(cluster.getItem(), item);
                if (result.getResult()) {
                    cluster.update(item.getId(), result);
                }
            }
        } else {
            final int maxStart = linkage.getMaxClusterableStartingPosition(item);
            final ActiveCluster cluster = new ActiveCluster(id, item, new ArrayList<>(), new ArrayList<>(),
                    new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), maxStart);
            for (final SVCallRecord truthItem : truthIdToItemMap.values()) {
                final CanonicalSVLinkage.CanonicalLinkageResult result = linkage.areClusterable(item, truthItem);
                if (result.getResult()) {
                    cluster.update(truthItem.getId(), result);
                }
            }
            idToClusterMap.put(id, cluster);
        }
    }

    public static Double computeLogAlleleFrequencyDifference(final SVCallRecord a, final SVCallRecord b,
                                                             final double maxAlleleNumberA,
                                                             final double maxAlleleNumberB) {
        final Map<String, Object> attrA = a.getAttributes();
        final Map<String, Object> attrB = b.getAttributes();
        if (attrA.containsKey(VCFConstants.ALLELE_FREQUENCY_KEY)
                && (attrA.get(VCFConstants.ALLELE_FREQUENCY_KEY) != null)
                && attrB.containsKey(VCFConstants.ALLELE_FREQUENCY_KEY)
                && attrB.get(VCFConstants.ALLELE_FREQUENCY_KEY) != null) {
            double alleleFrequencyA = (Double) attrA.get(VCFConstants.ALLELE_FREQUENCY_KEY);
            double alleleFrequencyB = (Double) attrB.get(VCFConstants.ALLELE_FREQUENCY_KEY);
            // handle zeroes before taking log: if AF is 0, set to just below the minimum cohort AF
            if (alleleFrequencyA == 0) {
                alleleFrequencyA = 1 / (maxAlleleNumberA + 1);
            }
            if (alleleFrequencyB == 0) {
                alleleFrequencyB = 1 / (maxAlleleNumberB + 1);
            }
            return Math.abs(Math.log10(alleleFrequencyA) - Math.log10(alleleFrequencyB));
        } else {
            return null;
        }
    }

    public static class ActiveCluster {
        final Long itemId;
        final SVCallRecord item;
        final List<String> matchVids;
        final List<CanonicalSVLinkage.CanonicalLinkageResult> linkageResults;
        final List<Integer> matchAlleleCounts;
        final List<Integer> matchAlleleNumbers;
        final List<Double> matchAlleleFrequencies;
        final int maxClusterableStartingPosition;

        public ActiveCluster(final Long itemId, final SVCallRecord item, final List<String> matchVids,
                             final List<CanonicalSVLinkage.CanonicalLinkageResult> linkageResults,
                             final List<Integer> matchAlleleCounts, final List<Integer> matchAlleleNumbers,
                             final List<Double> matchAlleleFrequencies, int maxClusterableStartingPosition) {
            this.itemId = itemId;
            this.item = item;
            this.matchVids = matchVids;
            this.linkageResults = linkageResults;
            this.matchAlleleCounts = matchAlleleCounts;
            this.matchAlleleNumbers = matchAlleleNumbers;
            this.matchAlleleFrequencies = matchAlleleFrequencies;
            this.maxClusterableStartingPosition = maxClusterableStartingPosition;
        }

        void update(String matchVid, CanonicalSVLinkage.CanonicalLinkageResult linkageResult) {
            matchVids.add(matchVid);
            linkageResults.add(linkageResult);
        }

        Long getItemId() {
            return itemId;
        }

        SVCallRecord getItem() {
            return item;
        }

        List<String> getMatchVids() {
            return matchVids;
        }

        List<CanonicalSVLinkage.CanonicalLinkageResult> getLinkageResults() {
            return linkageResults;
        }

        CanonicalSVLinkage.CanonicalLinkageResult getOneLinkageResult() {
            if (getLinkageResults().isEmpty()) {
                return null;
            } else {
                return getLinkageResults().get(0);
            }
        }

        List<Integer> getMatchAlleleCounts() { return matchAlleleCounts; }

        List<Integer> getMatchAlleleNumbers() { return matchAlleleNumbers; }

        List<Double> getMatchAlleleFrequencies() { return matchAlleleFrequencies; }

        int getMaxClusterableStartingPosition() {
            return maxClusterableStartingPosition;
        }
    }

//    public static class FederationLinkageResult extends CanonicalSVLinkage.CanonicalLinkageResult {
//        private final Double logAlleleFrequencyDifference;
//
//        public FederationLinkageResult(final boolean result, final Double reciprocalOverlap,
//                                       final Double sizeSimilarity,
//                                       final Integer breakpointDistance1,
//                                       final Integer breakpointDistance2,
//                                       final Double logAlleleFrequencyDifference) {
//            super(result, reciprocalOverlap, sizeSimilarity, breakpointDistance1, breakpointDistance2);
//            this.logAlleleFrequencyDifference = logAlleleFrequencyDifference;
//        }
//
//        public FederationLinkageResult(final CanonicalSVLinkage.CanonicalLinkageResult result,
//                                       final Double logAlleleFrequencyDifference) {
//            super(result.getResult(), result.getReciprocalOverlap(), result.getSizeSimilarity(), result.getBreakpointDistance1(), result.getBreakpointDistance2());
//            this.logAlleleFrequencyDifference = logAlleleFrequencyDifference;
//        }
//
//        public Double getLogAlleleFrequencyDifference() { return this.logAlleleFrequencyDifference; }
//    }
}
