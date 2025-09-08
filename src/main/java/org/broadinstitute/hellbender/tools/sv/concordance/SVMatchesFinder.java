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

public class SVMatchesFinder implements CorrespondingSVSelector {
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
                    cluster.update(item.getId(), result, item);
                }
            }
        } else {
            final int maxStart = linkage.getMaxClusterableStartingPosition(item);
            final ActiveCluster cluster = new ActiveCluster(id, item, new ArrayList<>(), new ArrayList<>(),
                    new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), maxStart);
            for (final SVCallRecord truthItem : truthIdToItemMap.values()) {
                final CanonicalSVLinkage.CanonicalLinkageResult result = linkage.areClusterable(item, truthItem);
                if (result.getResult()) {
                    cluster.update(truthItem.getId(), result, truthItem);
                }
            }
            idToClusterMap.put(id, cluster);
        }
    }

    public static class ActiveCluster {
        final Long itemId;
        final SVCallRecord item;
        final List<String> matchVids;
        final List<CanonicalSVLinkage.CanonicalLinkageResult> linkageResults;
        final List<Object> matchAlleleCounts;
        final List<Object> matchAlleleNumbers;
        final List<Object> matchAlleleFrequencies;
        final int maxClusterableStartingPosition;

        public ActiveCluster(final Long itemId, final SVCallRecord item, final List<String> matchVids,
                             final List<CanonicalSVLinkage.CanonicalLinkageResult> linkageResults,
                             final List<Object> matchAlleleCounts, final List<Object> matchAlleleNumbers,
                             final List<Object> matchAlleleFrequencies, int maxClusterableStartingPosition) {
            this.itemId = itemId;
            this.item = item;
            this.matchVids = matchVids;
            this.linkageResults = linkageResults;
            this.matchAlleleCounts = matchAlleleCounts;
            this.matchAlleleNumbers = matchAlleleNumbers;
            this.matchAlleleFrequencies = matchAlleleFrequencies;
            this.maxClusterableStartingPosition = maxClusterableStartingPosition;
        }

        void update(final String matchVid,
                    final CanonicalSVLinkage.CanonicalLinkageResult linkageResult,
                    final SVCallRecord matchRecord) {
            matchVids.add(matchVid);
            linkageResults.add(linkageResult);
            final Map<String, Object> matchAttr = matchRecord.getAttributes();
            matchAlleleCounts.add(matchAttr.get(VCFConstants.ALLELE_COUNT_KEY));
            matchAlleleNumbers.add(matchAttr.get(VCFConstants.ALLELE_NUMBER_KEY));
            matchAlleleFrequencies.add(matchAttr.get(VCFConstants.ALLELE_FREQUENCY_KEY));
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

        List<Object> getMatchAlleleCounts() { return matchAlleleCounts; }

        List<Object> getMatchAlleleNumbers() { return matchAlleleNumbers; }

        List<Object> getMatchAlleleFrequencies() { return matchAlleleFrequencies; }

        int getMaxClusterableStartingPosition() {
            return maxClusterableStartingPosition;
        }
    }

}
