package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.walkers.sv.CrossReferenceSVGenotypes;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class CrossRefSVClusterEngine<T extends CrossReferenceSVGenotypes.CrossRefSVCallRecord> extends SVClusterEngine<T, CrossReferenceSVGenotypes.CrossRefCluster, CrossReferenceSVGenotypes.CrossRefOutputCluster<T>> {

    protected final Map<Long, T> testIdToTestItemMap; // Active primary items to be annotated
    protected final Map<Long, Integer> testIdToClusterIdMap;

    public CrossRefSVClusterEngine(final SVClusterLinkage<T> linkage,
                                   final SAMSequenceDictionary dictionary) {
        super(linkage, dictionary);
        testIdToTestItemMap = new HashMap<>();
        testIdToClusterIdMap = new HashMap<>();
    }

    @Override
    protected CrossReferenceSVGenotypes.CrossRefOutputCluster<T> createOutputCluster(final CrossReferenceSVGenotypes.CrossRefCluster cluster) {
        return new CrossReferenceSVGenotypes.CrossRefOutputCluster<>(cluster.getMembers().stream().collect(Collectors.toMap(id -> id, this::getItem)), cluster.getTestItem(), getItem(cluster.getTestItem()));
    }

    @Override
    protected void putNewItem(final T item, final Long itemId) {
        if (item.isTestVariant()) {
            testIdToTestItemMap.put(itemId, item);
        } else {
            super.putNewItem(item, itemId);
        }
    }

    @Override
    protected final T getItem(final Long id) {
        if (testIdToTestItemMap.containsKey(id)) {
            return testIdToTestItemMap.get(id);
        }
        return super.getItem(id);
    }

    @Override
    protected void removeItem(final Long id) {
        if (testIdToTestItemMap.containsKey(id)) {
            testIdToTestItemMap.remove(id);
            testIdToClusterIdMap.remove(id);
        } else {
            super.removeItem(id);
        }
    }

    @Override
    protected List<CrossReferenceSVGenotypes.CrossRefOutputCluster<T>> hardFlush() {
        final List<CrossReferenceSVGenotypes.CrossRefOutputCluster<T>> output = super.hardFlush();
        testIdToClusterIdMap.clear();
        testIdToTestItemMap.clear();
        return output;
    }

    @Override
    public final Set<Long> getItemIds() {
        final Set<Long> set = new HashSet<>(super.getItemIds());
        set.addAll(testIdToTestItemMap.keySet());
        return set;
    }

    @Override
    protected void cluster(final Long itemId) {
        final T item = getItem(itemId);
        if (item.isTestVariant()) {
            final List<Long> linkedItemIds = getItemIds().stream()
                    .filter(other -> !getItem(other).isTestVariant() && linkage.areClusterable(item, getItem(other)))
                    .collect(Collectors.toList());
            final int maxPos;
            if (linkedItemIds.isEmpty()) {
                maxPos = linkage.getMaxClusterableStartingPosition(item);
            } else {
                maxPos = Math.max(linkage.getMaxClusterableStartingPosition(item), getMaxClusterableStartingPositionByIds(linkedItemIds));
            }
            final CrossReferenceSVGenotypes.CrossRefCluster cluster = new CrossReferenceSVGenotypes.CrossRefCluster(linkedItemIds, getCurrentContig(), maxPos, itemId);
            final Integer clusterId = registerCluster(cluster);
            testIdToClusterIdMap.put(itemId, clusterId);
        } else {
            final int itemMaxClusterableStart = linkage.getMaxClusterableStartingPosition(item);
            final String itemContig = item.getContigA();
            getClusters().stream()
                    .map(CrossReferenceSVGenotypes.CrossRefCluster::getTestItem)
                    .distinct()
                    .filter(other -> !other.equals(itemId) && linkage.areClusterable(item, getTestItem(other)))
                    .map(testIdToClusterIdMap::get)
                    .map(this::getCluster)
                    .forEach(cluster -> cluster.addMember(itemId, itemContig, itemMaxClusterableStart));
        }
    }

    private T getTestItem(final Long id) {
        Utils.validateArg(testIdToTestItemMap.containsKey(id), "Test item ID " + id + " does not exist.");
        return testIdToTestItemMap.get(id);
    }

}
