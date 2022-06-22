package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class PartitionedSVClusterEngine<T extends PartitionedSVCallRecord> extends SVClusterEngine<T, PartitionedCluster, PartitionedOutputCluster<T>> {

    protected final Map<Long, T> primaryIdToItemMap; // Active primary items
    protected final Map<Integer, Integer> primaryIdToClusterIdMap;

    public PartitionedSVClusterEngine(final SVClusterLinkage<T> linkage,
                                      final SAMSequenceDictionary dictionary) {
        super(linkage, dictionary);
        primaryIdToItemMap = new HashMap<>();
        primaryIdToClusterIdMap = new HashMap<>();
    }

    @Override
    protected PartitionedOutputCluster<T> createOutputCluster(final PartitionedCluster cluster) {
        return new PartitionedOutputCluster<>(cluster.getMembers().stream().collect(Collectors.toMap(id -> id, this::getItem)), getItem(cluster.getPrimaryItem()));
    }

    private T getPrimaryItem(final long id) {
        Utils.validateArg(primaryIdToItemMap.containsKey(id), "Primary item ID " + id + " does not exist.");
        return primaryIdToItemMap.get(id);
    }

    @Override
    protected void cluster(final Long itemId) {
        final T item = getItem(itemId);
        if (item.isPrimaryVariant()) {
            final List<Long> clusterItemIds = getItemIds().stream()
                    .filter(other -> linkage.areClusterable(item, getItem(other)))
                    .collect(Collectors.toList());
            final int maxPos = Math.max(linkage.getMaxClusterableStartingPosition(item), getMaxClusterableStartingPositionByIds(clusterItemIds));
            final int minPos = Math.min(item.getPositionA(), getMinStartingPositionByIds(clusterItemIds));
            final PartitionedCluster cluster = new PartitionedCluster(clusterItemIds, getCurrentContig(), maxPos, itemId);
            registerCluster(cluster);
        } else {
            final int itemMaxClusterableStart = linkage.getMaxClusterableStartingPosition(item);
            final String itemContig = item.getContigA();
            final int itemStart = item.getPositionA();
            getClusters().stream()
                    .map(PartitionedCluster::getPrimaryItem)
                    .distinct()
                    .filter(other -> !other.equals(itemId) && linkage.areClusterable(item, getPrimaryItem(other)))
                    .map(primaryIdToClusterIdMap::get)
                    .map(this::getCluster)
                    .forEach(cluster -> cluster.addMember(itemId, itemContig, itemMaxClusterableStart));
        }
    }

    protected void putNewItem(final long itemId, final T item) {
        if (item.isPrimaryVariant()) {
            primaryIdToItemMap.put(itemId, item);
        } else {
            super.putNewItem(item, itemId);
        }
    }

}
