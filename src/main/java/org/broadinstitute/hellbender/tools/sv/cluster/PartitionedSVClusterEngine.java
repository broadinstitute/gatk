package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class PartitionedSVClusterEngine<T extends PartitionedSVCallRecord> extends SVClusterEngine<T, PartitionedCluster, PartitionedOutputCluster<T>> {

    protected final Map<Integer, T> primaryIdToItemMap; // Active primary items
    protected final Map<Integer, Integer> primaryIdToClusterIdMap;

    public PartitionedSVClusterEngine(final SVClusterLinkage<T> linkage,
                                      final SAMSequenceDictionary dictionary) {
        super(linkage, dictionary);
        primaryIdToItemMap = new HashMap<>();
        primaryIdToClusterIdMap = new HashMap<>();
    }

    @Override
    protected ClusterSortingBuffer<PartitionedOutputCluster<T>> createClusterSortingBuffer() {
        Utils.nonNull(itemComparator);
        final Comparator<PartitionedOutputCluster<T>> comparator = (o1, o2) -> itemComparator.compare(o1.getPrimaryItem(), o2.getPrimaryItem());
        return new ClusterSortingBuffer<>(comparator);
    }

    @Override
    protected PartitionedOutputCluster<T> createOutputCluster(final PartitionedCluster cluster) {
        return new PartitionedOutputCluster<>(cluster.getMembers().stream().map(this::getItem).collect(Collectors.toList()), getItem(cluster.getPrimaryItem()));
    }

    @Override
    protected PartitionedOutputCluster<T> createSingleItemOutputClusterForBuffer(final T item) {
        return new PartitionedOutputCluster<>(Collections.emptyList(), item);
    }

    private final T getPrimaryItem(final int id) {
        Utils.validateArg(primaryIdToItemMap.containsKey(id), "Primary item ID " + id + " does not exist.");
        return primaryIdToItemMap.get(id);
    }

    @Override
    protected void cluster(final Integer itemId) {
        final T item = getItem(itemId);
        if (item.isPrimaryVariant()) {
            final List<Integer> clusterItemIds = getItemIds().stream()
                    .filter(other -> linkage.areClusterable(item, getItem(other)))
                    .collect(Collectors.toList());
            final int maxPos = Math.max(linkage.getMaxClusterableStartingPosition(item), getMaxClusterableStartingPositionByIds(clusterItemIds));
            final int minPos = Math.min(item.getPositionA(), getMinStartingPositionByIds(clusterItemIds));
            final PartitionedCluster cluster = new PartitionedCluster(clusterItemIds, minPos, maxPos, itemId);
            registerCluster(cluster);
        } else {
            final int itemMaxClusterableStart = linkage.getMaxClusterableStartingPosition(item);
            final int itemStart = item.getPositionA();
            getClusters().stream()
                    .map(PartitionedCluster::getPrimaryItem)
                    .distinct()
                    .filter(other -> !other.equals(itemId) && linkage.areClusterable(item, getPrimaryItem(other)))
                    .map(primaryIdToClusterIdMap::get)
                    .map(this::getCluster)
                    .forEach(cluster -> cluster.addMember(itemId, itemStart, itemMaxClusterableStart));
        }
    }

    @Override
    protected void putNewItem(final Integer itemId, final T item) {
        if (item.isPrimaryVariant()) {
            primaryIdToItemMap.put(itemId, item);
        } else {
            super.putNewItem(itemId, item);
        }
    }

}
