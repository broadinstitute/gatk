package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.walkers.sv.SVOverlap;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class PartitionedSVClusterEngine<T extends SVOverlap.PartitionedSVCallRecord> extends SVClusterEngine<T, PartitionedSVClusterEngine.PartitionedCluster> {

    protected final Map<Integer, T> primaryIdToItemMap; // Active primary items
    protected final Map<Integer, Integer> primaryIdToClusterIdMap;

    public PartitionedSVClusterEngine(final SVClusterLinkage<T> linkage,
                                      final SAMSequenceDictionary dictionary) {
        super(linkage, dictionary);
        primaryIdToItemMap = new HashMap<>();
        primaryIdToClusterIdMap = new HashMap<>();
    }

    @Override
    protected void createItemClusterComparator() {
        Utils.nonNull(itemComparator);
        final Comparator<PartitionedOutputCluster> comparator = (o1, o2) -> itemComparator.compare(o1.getPrimaryItem(), o2.getPrimaryItem());
        buffer = new ClusterSortingBuffer(comparator);
    }

    @Override
    protected PartitionedCluster createCluster(final List<Integer> itemIds) {
        final List<Integer> primaryItems = itemIds.stream().filter(id -> primaryIdToItemMap.containsKey(id)).collect(Collectors.toList());
        final List<Integer> nonPrimaryItemIds = itemIds.stream().filter(id -> !primaryIdToItemMap.containsKey(id)).collect(Collectors.toList());
        final List<T> items = itemIds.stream().map(this::getItem).collect(Collectors.toList());
        Utils.validate(primaryItems.size() == 1, "There are " + primaryItems.size() + " primary items in this cluster, but there should only be 1");
        return new PartitionedCluster(linkage.getMaxClusterableStartingPosition(items), nonPrimaryItemIds, primaryItems.get(0));
    }

    private final T getPrimaryItem(final int id) {
        Utils.validateArg(primaryIdToItemMap.containsKey(id), "Primary item ID " + id + " does not exist.");
        return primaryIdToItemMap.get(id);
    }

    protected Set<Integer> getClusterableItemIds() {
        final Set<Integer> set = new HashSet<>(getItemIds());
        set.addAll(primaryIdToItemMap.keySet());
        return set;
    }

    @Override
    protected void cluster(final Integer itemId) {
        final T item = getItem(itemId);
        if (item.isPrimaryVariant()) {
            final List<Integer> clusterItemIds = getItemIds().stream()
                    .filter(other -> linkage.areClusterable(item, getItem(other)))
                    .collect(Collectors.toList());
            final int maxPos = Math.max(linkage.getMaxClusterableStartingPosition(item), getMaxClusterableStartingPositionByIds(clusterItemIds));
            final PartitionedCluster cluster = new PartitionedCluster(maxPos, clusterItemIds, itemId);
            registerCluster(cluster);
        } else {
            final int itemMaxClusterableStart = linkage.getMaxClusterableStartingPosition(item);
            getClusters().stream()
                    .map(PartitionedCluster::getPrimaryItem)
                    .distinct()
                    .filter(other -> !other.equals(itemId) && linkage.areClusterable(item, getPrimaryItem(other)))
                    .map(primaryIdToClusterIdMap::get)
                    .map(this::getCluster)
                    .forEach(cluster -> cluster.addMember(itemId, itemMaxClusterableStart));
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

    public class PartitionedCluster extends SVClusterEngine.Cluster {

        final Integer primaryItem;

        public PartitionedCluster(final int maxClusterableStart, final List<Integer> members, final Integer primaryItem) {
            super(maxClusterableStart, members);
            Utils.nonNull(primaryItem);
            this.primaryItem = primaryItem;
        }

        public Integer getPrimaryItem() {
            return primaryItem;
        }
    }

    public class PartitionedOutputCluster extends SVClusterEngine.OutputCluster {

        final T primaryItem;

        public PartitionedOutputCluster(final List<T> members, final T primaryItem) {
            super(members);
            Utils.nonNull(primaryItem);
            this.primaryItem = primaryItem;
        }

        public T getPrimaryItem() {
            return primaryItem;
        }
    }

}
