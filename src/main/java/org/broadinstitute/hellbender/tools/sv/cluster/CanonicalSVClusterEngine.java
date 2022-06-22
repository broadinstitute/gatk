package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;

import java.util.*;
import java.util.stream.Collectors;

public class CanonicalSVClusterEngine<T extends SVLocatable> extends SVClusterEngine<T, BasicCluster, BasicOutputCluster<T>> {

    protected final CLUSTERING_TYPE clusteringType;

    /**
     * @param clusteringType algorithm choice
     * @param linkage
     * @param dictionary
     */
    public CanonicalSVClusterEngine(CLUSTERING_TYPE clusteringType, SVClusterLinkage<T> linkage, SAMSequenceDictionary dictionary) {
        super(linkage, dictionary);
        this.clusteringType = clusteringType;
    }

    protected BasicCluster createCluster(final List<Long> itemIds) {
        return new BasicCluster(itemIds, getCurrentContig(), getMaxClusterableStartingPositionByIds(itemIds));
    }

    @Override
    protected BasicOutputCluster<T> createOutputCluster(final BasicCluster cluster) {
        return new BasicOutputCluster<>(cluster.getMembers().stream().collect(Collectors.toMap(id -> id, this::getItem)));
    }

    /**
     * Add a new {@param <T>} to the current clusters and determine which are complete
     * @param itemId id of registered item to add
     */
    @Override
    protected void cluster(final Long itemId) {
        final T item = getItem(itemId);
        // Get list of item IDs from active clusters that cluster with this item
        final Set<Long> linkedItems = getClusters().stream().map(BasicCluster::getMembers)
                .flatMap(Set::stream)
                .distinct()
                .filter(other -> !other.equals(itemId) && linkage.areClusterable(item, getItem(other)))
                .collect(Collectors.toCollection(LinkedHashSet::new));

        // Clusters to which we simply add the item
        final List<Integer> clustersToAugment = new ArrayList<>();
        // New clusters, formed from subsets of currently active clusters, to which we will add the item
        final Set<List<Long>> clustersToSeedWith = new HashSet<>();    // Use set to prevent creating duplicate clusters
        for (final Integer clusterIndex : getClusterIds()) {
            final BasicCluster cluster = getCluster(clusterIndex);
            final Set<Long> clusterItems = cluster.getMembers();
            if (item.getPositionA() <= cluster.getMaxClusterableStart()) {
                if (clusteringType.equals(CLUSTERING_TYPE.MAX_CLIQUE)) {
                    final List<Long> linkedClusterItems = clusterItems.stream().filter(linkedItems::contains).collect(Collectors.toList());
                    final int numLinkedItems = linkedClusterItems.size();
                    if (numLinkedItems == clusterItems.size()) {
                        clustersToAugment.add(clusterIndex);
                    } else if (numLinkedItems > 0) {
                        clustersToSeedWith.add(linkedClusterItems);
                    }
                } else if (clusteringType.equals(CLUSTERING_TYPE.SINGLE_LINKAGE)) {
                    final boolean matchesCluster = clusterItems.stream().anyMatch(linkedItems::contains);
                    if (matchesCluster) {
                        clustersToAugment.add(clusterIndex);
                    }
                } else {
                    throw new IllegalArgumentException("Clustering algorithm for type " + clusteringType.name() + " not implemented");
                }
            }
        }

        // Create new clusters from subsets (max-clique only)
        if (!clustersToSeedWith.isEmpty()) {
            // Currently existing clusters to which we will add the current item
            final List<Set<Long>> augmentedClusterItemLists = clustersToAugment.stream()
                    .map(this::getCluster)
                    .map(BasicCluster::getMembers)
                    .map(HashSet::new)
                    .collect(Collectors.toList());
            final List<Set<Long>> triggeredClusterItemSets = new ArrayList<>(clustersToSeedWith.size() + augmentedClusterItemLists.size());
            // New clusters formed from subsets of the currently active clusters
            triggeredClusterItemSets.addAll(clustersToSeedWith.stream().map(HashSet::new).collect(Collectors.toList()));
            triggeredClusterItemSets.addAll(augmentedClusterItemLists);
            triggeredClusterItemSets.sort(Comparator.comparingInt(Set::size));
            for (int i = 0; i < triggeredClusterItemSets.size(); i++) {
                final Set<Long> seedItems = triggeredClusterItemSets.get(i);
                // Check that this cluster is not a sub-cluster of any of the others being created
                boolean isSubset = false;
                for (int j = i + 1; j < triggeredClusterItemSets.size(); j++) {
                    if (triggeredClusterItemSets.get(j).containsAll(seedItems)) {
                        isSubset = true;
                        break;
                    }
                }
                if (!isSubset) {
                    seedWithExistingCluster(itemId, seedItems);
                }
            }
        }

        // Add to or merge existing clusters
        if (clusteringType.equals(CLUSTERING_TYPE.SINGLE_LINKAGE)) {
            if (!clustersToAugment.isEmpty()) {
                combineClusters(clustersToAugment, itemId);
            }
        } else {
            for (final Integer clusterId : clustersToAugment) {
                addToCluster(clusterId, itemId);
            }
        }

        // If there weren't any matches, create a new singleton cluster
        if (clustersToAugment.isEmpty() && clustersToSeedWith.isEmpty()) {
            seedCluster(itemId);
        }
    }

    private void seedCluster(final Long itemId) {
        final List<Long> members = new ArrayList<>(1);
        members.add(itemId);
        registerCluster(createCluster(members));
    }

    /**
     * Add the item specified by {@param itemId} to the cluster specified by {@param clusterIndex}
     * and expand the clustering interval
     * @param clusterId
     * @param itemId
     */
    private final void addToCluster(final int clusterId, final Long itemId) {
        final T item = getItem(itemId);
        getCluster(clusterId).addMember(itemId, item.getContigA(), linkage.getMaxClusterableStartingPosition(item));
    }

    /**
     * Creates a new cluster by agglomerating clusters with the given ids together (and deleting them from the currently
     * active set), along with an additional item.
     * @param clusterIds ids of clusters to combine
     * @param itemId id of item to add to new cluster
     */
    protected final void combineClusters(final Collection<Integer> clusterIds, final Long itemId) {
        final Collection<BasicCluster> clusters = removeClusters(clusterIds);
        final List<Long> clusterItems = clusters.stream()
                .map(BasicCluster::getMembers)
                .flatMap(Set::stream)
                .distinct()
                .collect(Collectors.toList());
        final List<Long> newClusterItems = new ArrayList<>(clusterItems.size() + 1);
        newClusterItems.addAll(clusterItems);
        newClusterItems.add(itemId);
        registerCluster(createCluster(newClusterItems));
    }

    /**
     * Create a new cluster containing a new item and a set of existing items, and add it to the set of active clusters.
     * Note that if there exists an identical activate cluster, it will not be added.
     * @param item new item (assumed registered)
     * @param seedItems existing items
     */
    protected final void seedWithExistingCluster(final Long item, final Collection<Long> seedItems) {
        final List<Long> newClusterItems = new ArrayList<>(1 + seedItems.size());
        newClusterItems.addAll(seedItems);
        newClusterItems.add(item);
        registerCluster(createCluster(newClusterItems));
    }
}
