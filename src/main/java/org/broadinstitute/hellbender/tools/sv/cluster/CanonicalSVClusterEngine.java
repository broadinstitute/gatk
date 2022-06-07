package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class CanonicalSVClusterEngine<T extends SVLocatable> extends SVClusterEngine<T, SVClusterEngine.Cluster> {

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

    @Override
    protected void createItemClusterComparator() {
        Utils.nonNull(itemComparator);
        final Comparator<OutputCluster> comparator = (o1, o2) -> {
            final T min1 = Collections.min(o1.getMembers(), itemComparator);
            final T min2 = Collections.min(o2.getMembers(), itemComparator);
            return itemComparator.compare(min1, min2);
        };
        buffer = new ClusterSortingBuffer(comparator);
    }

    @Override
    protected Cluster createCluster(final List<Integer> itemIds) {
        final List<T> items = itemIds.stream().map(this::getItem).collect(Collectors.toList());
        return new Cluster(linkage.getMaxClusterableStartingPosition(items), itemIds);
    }


    /**
     * Add a new {@param <T>} to the current clusters and determine which are complete
     * @param itemId id of registered item to add
     */
    @Override
    protected void cluster(final Integer itemId) {
        final T item = getItem(itemId);
        // Get list of item IDs from active clusters that cluster with this item
        final Set<Integer> linkedItems = getClusters().stream().map(Cluster::getMembers)
                .flatMap(List::stream)
                .distinct()
                .filter(other -> !other.equals(itemId) && linkage.areClusterable(item, getItem(other)))
                .collect(Collectors.toCollection(LinkedHashSet::new));

        // Find clusters to which this item belongs, and which active clusters we're definitely done with
        final List<Integer> clusterIdsToProcess = new ArrayList<>();
        // Clusters to which we simply add the item
        final List<Integer> clustersToAugment = new ArrayList<>();
        // New clusters, formed from subsets of currently active clusters, to which we will add the item
        final Set<List<Integer>> clustersToSeedWith = new HashSet<>();    // Use set to prevent creating duplicate clusters
        for (final Integer clusterIndex : getClusterIds()) {
            final Cluster cluster = getCluster(clusterIndex);
            final List<Integer> clusterItems = cluster.getMembers();
            if (item.getPositionA() <= cluster.getMaxClusterableStart()) {
                if (clusteringType.equals(CLUSTERING_TYPE.MAX_CLIQUE)) {
                    final List<Integer> linkedClusterItems = clusterItems.stream().filter(linkedItems::contains).collect(Collectors.toList());
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
            final List<Set<Integer>> augmentedClusterItemLists = clustersToAugment.stream()
                    .map(this::getCluster)
                    .map(Cluster::getMembers)
                    .map(HashSet::new)
                    .collect(Collectors.toList());
            final List<Set<Integer>> triggeredClusterItemSets = new ArrayList<>(clustersToSeedWith.size() + augmentedClusterItemLists.size());
            // New clusters formed from subsets of the currently active clusters
            triggeredClusterItemSets.addAll(clustersToSeedWith.stream().map(HashSet::new).collect(Collectors.toList()));
            triggeredClusterItemSets.addAll(augmentedClusterItemLists);
            triggeredClusterItemSets.sort(Comparator.comparingInt(Set::size));
            for (int i = 0; i < triggeredClusterItemSets.size(); i++) {
                final Set<Integer> seedItems = triggeredClusterItemSets.get(i);
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

    /**
     * Add the item specified by {@param itemId} to the cluster specified by {@param clusterIndex}
     * and expand the clustering interval
     * @param clusterId
     * @param itemId
     */
    private final void addToCluster(final int clusterId, final Integer itemId) {
        final Cluster cluster = getCluster(clusterId);
        final List<Integer> clusterItems = cluster.getMembers();
        clusterItems.add(itemId);
        final T item = getItem(itemId);
        final int itemClusterableStartPosition = linkage.getMaxClusterableStartingPosition(item);
        cluster.setMaxClusterableStart(Math.max(cluster.getMaxClusterableStart(), itemClusterableStartPosition));
    }

    /**
     * Creates a new cluster by agglomerating clusters with the given ids together (and deleting them from the currently
     * active set), along with an additional item.
     * @param clusterIds ids of clusters to combine
     * @param itemId id of item to add to new cluster
     */
    protected final void combineClusters(final Collection<Integer> clusterIds, final Integer itemId) {
        final Collection<Cluster> clusters = removeClusters(clusterIds);
        final List<Integer> clusterItems = clusters.stream()
                .map(Cluster::getMembers)
                .flatMap(List::stream)
                .distinct()
                .collect(Collectors.toList());
        final List<Integer> newClusterItems = new ArrayList<>(clusterItems.size() + 1);
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
    protected final void seedWithExistingCluster(final Integer item, final Collection<Integer> seedItems) {
        final List<Integer> newClusterItems = new ArrayList<>(1 + seedItems.size());
        newClusterItems.addAll(seedItems);
        newClusterItems.add(item);
        registerCluster(createCluster(newClusterItems));
    }
}
