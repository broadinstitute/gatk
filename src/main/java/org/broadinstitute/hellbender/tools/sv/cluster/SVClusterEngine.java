package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>Base class for clustering items that possess start/end genomic coordinates. Efficient algorithms are implemented for
 * single-linkage and max-clique clustering that leverage the low dimensionality of position-based clustering criteria.
 * These algorithms are suitable for clustering problems testing for overlapping events in coordinate-sorted order, with
 * additional possible criteria, when the maximum feasible starting position for an item can be easily estimated (e.g.
 * reciprocal overlap, end-point distance).</p>
 *
 * <p>Based on the method from Marschall T, Costa VG, Canzar S, et al. CLEVER: Clique-enumerating variant finder.
 * Bioinformatics. 2012;28(22):2875-2882.</p>
 *
 * <p>NOTE: precise implementation of {@link SVClusterLinkage#getMaxClusterableStartingPosition(SVLocatable)}
 * is important for efficiency because it determines when a cluster can be finalized and omitted from further clustering tests.</p>
 *
 * @param <T> class of items to cluster
 */
public class SVClusterEngine<T extends SVLocatable> {

    /**
     * Available clustering algorithms
     */
    public enum CLUSTERING_TYPE {
        SINGLE_LINKAGE,
        MAX_CLIQUE
    }

    private final SVCollapser<T> collapser; // Flattens clusters into a single representative item for output
    private final SVClusterLinkage<T> linkage;
    private Map<Integer, Cluster> idToClusterMap; // Active clusters
    private final Map<Integer, T> idToItemMap; // Active items
    private final List<T> outputBuffer;
    protected final CLUSTERING_TYPE clusteringType;
    private String currentContig;
    private int nextItemId;
    private int nextClusterId;
    private int lastStart;

    /**
     * @param clusteringType algorithm choice
     * @param collapser function that ingests a collection of clustered items and returns a single representative item
     */
    public SVClusterEngine(final CLUSTERING_TYPE clusteringType,
                           final SVCollapser collapser,
                           final SVClusterLinkage<T> linkage) {
        this.clusteringType = clusteringType;
        this.collapser = Utils.nonNull(collapser);
        this.linkage = Utils.nonNull(linkage);
        idToClusterMap = new HashMap<>();
        outputBuffer = new ArrayList<>();
        currentContig = null;
        idToItemMap = new HashMap<>();
        nextItemId = 0;
        nextClusterId = 0;
        lastStart = 0;
    }


    /**
     * Flushes all active clusters, adding them to the output buffer. Results from the output buffer are then copied out
     * and the buffer is cleared. This should be called between contigs to save memory.
     */
    public final List<T> getOutput() {
        flushClusters();
        final List<T> output = new ArrayList<>(outputBuffer);
        outputBuffer.clear();
        return output;
    }

    @VisibleForTesting
    public SVCollapser<T> getCollapser() {
        return collapser;
    }

    @VisibleForTesting
    public SVClusterLinkage<T> getLinkage() {
        return linkage;
    }

    /**
     * Returns true if there are any active or finalized clusters.
     */
    public final boolean isEmpty() {
        return idToClusterMap.isEmpty() && outputBuffer.isEmpty();
    }

    /**
     * Adds and clusters the given item. Note that items must be added in order of increasing start position.
     * @param item item to cluster
     */
    public final void add(final T item) {
        // Start a new cluster if on a new contig
        if (!item.getContigA().equals(currentContig)) {
            flushClusters();
            currentContig = item.getContigA();
            lastStart = 0;
            seedCluster(registerItem(item));
            return;
        }
        final int itemId = registerItem(item);
        final List<Integer> clusterIdsToProcess = cluster(itemId);
        processClusters(clusterIdsToProcess);
    }

    private final int registerItem(final T item) {
        Utils.validate(item.getPositionA() >= lastStart, "Items must be added in order of increasing start coordinate");
        lastStart = item.getPositionA();
        final int itemId = nextItemId++;
        idToItemMap.put(itemId, item);
        return itemId;
    }

    private final int getMaxClusterableStartingPositionByIds(final Collection<Integer> itemIds) {
        Utils.nonNull(itemIds);
        Utils.nonEmpty(itemIds);
        return linkage.getMaxClusterableStartingPosition(itemIds.stream().map(this::getItem).collect(Collectors.toList()));
    }

    /**
     * Add a new {@param <T>} to the current clusters and determine which are complete
     * @param itemId id of registered item to add
     * @return the IDs for clusters that are complete and ready for processing
     */
    private final List<Integer> cluster(final Integer itemId) {
        final T item = getItem(itemId);
        // Get list of item IDs from active clusters that cluster with this item
        final Set<Integer> linkedItems = idToClusterMap.values().stream().map(Cluster::getItemIds)
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
        for (final Map.Entry<Integer, Cluster> entry : idToClusterMap.entrySet()) {
            final Integer clusterIndex = entry.getKey();
            final Cluster cluster = entry.getValue();
            final List<Integer> clusterItems = cluster.getItemIds();
            if (item.getPositionA() > cluster.getMaxClusterableStart()) {
                clusterIdsToProcess.add(clusterIndex);  //this cluster is complete -- process it when we're done
            } else {
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
                    .map(idToClusterMap::get)
                    .map(Cluster::getItemIds)
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
        return clusterIdsToProcess;
    }

    /**
     * Creates a new cluster by agglomerating clusters with the given ids together (and deleting them from the currently
     * active set), along with an additional item.
     * @param clusterIds ids of clusters to combine
     * @param itemId id of item to add to new cluster
     */
    private final void combineClusters(final Collection<Integer> clusterIds, final Integer itemId) {
        final List<Cluster> clusters = clusterIds.stream().map(this::getCluster).collect(Collectors.toList());
        clusterIds.stream().forEach(idToClusterMap::remove);
        final List<Integer> clusterItems = clusters.stream()
                .map(Cluster::getItemIds)
                .flatMap(List::stream)
                .distinct()
                .collect(Collectors.toList());
        final List<Integer> newClusterItems = new ArrayList<>(clusterItems.size() + 1);
        newClusterItems.addAll(clusterItems);
        newClusterItems.add(itemId);
        idToClusterMap.put(nextClusterId++, new Cluster(getMaxClusterableStartingPositionByIds(newClusterItems), newClusterItems));
    }

    /**
     * Finalizes a single cluster, removing it from the currently active set and adding it to the output buffer.
     */
    private final void processCluster(final int clusterIndex) {
        final Cluster cluster = getCluster(clusterIndex);
        idToClusterMap.remove(clusterIndex);
        outputBuffer.add(collapser.collapse(cluster.getItemIds().stream().map(idToItemMap::get).collect(Collectors.toList())));
    }

    /**
     * Finalizes a set of clusters.
     */
    private final void processClusters(final List<Integer> clusterIdsToProcess) {
        for (final Integer clusterId : clusterIdsToProcess) {
            processCluster(clusterId);
        }
    }

    /**
     * Finalizes all active clusters and adds them to the output buffer. Also clears the currently active set of clusters
     * and items.
     */
    private final void flushClusters() {
        final List<Integer> clustersToFlush = new ArrayList<>(idToClusterMap.keySet());
        for (final Integer clusterId : clustersToFlush) {
            processCluster(clusterId);
        }
        idToItemMap.clear();
        nextItemId = 0;
        nextClusterId = 0;
    }

    /**
     * Creates a new singleton cluster containing the given item.
     */
    private final void seedCluster(final Integer item) {
        final List<Integer> newClusters = new ArrayList<>(1);
        newClusters.add(item);
        idToClusterMap.put(nextClusterId++, new Cluster(linkage.getMaxClusterableStartingPosition(getItem(item)), newClusters));
    }

    /**
     * Create a new cluster containing a new item and a set of existing items, and add it to the set of active clusters.
     * Note that if there exists an identical activate cluster, it will not be added.
     * @param item new item (assumed registered)
     * @param seedItems existing items
     */
    private final void seedWithExistingCluster(final Integer item, final Collection<Integer> seedItems) {
        final List<Integer> newClusterItems = new ArrayList<>(1 + seedItems.size());
        newClusterItems.addAll(seedItems);
        newClusterItems.add(item);
        idToClusterMap.put(nextClusterId++,
                new Cluster(getMaxClusterableStartingPositionByIds(newClusterItems), newClusterItems));
    }

    private final Cluster getCluster(final int id) {
        Utils.validateArg(idToClusterMap.containsKey(id), "Cluster ID " + id + " does not exist.");
        return idToClusterMap.get(id);
    }

    private final T getItem(final int id) {
        Utils.validateArg(idToItemMap.containsKey(id), "Item ID " + id + " does not exist.");
        return idToItemMap.get(id);
    }

    /**
     * Add the item specified by {@param itemId} to the cluster specified by {@param clusterIndex}
     * and expand the clustering interval
     * @param clusterId
     * @param itemId
     */
    private final void addToCluster(final int clusterId, final Integer itemId) {
        final Cluster cluster = getCluster(clusterId);
        final List<Integer> clusterItems = cluster.getItemIds();
        clusterItems.add(itemId);
        final T item = getItem(itemId);
        final int itemClusterableStartPosition = linkage.getMaxClusterableStartingPosition(item);
        cluster.setMaxClusterableStart(Math.max(cluster.getMaxClusterableStart(), itemClusterableStartPosition));
    }

    /**
     * Container class for clustered items
     */
    private static final class Cluster {
        private int maxClusterableStart;
        private final List<Integer> itemIds;

        public Cluster(final int maxClusterableStart, final List<Integer> itemIds) {
            Utils.nonNull(itemIds);
            this.maxClusterableStart = maxClusterableStart;
            this.itemIds = itemIds;
        }

        public int getMaxClusterableStart() {
            return maxClusterableStart;
        }

        public void setMaxClusterableStart(final int position) {
            maxClusterableStart = position;
        }

        public List<Integer> getItemIds() {
            return itemIds;
        }

        /**
         * Note we do not check for equality on max clusterable start position, which could be dependent on the
         * state of the engine.
         */
        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Cluster)) return false;
            Cluster cluster = (Cluster) o;
            return itemIds.equals(cluster.itemIds);
        }

        @Override
        public int hashCode() {
            return Objects.hash(itemIds);
        }
    }
}
