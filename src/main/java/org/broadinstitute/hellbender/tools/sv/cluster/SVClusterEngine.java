package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
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
public abstract class SVClusterEngine<T extends SVLocatable, C extends BasicCluster, R extends BasicOutputCluster<T>> {

    /**
     * Available clustering algorithms
     */
    public enum CLUSTERING_TYPE {
        SINGLE_LINKAGE,
        MAX_CLIQUE
    }

    protected final SVClusterLinkage<T> linkage;
    protected final ClusterSortingBuffer<R> buffer;
    protected final Comparator<T> itemComparator;

    private Map<Integer, C> idToClusterMap; // Active clusters
    private final Map<Integer, T> idToItemMap; // Active items
    private String currentContig;
    private int nextItemId;
    private int nextClusterId;
    private int lastStart;
    private Integer minActiveStartingPositionItemId;

    public SVClusterEngine(final SVClusterLinkage<T> linkage,
                           final SAMSequenceDictionary dictionary) {
        this.linkage = Utils.nonNull(linkage);
        currentContig = null;
        idToItemMap = new HashMap<>();
        itemComparator = SVCallRecordUtils.getSVLocatableComparator(dictionary);
        nextItemId = 0;
        nextClusterId = 0;
        lastStart = 0;
        minActiveStartingPositionItemId = null;
        idToClusterMap = new HashMap<>();
        buffer = createClusterSortingBuffer();
    }

    protected abstract ClusterSortingBuffer<R> createClusterSortingBuffer();
    protected abstract void cluster(Integer itemId);
    protected abstract R createOutputCluster(final C cluster);
    protected abstract R createSingleItemOutputClusterForBuffer(final T item);

    protected final int getMaxClusterableStartingPositionByIds(final Collection<Integer> itemIds) {
        Utils.nonNull(itemIds);
        Utils.nonEmpty(itemIds);
        return linkage.getMaxClusterableStartingPosition(itemIds.stream().map(this::getItem).collect(Collectors.toList()));
    }

    protected final int getMinStartingPositionByIds(final Collection<Integer> itemIds) {
        Utils.nonNull(itemIds);
        Utils.nonEmpty(itemIds);
        return itemIds.stream().map(this::getItem).mapToInt(T::getPositionA).min().getAsInt();
    }


    protected boolean itemIsLessThanMinActiveStart(T item) {
        return item.getPositionA() < getMinActiveStartingPositionItem().getPositionA();
    }

    /**
     * Flushes all active clusters, adding them to the output buffer. Results from the output buffer are then copied out
     * and the buffer is cleared. This should be called between contigs to save memory.
     */
    public final List<R> forceFlush() {
        flushClusters();
        return buffer.forceFlush();
    }


    /**
     * Gets any available finalized clusters.
     */
    public final List<R> flush() {
        if (minActiveStartingPositionItemId == null) {
            return Collections.emptyList();
        }
        final T minActiveItem = getItem(minActiveStartingPositionItemId);
        final R minActiveComparisonCluster = createSingleItemOutputClusterForBuffer(minActiveItem);
        return buffer.flush(minActiveComparisonCluster);
    }

    @VisibleForTesting
    public SVClusterLinkage<T> getLinkage() {
        return linkage;
    }

    public T getMinActiveStartingPositionItem() {
        Utils.validate(minActiveStartingPositionItemId == null || idToItemMap.containsKey(minActiveStartingPositionItemId),
                "Unregistered item id " + minActiveStartingPositionItemId);
        return idToItemMap.get(minActiveStartingPositionItemId);
    }

    /**
     * Returns true if there are any active or finalized clusters.
     */
    public final boolean isEmpty() {
        return idToClusterMap.isEmpty() && buffer.isEmpty();
    }

    /**
     * Adds and clusters the given item. Note that items must be added in order of increasing start position.
     * @param item item to cluster
     */
    public final void add(final T item) {
        final int itemId = registerItem(item);
        cluster(itemId);
        finalizeClusters(item);
    }

    /**
     * May override
     * @param itemId
     * @param item
     */
    protected void putNewItem(final Integer itemId, final T item) {
        idToItemMap.put(itemId, item);
    }

    private final int registerItem(final T item) {
        // Start a new cluster if on a new contig
        if (!item.getContigA().equals(currentContig)) {
            flushClusters();
            currentContig = item.getContigA();
            lastStart = 0;
        } else {
            Utils.validate(item.getPositionA() >= lastStart, "Items must be added in order of increasing start coordinate");
        }
        lastStart = item.getPositionA();
        final int itemId = nextItemId++;
        putNewItem(itemId, item);
        if (minActiveStartingPositionItemId == null || itemIsLessThanMinActiveStart(item)) {
            minActiveStartingPositionItemId = itemId;
        }
        return itemId;
    }

    protected final Collection<C> removeClusters(final Collection<Integer> clusterIds) {
        final Collection<C> clusters = clusterIds.stream().map(this::getCluster).collect(Collectors.toList());
        clusterIds.stream().forEach(idToClusterMap::remove);
        return clusters;
    }

    /**
     * Finalizes a single cluster, removing it from the currently active set and adding it to the output buffer.
     * returns finalized item IDs
     */
    private final Collection<Integer> finalizeCluster(final int clusterIndex) {
        final C cluster = getCluster(clusterIndex);
        idToClusterMap.remove(clusterIndex);
        final Set<Integer> clusterItemIds = cluster.getMembers();
        buffer.add(createOutputCluster(cluster));
        // Clean up item id map
        if (clusterItemIds.size() == 1) {
            // Singletons won't be present in any other clusters
            idToItemMap.remove(clusterItemIds.iterator().next());
            return clusterItemIds;
        } else {
            // Need to check that items aren't present in any other clusters
            final Set<Integer> activeItemIds = idToClusterMap.values().stream()
                    .map(BasicCluster::getMembers)
                    .flatMap(Set::stream)
                    .collect(Collectors.toSet());
            final List<Integer> itemsToRemove = idToItemMap.keySet().stream().filter(i -> !activeItemIds.contains(i))
                    .collect(Collectors.toList());
            for (final Integer i : itemsToRemove) {
                idToItemMap.remove(i);
            }
            return itemsToRemove;
        }
    }

    /**
     * Scans active items for the current min active starting position.
     */
    private final void findAndSetMinActiveStart() {
        minActiveStartingPositionItemId = null;
        T minActiveStartingPositionItem = null;
        for (final Integer itemId : idToItemMap.keySet()) {
            final T item = idToItemMap.get(itemId);
            if (minActiveStartingPositionItemId == null || itemComparator.compare(item, minActiveStartingPositionItem) < 0) {
                minActiveStartingPositionItemId = itemId;
                minActiveStartingPositionItem = idToItemMap.get(itemId);
            }
        }
    }

    /**
     * Finalizes clusters given the current item. If passed null, finalizes all clusters.
     */
    protected final void finalizeClusters(final T currentItem) {
        final Set<Integer> finalizedClusterIds;
        if (currentItem == null) {
            // Copy ids to avoid ConcurrentModificationException
            finalizedClusterIds = new HashSet<>(idToClusterMap.keySet());
        } else {
            finalizedClusterIds = idToClusterMap.entrySet().stream()
                    .filter(e -> currentItem.getPositionA() > e.getValue().getMaxClusterableStart())
                    .map(Map.Entry::getKey)
                    .collect(Collectors.toSet());
        }
        final Set<Integer> finalizedItemIds = finalizedClusterIds.stream().map(this::finalizeCluster)
                .flatMap(Collection::stream).collect(Collectors.toSet());
        // Update min active start position
        if (finalizedItemIds.contains(minActiveStartingPositionItemId)) {
            findAndSetMinActiveStart();
        }
    }

    /**
     * Finalizes all active clusters and adds them to the output buffer. Also clears the currently active set of clusters
     * and items.
     */
    protected final void flushClusters() {
        finalizeClusters(null);
        idToItemMap.clear();
        minActiveStartingPositionItemId = null;
        nextItemId = 0;
        nextClusterId = 0;
    }

    protected final Set<Integer> getClusterIds() {
        return idToClusterMap.keySet();
    }

    protected final Collection<C> getClusters() {
        return idToClusterMap.values();
    }

    protected final C getCluster(final int id) {
        Utils.validateArg(idToClusterMap.containsKey(id), "Cluster ID " + id + " does not exist.");
        return idToClusterMap.get(id);
    }

    protected final T getItem(final int id) {
        Utils.validateArg(idToItemMap.containsKey(id), "Item ID " + id + " does not exist.");
        return idToItemMap.get(id);
    }

    protected final Set<Integer> getItemIds() {
        return idToItemMap.keySet();
    }

    protected final void registerCluster(final C cluster) {
        idToClusterMap.put(nextClusterId++, cluster);
    }

}
