package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.BoundType;
import com.google.common.collect.SortedMultiset;
import com.google.common.collect.TreeMultiset;
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
public abstract class SVClusterEngine<T extends SVLocatable, C extends SVClusterEngine.Cluster> {

    /**
     * Available clustering algorithms
     */
    public enum CLUSTERING_TYPE {
        SINGLE_LINKAGE,
        MAX_CLIQUE
    }

    protected final SVClusterLinkage<T> linkage;
    protected ClusterSortingBuffer buffer;
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
        createItemClusterComparator();
    }

    protected abstract void createItemClusterComparator();
    protected abstract void cluster(Integer itemId);
    protected abstract C createCluster(List<Integer> itemIds);


    protected boolean itemIsMinActiveStart(T item) {
        return item.getPositionA() < getMinActiveStartingPositionItem().getPositionA();
    }

    /**
     * Flushes all active clusters, adding them to the output buffer. Results from the output buffer are then copied out
     * and the buffer is cleared. This should be called between contigs to save memory.
     */
    public final List<OutputCluster> forceFlush() {
        flushClusters();
        return buffer.forceFlush();
    }

    /**
     * Gets any available finalized clusters.
     */
    public final List<OutputCluster> flush() {
        if (minActiveStartingPositionItemId == null) {
            return Collections.emptyList();
        }
        final T minActiveItem = idToItemMap.get(minActiveStartingPositionItemId);
        Utils.nonNull(minActiveItem, "Unknown item id");
        final OutputCluster minActiveComparisonCluster = new OutputCluster(Collections.singletonList(minActiveItem));
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
        putNewItem(nextItemId++, item);
        if (minActiveStartingPositionItemId == null || itemIsMinActiveStart(item)) {
            minActiveStartingPositionItemId = itemId;
        }
        return itemId;
    }

    protected final Collection<C> removeClusters(final Collection<Integer> clusterIds) {
        final Collection<C> clusters = clusterIds.stream().map(this::getCluster).collect(Collectors.toList());
        clusterIds.stream().forEach(idToClusterMap::remove);
        return clusters;
    }

    protected final int getMaxClusterableStartingPositionByIds(final Collection<Integer> itemIds) {
        Utils.nonNull(itemIds);
        Utils.nonEmpty(itemIds);
        return linkage.getMaxClusterableStartingPosition(itemIds.stream().map(this::getItem).collect(Collectors.toList()));
    }

    /**
     * Finalizes a single cluster, removing it from the currently active set and adding it to the output buffer.
     */
    private final void finalizeCluster(final int clusterIndex) {
        final Cluster idCluster = getCluster(clusterIndex);
        idToClusterMap.remove(clusterIndex);
        final List<Integer> clusterItemIds = idCluster.getMembers();
        final OutputCluster itemCluster = new OutputCluster(idCluster.getMembers().stream().map(idToItemMap::get).collect(Collectors.toList()));
        buffer.add(itemCluster);
        // Clean up item id map
        if (clusterItemIds.size() == 1) {
            // Singletons won't be present in any other clusters
            idToItemMap.remove(clusterItemIds.get(0));
        } else {
            // Need to check that items aren't present in any other clusters
            final Set<Integer> activeItemIds = idToClusterMap.values().stream()
                    .map(Cluster::getMembers)
                    .flatMap(List::stream)
                    .collect(Collectors.toSet());
            final List<Integer> itemsToRemove = idToItemMap.keySet().stream().filter(i -> !activeItemIds.contains(i))
                    .collect(Collectors.toList());
            for (final Integer i : itemsToRemove) {
                idToItemMap.remove(i);
            }
        }
        // Update min active start position
        if (clusterItemIds.contains(minActiveStartingPositionItemId)) {
            findAndSetMinActiveStart();
        }
    }

    /**
     * Scans active items for the current min active starting position.
     */
    private final void findAndSetMinActiveStart() {
        minActiveStartingPositionItemId = null;
        T minActiveStartingPositionItem = null;
        for (final Integer itemId : getClusterableItemIds()) {
            final T item = idToItemMap.get(itemId);
            if (minActiveStartingPositionItemId == null || itemComparator.compare(item, minActiveStartingPositionItem) < 0) {
                minActiveStartingPositionItemId = itemId;
                minActiveStartingPositionItem = idToItemMap.get(itemId);
            }
        }
    }

    /**
     * Finalizes a set of clusters.
     */
    protected final void finalizeClusters(final T currentItem) {
        for (final Map.Entry<Integer, C> entry : idToClusterMap.entrySet()) {
            if (currentItem == null || currentItem.getPositionA() > entry.getValue().getMaxClusterableStart()) {
                finalizeCluster(entry.getKey());
            }
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

    /**
     * All item ids that can be clustered. May override when Cluster subclasses include additional item groups.
     * @return
     */
    protected Set<Integer> getClusterableItemIds() {
        return getItemIds();
    }

    /**
     * Item ids that are members of clusters
     * @return
     */
    protected final Set<Integer> getItemIds() {
        return idToItemMap.keySet();
    }

    protected final void registerCluster(final C cluster) {
        idToClusterMap.put(nextClusterId++, cluster);
    }

    /**
     * Container class for clustered items
     */
    public static class Cluster {
        private int maxClusterableStart;
        private final List<Integer> members;

        public Cluster(final int maxClusterableStart, final List<Integer> members) {
            Utils.nonNull(members);
            this.maxClusterableStart = maxClusterableStart;
            this.members = new ArrayList<>(members);
        }

        public int getMaxClusterableStart() {
            return maxClusterableStart;
        }

        public void setMaxClusterableStart(final int position) {
            maxClusterableStart = position;
        }

        public List<Integer> getMembers() {
            return members;
        }

        public void addMember(final Integer member, final int memberMaxClusterableStart) {
            members.add(member);
            maxClusterableStart = Math.max(maxClusterableStart, memberMaxClusterableStart);
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
            return members.equals(cluster.members);
        }
    }

    public class OutputCluster {
        private final List<T> members;

        public OutputCluster(final List<T> members) {
            Utils.nonNull(members);
            this.members = members;
        }

        public List<T> getMembers() {
            return members;
        }
    }

    protected final class ClusterSortingBuffer {
        private SortedMultiset<OutputCluster> buffer;

        public ClusterSortingBuffer(final Comparator<OutputCluster> clusterComparator) {
            this.buffer = TreeMultiset.create(clusterComparator);
        }

        public void add(final OutputCluster cluster) {
            buffer.add(cluster);
        }

        /**
         * Returns any records that can be safely flushed based on the current minimum starting position
         * of items still being actively clustered.
         */
        public List<OutputCluster> flush(final OutputCluster minActiveStartingPositionCluster) {
            if (buffer.isEmpty()) {
                return Collections.emptyList();
            }
            if (minActiveStartingPositionCluster == null) {
                forceFlush();
            }
            final SortedMultiset<OutputCluster> finalizedView = buffer.headMultiset(minActiveStartingPositionCluster, BoundType.CLOSED);
            final ArrayList<OutputCluster> finalizedClusters = new ArrayList<>(finalizedView);
            // Clearing a view of the buffer also clears the items from the buffer itself
            finalizedView.clear();
            return finalizedClusters;
        }

        /**
         * Returns all buffered records, regardless of any active clusters. To be used only when certain that no
         * active clusters can be clustered with any future inputs.
         */
        public List<OutputCluster> forceFlush() {
            final List<OutputCluster> result = new ArrayList<>(buffer);
            buffer.clear();
            return result;
        }

        public boolean isEmpty() {
            return buffer.isEmpty();
        }
    }
}
