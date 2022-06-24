package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
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
    // TODO move to outer, must sort on collapsed items
    protected final Comparator<SVLocatable> itemComparator;
    protected final Comparator<Locatable> intervalComparator;

    private Map<Integer, C> idToClusterMap; // Active clusters
    private final Map<Long, T> idToItemMap; // Active items
    private String currentContig;
    private int nextClusterId;
    private T lastItem;

    public SVClusterEngine(final SVClusterLinkage<T> linkage,
                           final SAMSequenceDictionary dictionary) {
        this.linkage = Utils.nonNull(linkage);
        itemComparator = SVCallRecordUtils.getSVLocatableComparator(dictionary);
        intervalComparator = IntervalUtils.getDictionaryOrderComparator(dictionary);
        idToItemMap = new HashMap<>();
        idToClusterMap = new HashMap<>();
        currentContig = null;
        nextClusterId = 0;
        lastItem = null;
    }

    protected abstract void cluster(Long itemId);
    protected abstract R createOutputCluster(C cluster);

    protected final int getMaxClusterableStartingPositionByIds(final Collection<Long> itemIds) {
        Utils.nonNull(itemIds);
        Utils.nonEmpty(itemIds);
        return linkage.getMaxClusterableStartingPosition(itemIds.stream().map(this::getItem).collect(Collectors.toList()));
    }

    /**
     * Gets any available finalized clusters.
     */
    public final List<R> flush(final boolean hard) {
        return hard ? hardFlush() : softFlush();
    }

    private boolean intervalMayClusterWith(final SimpleInterval interval, final C cluster) {
        final int maxStart = cluster.getMaxClusterableStart();
        return intervalComparator.compare(interval, new SimpleInterval(cluster.getContig(), maxStart, maxStart)) <= 0;
    }

    private boolean intervalMayClusterWith(final SimpleInterval interval, final T item) {
        final int maxStart = linkage.getMaxClusterableStartingPosition(item);
        return intervalComparator.compare(interval, new SimpleInterval(item.getContigA(), maxStart, maxStart)) <= 0;
    }

    private final List<R> softFlush() {
        if (lastItem == null) {
            return Collections.emptyList();
        }
        final SimpleInterval lastStart = new SimpleInterval(lastItem.getContigA(), lastItem.getPositionA(), lastItem.getPositionA());
        final List<Map.Entry<Integer, C>> finalizedClusters = getClusterIds().stream()
                .map(id -> new AbstractMap.SimpleEntry<>(id, getCluster(id)))
                .filter(e -> !intervalMayClusterWith(lastStart, e.getValue()))
                .collect(Collectors.toList());
        final List<R> output = finalizedClusters.stream().map(Map.Entry::getValue).map(this::createOutputCluster).collect(Collectors.toList());
        finalizedClusters.stream().map(Map.Entry::getKey).forEach(idToClusterMap::remove);
        final Set<Long> outputItemIds = output.stream().map(R::getAllIds).flatMap(Collection::stream).collect(Collectors.toSet());
        final Set<Long> remainingClusterItemIds = idToClusterMap.values().stream().map(C::getAllIds).flatMap(Collection::stream).collect(Collectors.toSet());
        getItemIds().stream().filter(id -> outputItemIds.contains(id) && !remainingClusterItemIds.contains(id)).forEach(this::removeItem);
        // Remove items not in clusters but also no longer clusterable with any new items
        getItemIds().stream().filter(id -> !intervalMayClusterWith(lastStart, getItem(id)) && !remainingClusterItemIds.contains(id)).forEach(this::removeItem);
        return output;
    }

    protected void removeItem(final Long id) {
        idToItemMap.remove(id);
    }

    protected List<R> hardFlush() {
        final List<R> output = getClusters().stream().map(this::createOutputCluster).collect(Collectors.toList());
        idToClusterMap.clear();
        idToItemMap.clear();
        currentContig = null;
        lastItem = null;
        return output;
    }

    @VisibleForTesting
    public SVClusterLinkage<T> getLinkage() {
        return linkage;
    }

    /**
     * Returns true if there are any active or finalized clusters.
     */
    public final boolean isEmpty() {
        return idToClusterMap.isEmpty();
    }

    /**
     * Adds and clusters the given item. Note that items must be added in order of increasing start position.
     * @param item item to cluster
     */
    public final void add(final T item, final Long id) {
        Utils.nonNull(id);
        Utils.validateArg(!idToItemMap.containsKey(id), "Item with id " + id + " already exists");
        registerItem(item, id);
        cluster(id);
    }

    /**
     * May override
     * @param itemId
     * @param item
     */
    protected void putNewItem(final T item, final Long itemId) {
        idToItemMap.put(itemId, item);
    }

    private final Long registerItem(final T item, final Long id) {
        // Start a new cluster if on a new contig
        if (!item.getContigA().equals(currentContig)) {
            currentContig = item.getContigA();
        } else {
            Utils.validate(item.getPositionA() >= lastItem.getPositionA(), "Items must be added in order of increasing start coordinate");
        }
        putNewItem(item, id);
        lastItem = item;
        return id;
    }

    protected final Collection<C> removeClusters(final Collection<Integer> clusterIds) {
        final Collection<C> clusters = clusterIds.stream().map(this::getCluster).collect(Collectors.toList());
        clusterIds.stream().forEach(idToClusterMap::remove);
        return clusters;
    }

    protected final Set<Integer> getClusterIds() {
        return idToClusterMap.keySet();
    }

    protected final Collection<C> getClusters() {
        return idToClusterMap.values();
    }

    protected final C getCluster(final Integer id) {
        Utils.validateArg(idToClusterMap.containsKey(id), "Cluster ID " + id + " does not exist.");
        return idToClusterMap.get(id);
    }

    protected T getItem(final Long id) {
        Utils.validateArg(idToItemMap.containsKey(id), "Item ID " + id + " does not exist.");
        return idToItemMap.get(id);
    }

    protected final String getCurrentContig() {
        return currentContig;
    }

    public Set<Long> getItemIds() {
        return idToItemMap.keySet();
    }

    protected Integer registerCluster(final C cluster) {
        final Integer id = nextClusterId++;
        idToClusterMap.put(id, cluster);
        return id;
    }

}
