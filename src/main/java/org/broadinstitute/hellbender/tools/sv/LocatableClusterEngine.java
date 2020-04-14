package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class LocatableClusterEngine<T extends Locatable> {

    public enum CLUSTERING_TYPE {
        SINGLE_LINKAGE,
        MAX_CLIQUE
    }

    protected final SAMSequenceDictionary dictionary;
    private final List<Tuple2<SimpleInterval, List<Long>>> currentClusters; // Pairs of cluster start interval with item IDs
    private final Map<Long,T> idToItemMap;
    private final List<T> outputBuffer;
    private final CLUSTERING_TYPE clusteringType;
    private long currentItemId;
    private String currentContig;


    public LocatableClusterEngine(final SAMSequenceDictionary dictionary, final CLUSTERING_TYPE clusteringType) {
        this.dictionary = dictionary;
        this.clusteringType = clusteringType;
        this.currentClusters = new LinkedList<>();
        this.idToItemMap = new HashMap<>();
        this.outputBuffer = new ArrayList<>();
        currentItemId = 0;
        currentContig = null;
    }

    abstract protected boolean clusterTogether(final T a, final T b);
    abstract protected SimpleInterval getClusteringInterval(final T item, final SimpleInterval currentClusterInterval);
    abstract protected T deduplicateIdenticalItems(final Collection<T> items);
    abstract protected boolean itemsAreIdentical(final T a, final T b);
    abstract protected T flattenCluster(final Collection<T> cluster);

    public List<T> getOutput() {
        flushClusters();
        final List<T> output = deduplicateItems(outputBuffer);
        outputBuffer.clear();
        return output;
    }

    private void resetItemIds() {
        Utils.validate(currentClusters.isEmpty(), "Current cluster collection not empty");
        currentItemId = 0;
        idToItemMap.clear();
    }

    public boolean isEmpty() {
        return currentContig == null;
    }

    public void add(final T item) {

        // Start a new cluster if on a new contig
        if (!item.getContig().equals(currentContig)) {
            flushClusters();
            currentContig = item.getContig();
            idToItemMap.put(currentItemId, item);
            seedCluster(currentItemId);
            currentItemId++;
            return;
        }

        // Keep track of a unique id for each item
        idToItemMap.put(currentItemId, item);
        final List<Integer> clusterIdsToProcess = cluster(item);
        processFinalizedClusters(clusterIdsToProcess);
        deleteRedundantClusters();
        currentItemId++;
    }

    public String getCurrentContig() {
        return currentContig;
    }

    public List<T> deduplicateItems(final List<T> items) {
        final List<T> sortedItems = IntervalUtils.sortLocatablesBySequenceDictionary(items, dictionary);
        final List<T> deduplicatedList = new ArrayList<>();
        int i = 0;
        while (i < sortedItems.size()) {
            final T record = sortedItems.get(i);
            int j = i + 1;
            final Collection<Integer> identicalItemIndexes = new ArrayList<>();
            while (j < sortedItems.size() && record.getStart() == sortedItems.get(j).getStart()) {
                final T other = sortedItems.get(j);
                if (itemsAreIdentical(record, other)) {
                    identicalItemIndexes.add(j);
                }
                j++;
            }
            if (identicalItemIndexes.isEmpty()) {
                deduplicatedList.add(record);
                i++;
            } else {
                identicalItemIndexes.add(i);
                final List<T> identicalItems = identicalItemIndexes.stream().map(sortedItems::get).collect(Collectors.toList());
                deduplicatedList.add(deduplicateIdenticalItems(identicalItems));
                i = j;
            }
        }
        return deduplicatedList;
    }

    /**
     * Add a new {@param <T>} to the current clusters and determine which are complete
     * @param item to be added
     * @return the IDs for clusters that are complete and ready for processing
     */
    private List<Integer> cluster(final T item) {
        // Get list of item IDs from active clusters that cluster with this item
        final Set<Long> linkedItemIds = idToItemMap.entrySet().stream()
                .filter(other -> other.getKey().intValue() != currentItemId && clusterTogether(item, other.getValue()))
                .map(Map.Entry::getKey)
                .collect(Collectors.toCollection(LinkedHashSet::new));

        // Find clusters to which this item belongs, and which active clusters we're definitely done with
        int clusterIndex = 0;
        final List<Integer> clusterIdsToProcess = new ArrayList<>();
        final List<Integer> clustersToAdd = new ArrayList<>();
        final List<Integer> clustersToSeedWith = new ArrayList<>();
        for (final Tuple2<SimpleInterval, List<Long>> cluster : currentClusters) {
            final SimpleInterval clusterInterval = cluster._1;
            final List<Long> clusterItemIds = cluster._2;
            if (item.getStart() > clusterInterval.getEnd()) {
                clusterIdsToProcess.add(clusterIndex);  //this cluster is complete -- process it when we're done
            } else {
                if (clusteringType.equals(CLUSTERING_TYPE.SINGLE_LINKAGE)) {
                    final int n = (int) clusterItemIds.stream().filter(linkedItemIds::contains).count();
                    if (n == clusterItemIds.size()) {
                        clustersToAdd.add(clusterIndex);
                    } else if (n > 0) {
                        clustersToSeedWith.add(clusterIndex);
                    }
                } else if (clusteringType.equals(CLUSTERING_TYPE.MAX_CLIQUE)) {
                    final boolean matchesCluster = clusterItemIds.stream().anyMatch(linkedItemIds::contains);
                    if (matchesCluster) {
                        clustersToAdd.add(clusterIndex);
                    }
                } else {
                    throw new IllegalArgumentException("Clustering algorithm for type " + clusteringType.name() + " not implemented");
                }
            }
            clusterIndex++;
        }

        // Add to item clusters
        for (final int index : clustersToAdd) {
            addToCluster(index, currentItemId);
        }
        // Create new clusters/cliques
        for (final int index : clustersToSeedWith) {
            seedWithExistingCluster(currentItemId, index, linkedItemIds);
        }
        // If there weren't any matches, create a new singleton cluster
        if (clustersToAdd.isEmpty() && clustersToSeedWith.isEmpty()) {
            seedCluster(currentItemId);
        }
        return clusterIdsToProcess;
    }

    private void processCluster(final int clusterIndex) {
        final Tuple2<SimpleInterval, List<Long>> cluster = validateClusterIndex(clusterIndex);
        final List<Long> clusterItemIds = cluster._2;
        currentClusters.remove(clusterIndex);
        final List<T> clusterItems = clusterItemIds.stream().map(idToItemMap::get).collect(Collectors.toList());
        outputBuffer.add(flattenCluster(clusterItems));
    }

    private void processFinalizedClusters(final List<Integer> clusterIdsToProcess) {
        final Set<Integer> activeClusterIds = IntStream.range(0, currentClusters.size()).boxed().collect(Collectors.toSet());
        activeClusterIds.removeAll(clusterIdsToProcess);
        final Set<Long> activeClusterItemIds = activeClusterIds.stream().flatMap(i -> currentClusters.get(i)._2.stream()).collect(Collectors.toSet());
        final Set<Long> finalizedItemIds = clusterIdsToProcess.stream()
                .flatMap(i -> currentClusters.get(i)._2.stream())
                .filter(i -> !activeClusterItemIds.contains(i))
                .collect(Collectors.toSet());
        for (int i = clusterIdsToProcess.size() - 1; i >= 0; i--) {
            processCluster(clusterIdsToProcess.get(i));
        }
        finalizedItemIds.stream().forEach(idToItemMap::remove);
    }

    private void deleteRedundantClusters() {
        final Set<Integer> redundantClusterSet = new HashSet<>();
        for (int i = 0; i < currentClusters.size(); i++) {
            final Set<Long> clusterSetA = new HashSet<>(currentClusters.get(i)._2);
            for (int j = 0; j < i; j++) {
                final Set<Long> clusterSetB = new HashSet<>(currentClusters.get(j)._2);
                if (clusterSetA.containsAll(clusterSetB)) {
                    redundantClusterSet.add(j);
                } else if (clusterSetA.size() != clusterSetB.size() && clusterSetB.containsAll(clusterSetA)) {
                    redundantClusterSet.add(i);
                }
            }
        }
        final List<Integer> redundantClustersList = new ArrayList<>(redundantClusterSet);
        redundantClustersList.sort(Comparator.naturalOrder());
        for (int i = redundantClustersList.size() - 1; i >= 0; i--) {
            currentClusters.remove((int)redundantClustersList.get(i));
        }
    }

    private void flushClusters() {
        while (!currentClusters.isEmpty()) {
            processCluster(0);
        }
        resetItemIds();
    }

    private void seedCluster(final long seedId) {
        final T seed = validateItemIndex(seedId);
        final List<Long> newCluster = new ArrayList<>(1);
        newCluster.add(seedId);
        currentClusters.add(new Tuple2<>(getClusteringInterval(seed, null), newCluster));
    }

    /**
     * Create a new cluster
     * @param seedId    itemId
     * @param existingClusterIndex
     * @param clusteringIds
     */
    private void seedWithExistingCluster(final Long seedId, final int existingClusterIndex, final Set<Long> clusteringIds) {
        final T seed = validateItemIndex(seedId);
        final List<Long> existingCluster = currentClusters.get(existingClusterIndex)._2;
        final List<Long> validClusterIds = existingCluster.stream().filter(clusteringIds::contains).collect(Collectors.toList());
        final List<Long> newCluster = new ArrayList<>(1 + existingCluster.size());
        newCluster.addAll(validClusterIds);
        newCluster.add(seedId);
        currentClusters.add(new Tuple2<>(getClusteringInterval(seed, currentClusters.get(existingClusterIndex)._1), newCluster));
    }

    private T validateItemIndex(final long index) {
        final T item = idToItemMap.get(index);
        if (item == null) {
            throw new IllegalArgumentException("Item id " + index + " not found in table");
        }
        if (!currentContig.equals(item.getContig())) {
            throw new IllegalArgumentException("Attempted to seed new cluster with item on contig " + item.getContig() + " but the current contig is " + currentContig);
        }
        return item;
    }

    private Tuple2<SimpleInterval, List<Long>> validateClusterIndex(final int index) {
        if (index < 0 || index >= currentClusters.size()) {
            throw new IllegalArgumentException("Specified cluster index " + index + " is out of range.");
        }
        final Tuple2<SimpleInterval, List<Long>> cluster = currentClusters.get(index);
        final List<Long> clusterItemIds = cluster._2;
        if (clusterItemIds.isEmpty()) {
            throw new IllegalArgumentException("Encountered empty cluster");
        }
        return cluster;
    }

    /**
     * Add the item specified by {@param itemId} to the cluster specified by {@param clusterIndex}
     * and expand the clustering interval
     * @param clusterIndex
     * @param itemId
     */
    private void addToCluster(final int clusterIndex, final long itemId) {
        final T item = idToItemMap.get(itemId);
        if (item == null) {
            throw new IllegalArgumentException("Item id " + item + " not found in table");
        }
        if (!currentContig.equals(item.getContig())) {
            throw new IllegalArgumentException("Attempted to add new item on contig " + item.getContig() + " but the current contig is " + currentContig);
        }
        if (clusterIndex >= currentClusters.size()) {
            throw new IllegalArgumentException("Specified cluster index " + clusterIndex + " is greater than the largest index.");
        }
        final Tuple2<SimpleInterval, List<Long>> cluster = currentClusters.get(clusterIndex);
        final SimpleInterval clusterInterval = cluster._1;
        final List<Long> clusterItems = cluster._2;
        clusterItems.add(itemId);
        final SimpleInterval clusteringStartInterval = getClusteringInterval(item, clusterInterval);
        if (clusteringStartInterval.getStart() != clusterInterval.getStart() || clusteringStartInterval.getEnd() != clusterInterval.getEnd()) {
            currentClusters.remove(clusterIndex);
            currentClusters.add(clusterIndex, new Tuple2<>(clusteringStartInterval, clusterItems));
        }
    }
}
