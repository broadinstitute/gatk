package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class CrossRefSVClusterEngine {

    protected final Map<Long, SVCallRecord> refIdToItemMap;
    protected final Map<Long, CrossRefCluster> idToClusterMap;
    private final SVClusterLinkage<SVCallRecord> linkage;
    private final Comparator<SVLocatable> locatableComparator;
    private final Function<CrossRefOutputCluster, SVCallRecord> collapser;

    private final PriorityQueue<SVCallRecord> outputBuffer;

    private Integer lastItemStart;
    private String lastItemContig;

    public CrossRefSVClusterEngine(final SVClusterLinkage<SVCallRecord> linkage,
                                   final Function<CrossRefOutputCluster, SVCallRecord> collapser,
                                   final SAMSequenceDictionary dictionary) {
        this.linkage = Utils.nonNull(linkage);
        this.collapser = Utils.nonNull(collapser);
        this.locatableComparator = SVCallRecordUtils.getSVLocatableComparator(dictionary);
        outputBuffer = new PriorityQueue<>(SVCallRecordUtils.getCallComparator(dictionary));
        refIdToItemMap = new HashMap<>();
        idToClusterMap = new HashMap<>();
        lastItemStart = null;
        lastItemContig = null;
    }

    public List<SVCallRecord> flush(final boolean force) {
        final List<SVCallRecord> collapsedRecords = flushClusters(force).stream()
                .map(c -> new CrossRefOutputCluster(c.getTestItem(), c.getMembers().values()))
                .map(collapser)
                .collect(Collectors.toList());
        outputBuffer.addAll(collapsedRecords);
        return flushBuffer(force);
    }

    private List<SVCallRecord> flushBuffer(final boolean force) {
        if (force) {
            final List<SVCallRecord> output = new ArrayList<>(outputBuffer.size());
            while (!outputBuffer.isEmpty()) {
                output.add(outputBuffer.poll());
            }
            return output;
        } else {
            final ArrayList<SVCallRecord> output = new ArrayList<>();
            final Integer minActiveStartPosition = minActiveStartPosition();
            while (!outputBuffer.isEmpty() &&
                    (minActiveStartPosition == null || outputBuffer.peek().getPositionA() <= minActiveStartPosition)) {
                output.add(outputBuffer.poll());
            }
            output.trimToSize();
            return output;
        }
    }

    private Integer minActiveStartPosition() {
        return idToClusterMap.isEmpty() ? null : idToClusterMap.values().stream().mapToInt(c -> c.getTestItem().getPositionA()).min().getAsInt();
    }

    private List<CrossRefCluster> flushClusters(final boolean force) {
        if (force) {
            final List<CrossRefCluster> output = new ArrayList<>(idToClusterMap.values());
            refIdToItemMap.clear();
            idToClusterMap.clear();
            lastItemStart = null;
            lastItemContig = null;
            return output;
        } else {
            // Find finalized ref items
            final List<Long> finalizedRefItems = refIdToItemMap.entrySet().stream()
                    .filter(e -> linkage.getMaxClusterableStartingPosition(e.getValue()) < lastItemStart)
                    .map(Map.Entry::getKey)
                    .collect(Collectors.toList());
            finalizedRefItems.forEach(refIdToItemMap::remove);
            // Find finalized clusters
            final List<Map.Entry<Long, CrossRefCluster>> finalizedClusters = idToClusterMap.entrySet().stream()
                    .filter(e -> e.getValue().getMaxClusterableStartingPosition() < lastItemStart)
                    .collect(Collectors.toList());
            finalizedClusters.forEach(e -> idToClusterMap.remove(e.getKey()));
            return finalizedClusters.stream().map(Map.Entry::getValue).collect(Collectors.toList());
        }
    }

    public void add(final Long itemId, final SVCallRecord item, final boolean isRefVariant) {
        Utils.validateArg(!refIdToItemMap.containsKey(itemId) && !idToClusterMap.containsKey(itemId), "Item id " + itemId + " already in use");
        Utils.validateArg(lastItemContig == null || lastItemContig.equals(item.getContigA()), "Attempted to add item on a new contig; please run a force flush beforehand");
        Utils.validateArg(lastItemStart == null || lastItemStart <= item.getPositionA(), "Items must be added in dictionary-sorted order");
        lastItemContig = item.getContigA();
        lastItemStart = item.getPositionA();
        if (isRefVariant) {
            idToClusterMap.values().stream()
                .filter(other -> linkage.areClusterable(other.getTestItem(), item))
                .forEach(cluster -> cluster.addMember(itemId, item));
        } else {
                final List<Map.Entry<Long, SVCallRecord>> linkedItems = refIdToItemMap.entrySet().stream()
                        .filter(other -> linkage.areClusterable(item, other.getValue()))
                        .collect(Collectors.toList());
                final int itemMaxClusterableStart = linkage.getMaxClusterableStartingPosition(item);
                final int maxPos;
                if (linkedItems.isEmpty()) {
                    maxPos = itemMaxClusterableStart;
                } else {
                    final int currentMaxPos = linkage.getMaxClusterableStartingPosition(linkedItems.stream().map(Map.Entry::getValue).collect(Collectors.toList()));
                    maxPos = Math.max(itemMaxClusterableStart, currentMaxPos);
                }
                final CrossRefCluster cluster = new CrossRefCluster(linkedItems.stream().collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)), itemId, item, maxPos);
                idToClusterMap.put(itemId, cluster);
        }
    }

    public class CrossRefOutputCluster {
        final SVCallRecord testItem;
        final Collection<SVCallRecord> items;
        public CrossRefOutputCluster(final SVCallRecord testItem, final Collection<SVCallRecord> items) {
            this.testItem = testItem;
            this.items = items;
        }

        public SVCallRecord getTestItem() {
            return testItem;
        }

        public Collection<SVCallRecord> getItems() {
            return items;
        }
    }

    public class CrossRefCluster {

        final Long testItemId;
        final SVCallRecord testItem;
        final Map<Long, SVCallRecord> members;
        final int maxClusterableStartingPosition;

        public CrossRefCluster(final Map<Long, SVCallRecord> members, final Long testItemId, final SVCallRecord testItem, final int maxClusterableStartingPosition) {
            this.members = Utils.nonNull(members);
            this.testItemId = Utils.nonNull(testItemId);
            this.testItem = Utils.nonNull(testItem);
            this.maxClusterableStartingPosition = maxClusterableStartingPosition;
        }

        public void addMember(final Long id, final SVCallRecord item) {
            Utils.validateArg(!members.containsKey(id), "Cluster already contains id " + id);
            members.put(id, item);
        }

        public Long getTestItemId() {
            return testItemId;
        }

        public SVCallRecord getTestItem() {
            return testItem;
        }

        public Map<Long, SVCallRecord> getMembers() {
            return members;
        }

        public int getMaxClusterableStartingPosition() {
            return maxClusterableStartingPosition;
        }

        public String getContig() {
            return testItem.getContigA();
        }
    }
}
