package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.tools.sv.SVLocus;

import java.util.*;
import java.util.stream.Collectors;

public final class SortingClusterCollapser<T extends SVLocatable, C extends BasicOutputCluster<T>> {

    private final Map<Long, SVLocus> activeItemLoci;
    private final PriorityQueue<Long> activeItemIds;
    private final SVCollapser<T, C> collapser;
    private final SVClusterEngine<T, ?, C> engine;
    private final Comparator<SVLocatable> locusComparator;
    private final PriorityQueue<T> outputBuffer;
    private long nextItemId;

    public SortingClusterCollapser(final SVClusterEngine<T, ?, C> engine,
                                   final SVCollapser<T, C> collapser,
                                   final SAMSequenceDictionary dictionary) {
        this.engine = engine;
        this.collapser = collapser;
        this.activeItemLoci = new HashMap<>();
        this.locusComparator = SVCallRecordUtils.getSVLocatableComparator(dictionary);
        final Comparator<Long> idComparator = (o1, o2) -> locusComparator.compare(activeItemLoci.get(o1), activeItemLoci.get(o2));
        this.activeItemIds = new PriorityQueue<>(idComparator);
        this.outputBuffer = new PriorityQueue<>(locusComparator);
        this.nextItemId = 0L;
    }

    public void add(final T item) {
        final SVLocus locus = item.getAsLocus();
        activeItemLoci.put(nextItemId, locus);
        activeItemIds.add(nextItemId);
        engine.add(item, nextItemId++);
    }

    /**
     * Returns any records that can be safely flushed that precede the current minimum starting position
     * of items still being actively clustered. Note this assumes that any future output clusters will not have
     * start position past this item's start position, which should be true for basic
     * {@link org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser.BreakpointSummaryStrategy} types.
     */
    public List<T> flush(final boolean force) {
        final List<C> output = engine.flush(force);
        final Set<Long> engineItemIds = engine.getItemIds();
        final Set<Long> finalizedItemIds = output.stream()
                .map(BasicOutputCluster::getAllIds)
                .flatMap(Collection::stream)
                .filter(id -> !engineItemIds.contains(id))  // Items could still be active in other clusters
                .collect(Collectors.toSet());
        for (final Long id : finalizedItemIds) {
            activeItemIds.remove(id);
            activeItemLoci.remove(id);
        }
        output.stream().map(collapser::collapse).forEachOrdered(outputBuffer::add);
        final SVLocus minActiveStart = activeItemLoci.get(activeItemIds.peek());
        final ArrayList<T> finishedItems = new ArrayList<>();
        while (!outputBuffer.isEmpty() && (minActiveStart == null || locusComparator.compare(outputBuffer.peek(), minActiveStart) < 0)) {
            finishedItems.add(outputBuffer.poll());
        }
        finishedItems.trimToSize();
        return finishedItems;
    }
}
