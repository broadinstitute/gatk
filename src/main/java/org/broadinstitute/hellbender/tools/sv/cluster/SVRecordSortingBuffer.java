package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVLocatable;

import java.util.*;
import java.util.stream.Collectors;

public final class SVRecordSortingBuffer<R extends SVLocatable> {
    private TreeMap<SVLocatable, R> buffer;

    public SVRecordSortingBuffer(final Comparator<SVLocatable> comparator) {
        this.buffer = new TreeMap<>(comparator);
    }

    public void add(final R item) {
        buffer.put(item.getAsLocus(), item);
    }

    /**
     * Returns any records that can be safely flushed that precede the current minimum starting position
     * of items still being actively clustered. Note this assumes that any future output clusters will not have
     * start position past this item's start position, which should be true for basic
     * {@link org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser.BreakpointSummaryStrategy} types.
     */
    public List<R> flush(final SVLocatable minActivePositionLocus) {
        if (buffer.isEmpty()) {
            return Collections.emptyList();
        }
        if (minActivePositionLocus == null) {
            return forceFlush();
        }
        final List<Map.Entry<SVLocatable, R>> finalizedItems = new ArrayList<>(buffer.headMap(minActivePositionLocus).entrySet());
        finalizedItems.stream().map(Map.Entry::getKey).forEach(buffer::remove);
        return finalizedItems.stream().map(Map.Entry::getValue).collect(Collectors.toList());
    }

    /**
     * Returns all buffered records, regardless of position. To be used only when certain that no future
     * items will come before any currently in the buffer.
     */
    private List<R> forceFlush() {
        final List<R> result = new ArrayList<>(buffer.values());
        buffer.clear();
        return result;
    }

    public boolean isEmpty() {
        return buffer.isEmpty();
    }
}
