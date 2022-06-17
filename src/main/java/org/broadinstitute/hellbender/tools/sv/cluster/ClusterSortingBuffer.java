package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.BoundType;
import com.google.common.collect.SortedMultiset;
import com.google.common.collect.TreeMultiset;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

final class ClusterSortingBuffer<R extends BasicOutputCluster> {
    private SortedMultiset<R> buffer;

    public ClusterSortingBuffer(final Comparator<R> clusterComparator) {
        this.buffer = TreeMultiset.create(clusterComparator);
    }

    public void add(final R cluster) {
        buffer.add(cluster);
    }

    /**
     * Returns any records that can be safely flushed based on the current minimum starting position
     * of items still being actively clustered.
     */
    public List<R> flush(final R minActiveStartingPositionCluster) {
        if (buffer.isEmpty()) {
            return Collections.emptyList();
        }
        if (minActiveStartingPositionCluster == null) {
            forceFlush();
        }
        final SortedMultiset<R> finalizedView = buffer.headMultiset(minActiveStartingPositionCluster, BoundType.CLOSED);
        final ArrayList<R> finalizedClusters = new ArrayList<>(finalizedView);
        // Clearing a view of the buffer also clears the items from the buffer itself
        finalizedView.clear();
        return finalizedClusters;
    }

    /**
     * Returns all buffered records, regardless of any active clusters. To be used only when certain that no
     * active clusters can be clustered with any future inputs.
     */
    public List<R> forceFlush() {
        final List<R> result = new ArrayList<>(buffer);
        buffer.clear();
        return result;
    }

    public boolean isEmpty() {
        return buffer.isEmpty();
    }
}
