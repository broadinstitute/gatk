package org.broadinstitute.hellbender.tools.sv.cluster;

import java.util.List;

public class PartitionedCluster extends BasicCluster {

    final int primaryItem;

    public PartitionedCluster(final List<Integer> members,
                              final int minStart,
                              final int maxClusterableStart,
                              final int primaryItem) {
        super(members, minStart, maxClusterableStart);
        this.primaryItem = primaryItem;
    }

    public Integer getPrimaryItem() {
        return primaryItem;
    }
}
