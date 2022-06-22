package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

public class PartitionedCluster extends BasicCluster {

    final Long primaryItem;

    public PartitionedCluster(final List<Long> members,
                              final String contig,
                              final int maxClusterableStart,
                              final Long primaryItem) {
        super(members, contig, maxClusterableStart);
        this.primaryItem = Utils.nonNull(primaryItem);
    }

    public Long getPrimaryItem() {
        return primaryItem;
    }
}
