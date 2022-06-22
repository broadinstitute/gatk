package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;

public class PartitionedOutputCluster<T extends SVLocatable> extends BasicOutputCluster<T> {

    final T primaryItem;

    public PartitionedOutputCluster(final Map<Long, T> members, final T primaryItem) {
        super(members);
        Utils.nonNull(primaryItem);
        this.primaryItem = primaryItem;
    }

    public T getPrimaryItem() {
        return primaryItem;
    }
}
