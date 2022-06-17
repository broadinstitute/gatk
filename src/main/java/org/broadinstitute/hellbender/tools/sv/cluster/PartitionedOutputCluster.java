package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

public class PartitionedOutputCluster<T extends SVLocatable> extends BasicOutputCluster<T> {

    final T primaryItem;

    public PartitionedOutputCluster(final List<T> members, final T primaryItem) {
        super(members);
        Utils.nonNull(primaryItem);
        this.primaryItem = primaryItem;
    }

    public T getPrimaryItem() {
        return primaryItem;
    }
}
