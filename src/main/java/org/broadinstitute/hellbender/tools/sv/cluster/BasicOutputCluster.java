package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

public class BasicOutputCluster<T extends SVLocatable> {
    private final List<T> members;

    public BasicOutputCluster(final List<T> members) {
        Utils.nonNull(members);
        this.members = members;
    }

    public List<T> getMembers() {
        return members;
    }
}
