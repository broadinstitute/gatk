package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

public class BasicOutputCluster<T extends SVLocatable> {
    private final Map<Long, T> members;

    public BasicOutputCluster(final Map<Long, T> members) {
        Utils.nonNull(members);
        this.members = members;
    }

    public Collection<T> getMembers() {
        return members.values();
    }

    public Set<Long> getMemberIds() {
        return members.keySet();
    }
}
