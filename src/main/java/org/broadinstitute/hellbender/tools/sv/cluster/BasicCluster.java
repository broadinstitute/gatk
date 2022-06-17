package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Container class for clustered items
 */
public class BasicCluster {
    private final Set<Integer> members;
    private int minStart;
    private int maxClusterableStart;

    public BasicCluster(final List<Integer> members, final int minStart, final int maxClusterableStart) {
        Utils.nonNull(members);
        this.members = new HashSet<>(members);
        this.minStart = minStart;
        this.maxClusterableStart = maxClusterableStart;
    }

    public void addMember(final int member, final int memberStart, final int memberMaxClusterableStart) {
        members.add(member);
        maxClusterableStart = Math.max(maxClusterableStart, memberMaxClusterableStart);
        minStart = Math.max(minStart, memberStart);
    }

    public Set<Integer> getMembers() {
        return members;
    }

    public int getMinStart() {
        return minStart;
    }

    public int getMaxClusterableStart() {
        return maxClusterableStart;
    }

    /**
     * Note we do not check for equality on max clusterable start position, which could be dependent on the
     * state of the engine.
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof BasicCluster)) return false;
        BasicCluster cluster = (BasicCluster) o;
        return members.equals(cluster.members);
    }
}
