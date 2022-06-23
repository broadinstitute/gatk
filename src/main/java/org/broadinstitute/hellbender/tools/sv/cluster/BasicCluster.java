package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Container class for clustered items
 */
public class BasicCluster {
    private final Set<Long> members;
    private String contig;
    private int maxClusterableStart;

    public BasicCluster(final List<Long> members, final String contig, final int maxClusterableStart) {
        Utils.nonNull(members);
        this.members = Utils.nonNull(new HashSet<>(members));
        this.contig = Utils.nonNull(contig);
        this.maxClusterableStart = maxClusterableStart;
    }

    public void addMember(final Long member, final String memberContig, final int memberMaxClusterableStart) {
        Utils.validate(contig.equals(memberContig), "Attempted to add member on contig " + memberContig  + " to cluster on contig " + contig);
        members.add(member);
        maxClusterableStart = Math.max(maxClusterableStart, memberMaxClusterableStart);
    }

    public Set<Long> getMembers() {
        return members;
    }

    public Set<Long> getAllIds() {
        return getMembers();
    }

    public String getContig() {
        return contig;
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
