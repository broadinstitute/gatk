package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public final class SVGraphNode {

    private final int contig;
    private final int position;
    private final List<IndexedSVGraphEdge> outEdges;
    private final List<IndexedSVGraphEdge> inEdges;

    public SVGraphNode(final int contig, final int position) {
        this.contig = contig;
        this.position = position;
        this.outEdges = new ArrayList<>(1);
        this.inEdges = new ArrayList<>(1);
    }

    public String toString() {
        return contig + ":" + position;
    }

    public int getContig() {
        return contig;
    }

    public int getPosition() {
        return position;
    }

    public List<IndexedSVGraphEdge> getOutEdges() {
        return outEdges;
    }

    public List<IndexedSVGraphEdge> getInEdges() {
        return inEdges;
    }

    public int countReferenceInEdges() {
        return (int) inEdges.stream().filter(IndexedSVGraphEdge::isReference).count();
    }

    public int countReferenceOutEdges() {
        return (int) outEdges.stream().filter(IndexedSVGraphEdge::isReference).count();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SVGraphNode)) return false;
        SVGraphNode that = (SVGraphNode) o;
        return contig == that.contig &&
                position == that.position &&
                Objects.equals(outEdges, that.outEdges) &&
                Objects.equals(inEdges, that.inEdges);
    }

    @Override
    public int hashCode() {

        return Objects.hash(contig, position, outEdges, inEdges);
    }
}
