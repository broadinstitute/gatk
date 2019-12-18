package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.Objects;

/**
 * Graph path class with edges defined in coordinates
 */
public final class CoordinateSVGraphPath {
    private final List<CoordinateSVGraphEdge> edges;

    public CoordinateSVGraphPath(final List<CoordinateSVGraphEdge> edges) {
        Utils.nonNull(edges, "Edges list cannot be null");
        this.edges = edges;
    }

    public List<CoordinateSVGraphEdge> getEdges() {
        return edges;
    }

    /**
     * String format is a series of segments specified as:
     * +[contig:start-end]  -   Forward strand reference
     * +(contig:start-end)  -   Forward strand breakpoint
     * -[contig:start-end]  -   Reverse strand reference
     * -(contig:start-end)  -   Reverse strand breakpoint
     * Note start is always less than end, so the direction of each segment must be inferred in context.
     */
    @Override
    public String toString() {
        final StringBuilder stringBuilder = new StringBuilder(edges.size());
        for (final CoordinateSVGraphEdge edge : edges) {
            if (edge.isInverted()) {
                stringBuilder.append("-");
            } else {
                stringBuilder.append("+");
            }
            if (edge.isReference()) {
                stringBuilder.append("[");
            } else {
                stringBuilder.append("(");
            }
            stringBuilder.append(edge.getContigA());
            stringBuilder.append(":");
            stringBuilder.append(edge.getNodeAPosition());
            stringBuilder.append("-");
            stringBuilder.append(edge.getNodeBPosition());
            if (edge.isReference()) {
                stringBuilder.append("]");
            } else {
                stringBuilder.append(")");
            }
        }
        return stringBuilder.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof CoordinateSVGraphPath)) return false;
        CoordinateSVGraphPath that = (CoordinateSVGraphPath) o;
        return Objects.equals(edges, that.edges);
    }

    @Override
    public int hashCode() {
        return Objects.hash(edges);
    }
}
