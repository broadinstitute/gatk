package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public final class IndexedSVGraphPath implements Comparable<IndexedSVGraphPath> {

    private final List<IndexedSVGraphEdge> edges;

    public IndexedSVGraphPath(final List<IndexedSVGraphEdge> edges) {
        this.edges = edges;
    }

    public List<IndexedSVGraphEdge> getEdges() {
        return edges;
    }

    public CoordinateSVGraphPath convertToCoordinatePath(final SVGraph graph) {
        final List<CoordinateSVGraphEdge> coordinateEdges = new ArrayList<>(edges.size());
        for (final IndexedSVGraphEdge edge : edges) {
            coordinateEdges.add(new CoordinateSVGraphEdge(edge, graph.getNodes()));
        }
        return new CoordinateSVGraphPath(coordinateEdges);
    }

    public boolean isReference() {
        for (final IndexedSVGraphEdge edge : edges) {
            if (!edge.isReference()) {
                return false;
            }
        }
        return true;
    }

    public int size() {
        return edges.size();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof IndexedSVGraphPath)) return false;
        IndexedSVGraphPath that = (IndexedSVGraphPath) o;
        return Objects.equals(edges, that.edges);
    }

    @Override
    public int hashCode() {
        return Objects.hash(edges);
    }

    @Override
    public int compareTo(final IndexedSVGraphPath other) {
        if (this.size() > other.size()) return 1;
        if (this.size() < other.size()) return -1;
        for (int i = 0; i < edges.size(); i++) {
            if (i >= other.size()) return 1;
            int result = edges.get(i).compareTo(other.getEdges().get(i));
            if (result != 0) {
                return result;
            }
        }
        return 0;
    }
}
