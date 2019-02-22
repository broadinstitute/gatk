package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.*;

/**
 * Enumerates genotypes on an SV graph
 */
public final class SVGraphGenotyper {

    final SVGraph graph;

    public SVGraphGenotyper(final SVGraph graph) {
        this.graph = graph;
    }

    /**
     * Edge visit for graph traversal
     */
    private static final void visitEdge(final IndexedSVGraphEdge edge, final SVGraphPathState path, final Queue<SVGraphPathState> queue, final int nodeIndex, final int maxEdgeVisits) {
        final IndexedSVGraphEdge augmentingEdge;
        if (edge.causesStrandSwitch()) {
            augmentingEdge = path.isCurrentlyInverted() ? edge.invertedCopy().copyWithUpstreamNodeAs(nodeIndex) : edge.copyWithUpstreamNodeAs(nodeIndex);
        } else {
            augmentingEdge = path.isCurrentlyInverted() ? edge.invertedCopy() : edge.nonInvertedCopy();
        }
        final SVGraphPathState augmentedPath = path.copy();
        if (augmentedPath.addEdge(augmentingEdge, maxEdgeVisits)) {
            if (augmentingEdge.causesStrandSwitch()) {
                augmentedPath.invert();
            }
            queue.add(augmentedPath);
        }
    }


    /**
     * Enumerates genotypes with breadth-first search
     */
    public Collection<IndexedSVGraphPath> enumerate(final double maxPathLengthFactor, final int maxEdgeVisits, final int maxQueueSize) {

        final List<SVGraphNode> nodes = graph.generateNodes();
        final Set<Integer> startingNodes = graph.getStartingNodes();
        final Set<Integer> endingNodes = graph.getEndingNodes();

        final int maxPathLength = (int) (maxPathLengthFactor * graph.getReferenceEdges().size());

        final Set<IndexedSVGraphPath> paths = new HashSet<>();
        final Queue<SVGraphPathState> queue = new ArrayDeque<>();
        for (final Integer startingNode : startingNodes) {
            queue.add(new SVGraphPathState(startingNode, graph.getEdges().size()));
        }
        while (!queue.isEmpty()) {
            if (queue.size() > maxQueueSize) return null;
            final SVGraphPathState path = queue.poll();
            final Integer nodeIndex = path.getCurrentNode();

            //Cap max path length (extremely long paths caused by degenerate loops)
            if (path.getPath().size() > maxPathLength) continue;

            //Finish haplotype if we reach the other end of the event in an uninverted state
            if (!path.isCurrentlyInverted() && endingNodes.contains(nodeIndex)) {
                paths.add(new IndexedSVGraphPath(path.getPath()));
            }

            //Visit neighbors
            final boolean lastEdgeWasReference = path.getPath().isEmpty() || path.getPath().get(path.getPath().size() - 1).isReference();
            final SVGraphNode node = nodes.get(nodeIndex);
            final List<IndexedSVGraphEdge> neighborEdges = path.isCurrentlyInverted() ? node.getInEdges() : node.getOutEdges();
            for (final IndexedSVGraphEdge edge : neighborEdges) {
                //Breakpoint double-jumps are invalid
                if (edge.isReference() || lastEdgeWasReference) {
                    visitEdge(edge, path, queue, nodeIndex, maxEdgeVisits);
                }
            }
        }
        return paths;
    }

    /**
     * Container for the state of a genotype during graph traversal
     */
    private static final class SVGraphPathState {
        private final List<IndexedSVGraphEdge> path;
        private final int[] edgeCopyNumberStates;
        private int currentNode;
        boolean currentlyInverted;

        public SVGraphPathState(final int initialNode, final int numEdges) {
            this.path = new ArrayList<>();
            this.edgeCopyNumberStates = new int[numEdges];
            this.currentNode = initialNode;
            currentlyInverted = false;
        }

        private SVGraphPathState(final int currentNode,
                                 final int[] edgeCopyNumberStates,
                                 final boolean currentlyInverted,
                                 final List<IndexedSVGraphEdge> path) {
            this.currentNode = currentNode;
            this.path = path;
            this.edgeCopyNumberStates = Arrays.copyOf(edgeCopyNumberStates, edgeCopyNumberStates.length);
            this.currentlyInverted = currentlyInverted;
        }

        public boolean isCurrentlyInverted() {
            return currentlyInverted;
        }

        public void invert() {
            currentlyInverted = !currentlyInverted;
        }

        public List<IndexedSVGraphEdge> getPath() {
            return path;
        }

        public int getCurrentNode() {
            return currentNode;
        }

        public boolean addEdge(final IndexedSVGraphEdge edge, final double maxState) {
            getPath().add(edge);
            if (edge.isReference()) {
                edgeCopyNumberStates[edge.getIndex()] += 1;
                if (edgeCopyNumberStates[edge.getIndex()] > maxState) return false;
            }
            currentNode = edge.getDownstreamNodeIndex();
            if (currentNode == -1) {
                throw new GATKException("Attempted to add edge without a set upstream node");
            }
            return true;
        }

        public SVGraphPathState copy() {
            final List<IndexedSVGraphEdge> pathCopy = new ArrayList<>(path.size());
            for (final IndexedSVGraphEdge edge : path) {
                pathCopy.add(edge);
            }
            return new SVGraphPathState(currentNode, edgeCopyNumberStates, currentlyInverted, pathCopy);
        }
    }
}
