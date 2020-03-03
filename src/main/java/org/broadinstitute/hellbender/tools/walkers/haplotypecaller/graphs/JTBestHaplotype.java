package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.JunctionTreeLinkedDeBruijnGraph;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A best haplotype object for being used with junction trees.
 *
 * Each path holds a list of all the junction trees describing its current path. This list consists of pointers to nodes
 * in the junction trees corresponding to the paths that have already been taken by the JTBestHaplotype object.
 *
 * In order to invoke the junction trees simply call {@link #getApplicableNextEdgesBasedOnJunctionTrees} which will return a list
 * of cloned path objects corresponding to each path present in the eldest tree. This method handles popping old trees with insufficient
 * data off of the list as well as incrementing all of the trees in the list to point at the next element based on the chosen path.
 */
public class JTBestHaplotype<V extends BaseVertex, E extends BaseEdge> extends KBestHaplotype<V, E> {
    private JunctionTreeManager junctionTreeManager; // An object for storing and managing operations on the queue of junction trees active for this path
    private int decisionEdgesTakenSinceLastJunctionTreeEvidence;
    private int maxReferenceSpan;

    public boolean isWasPoorlyRecovered() {
        return wasPoorlyRecovered;
    }

    private boolean wasPoorlyRecovered = false;

    // NOTE, this constructor is used by JunctionTreeKBestHaplotypeFinder, in both cases paths are chosen by non-junction tree paths
    public JTBestHaplotype(final JTBestHaplotype<V, E> previousPath, final List<E> edgesToExtend, final double edgePenalty) {
        super(previousPath, edgesToExtend, edgePenalty);
        junctionTreeManager = new JunctionTreeManager(previousPath.junctionTreeManager);
        decisionEdgesTakenSinceLastJunctionTreeEvidence = junctionTreeManager.hasJunctionTreeEvidence() ? 0 : previousPath.decisionEdgesTakenSinceLastJunctionTreeEvidence;
    }

    // Constructor to be used for internal calls from {@link #getApplicableNextEdgesBasedOnJunctionTrees()}
    private JTBestHaplotype(final JTBestHaplotype<V, E> previousPath, final List<E> chain, final int edgeMultiplicity, final int totalOutgoingMultiplicity, final boolean thisPathBasedOnJT) {
        super(previousPath, chain, computeLogPenaltyScore( edgeMultiplicity, totalOutgoingMultiplicity));
        junctionTreeManager = new JunctionTreeManager(previousPath.junctionTreeManager);
        junctionTreeManager.traverseEdgeForAllTrees(chain.get(chain.size() - 1));
        // I'm aware that the chain is only an estimate of the proper length, especially if we got here due to being under the weight threshold for a given tree... the chain lenght is a heuristic as it is...
        decisionEdgesTakenSinceLastJunctionTreeEvidence = thisPathBasedOnJT ? 0 : previousPath.decisionEdgesTakenSinceLastJunctionTreeEvidence + 1;
    }

    // JTBestHaplotype constructor for construction an entirely new haplotype builder.
    public JTBestHaplotype(final V initialVertex, final BaseGraph<V,E> graph) {
        super(initialVertex, graph);
        junctionTreeManager = new JunctionTreeManager();
        decisionEdgesTakenSinceLastJunctionTreeEvidence = 0;
    }

    public boolean hasJunctionTreeEvidence() {
        return junctionTreeManager.hasJunctionTreeEvidence();
    }

    public boolean wasLastEdgeFollowedBasedOnJTEvidence() {
        return decisionEdgesTakenSinceLastJunctionTreeEvidence == 0;
    }

    // returns true if there is a symbolic edge pointing to the reference end or if there is insufficient node data
    public boolean hasStoppingEvidence(final int weightThreshold) {

        // Traverse the non-empty trees until we find one with evidence over our threshold. If we ever find a symbolic end vertex then we stop.
        for (JunctionTreeLinkedDeBruijnGraph.ThreadingNode tree : junctionTreeManager.removeEmptyNodesAndReturnIterator()) {
            int totalOut = getTotalOutForBranch(tree);

            // Are any of these vertexes symbolic stops?
            if (tree.getChildrenNodes().values().stream()
                    .anyMatch(JunctionTreeLinkedDeBruijnGraph.ThreadingNode::isSymbolicEnd)) {
                return true;
            }
            if ( totalOut >= weightThreshold) {
                return false;
            }
        }

        // None of our junction trees cover the stop vertex, close it
        return true;
    }

    // Tally the total outgoing weight for a particular branch
    private static int getTotalOutForBranch(final JunctionTreeLinkedDeBruijnGraph.ThreadingNode eldestTree) {
        return eldestTree == null ? 0 : eldestTree.getChildrenNodes().values().stream()
                .mapToInt(JunctionTreeLinkedDeBruijnGraph.ThreadingNode::getEvidenceCount).sum();
    }

    // Helper method for marking trees as visited
    public void markTreesAsVisited(final List<JunctionTreeLinkedDeBruijnGraph.ThreadingTree> trees) {
        junctionTreeManager.visitedTrees.addAll(trees);
    }


    /**
     * This method is the primary logic of deciding how to traverse junction paths and with what score.
     *
     * TODO this will likely change to use the eldest tree regardless of threshold passage
     * This method checks the list of junction tree nodes, looking first at the eldest tree to perform the following:
     *  - Checks the total outgoing weight, if its below weight threshold then the tree is popped and a new tree is considered
     *  - For each path in the oldest tree clones this path with the chain edges added, taking the edge target for each path present in the tree.
     *
     * @param chain List of edges to add between the current path and the junction tree edge
     * @param weightThreshold threshold of evidence under which old junction trees are discarded.
     * @return A list of new RTBestHaplotypeObjects corresponding to each path chosen from the exisitng junction trees,
     *         or an empty list if there is no path illuminated by junction trees.
     */
    //TODO for reviewer - is this the best way to structure this? I'm not sure how to decide about end nodes based on this, passing them back seesm wrong
    @SuppressWarnings({"unchecked"})
    public List<JTBestHaplotype<V, E>> getApplicableNextEdgesBasedOnJunctionTrees(final List<E> chain, final Set<E> outgoingEdges, final int weightThreshold) {
        Set<MultiSampleEdge> edgesAccountedForByJunctionTrees = new HashSet<>(); // Since we check multiple junction trees for paths, keep track of which paths we have taken to adding duplicate paths to the graph
        List<JTBestHaplotype<V, E>> output = new ArrayList<>();
        for ( JunctionTreeLinkedDeBruijnGraph.ThreadingNode tree : junctionTreeManager.removeEmptyNodesAndReturnIterator()) {
            int totalOut = getTotalOutForBranch(tree);

            // If the total evidence emerging from a given branch
            for (Map.Entry<MultiSampleEdge, JunctionTreeLinkedDeBruijnGraph.ThreadingNode> childNode : tree.getChildrenNodes().entrySet()) {
                if (!outgoingEdges.contains(childNode.getKey())) {
                    throw new GATKException("While constructing graph, there was an incongruity between a JunctionTree edge and the edge present on graph traversal");
                }

                // Don't add edges to the symbolic end vertex here at all, that's handled by {@link #hasStoppingEvidence()}, also don't add the same edge again if we pulled it in from a younger tree.
                if (!childNode.getValue().isSymbolicEnd() && // ignore symbolic end branches, those are handled elsewhere
                        !edgesAccountedForByJunctionTrees.contains(childNode.getKey())) {
                    edgesAccountedForByJunctionTrees.add(childNode.getKey());
                    JunctionTreeLinkedDeBruijnGraph.ThreadingNode child = childNode.getValue();
                    List<E> chainCopy = new ArrayList<>(chain);
                    chainCopy.add((E) childNode.getKey());
                    output.add(new JTBestHaplotype<>(this, chainCopy, child.getEvidenceCount(), totalOut, true));
                }
            }

            // If there isn't enough outgoing evidence, then we poll the next oldest tree for evidence
            // This is done to alleviate the problem that the oldest junction tree may have little evidence and drop connectivity
            // information better represented by one of the younger trees in the path.
            if (totalOut >= weightThreshold) {
                return output;
            }
        }
        // If we hit this point, then eldestTree == null, suggesting that none of the nodes exceeded our threshold for evidence (though some may have found evidence)
        // Standard behavior from the old GraphBasedKBestHaplotypeFinder, base our next path on the edge weights instead
        int totalOutgoingMultiplicity = 0;
        for (final BaseEdge edge : outgoingEdges) {
            totalOutgoingMultiplicity += edge.getMultiplicity();
        }

        // Add all valid edges to the graph
        for (final E edge : outgoingEdges) {
            // Don't traverse an edge if it only has reference evidence supporting it (unless there is no other evidence whatsoever)
            if (!edgesAccountedForByJunctionTrees.contains((MultiSampleEdge) edge) &&
                    totalOutgoingMultiplicity != 0 && edge.getMultiplicity() != 0) {

                // TODO better justify this to myself and others
                // only include ref edges with multiplicity of 1 (i.e. only the ref read spanned it) if there are no other choices at this site (from Junction trees or otherwise)
                if ((edge.isRef() && edge.getMultiplicity() == 1) &&
                        !(edgesAccountedForByJunctionTrees.isEmpty() && outgoingEdges.size() < 2)) { // no junction tree evidence and only one ref path to take
                    continue;
                }

                List<E> chainCopy = new ArrayList<>(chain);
                chainCopy.add(edge);
                output.add(new JTBestHaplotype<>(this, chainCopy, edge.getMultiplicity(), totalOutgoingMultiplicity, false));
            }
        }
        return output;
    }

    /**
     * Return an accounting of how many edge decision (where there is a choice) have consecutively been based on the raw graph weights
     * (as opposed to being based on junction tree evidence).
     *
     * @return number of decision edges
     */
    public int getDecisionEdgesTakenSinceLastJunctionTreeEvidence() {
        return decisionEdgesTakenSinceLastJunctionTreeEvidence;
    }

    /**
     * Add a junction tree (corresponding to the current vertex for traversal, note that if a tree has already been visited by this path then it is ignored)
     * @param junctionTreeForNode Junction tree to add
     */
    public void addJunctionTree(final JunctionTreeLinkedDeBruijnGraph.ThreadingTree junctionTreeForNode) {
        if (junctionTreeManager.addJunctionTree(junctionTreeForNode)) {
            decisionEdgesTakenSinceLastJunctionTreeEvidence = 0;
        }
    }

    /**
     * Add a flag of graph that based on this haplotype we think we should expand the kmer size
     */
    public void setWasPoorlyRecovered(final boolean b) {
        this.wasPoorlyRecovered = b;
    }

    /**
     * A helper class for managing the various junction tree operations that need to be done by JTBestHaplotypeFinder
     *
     * This class tracks traversing all active trees simultaneously so the right corresponding nodes are accounted for in every case.
     * It also keeps track of the previously visited trees so we can save ourselves from double-counting evidence from a particular tree.
     */
    private class JunctionTreeManager {
        Set<JunctionTreeLinkedDeBruijnGraph.ThreadingTree> visitedTrees;
        List<JunctionTreeLinkedDeBruijnGraph.ThreadingNode> activeNodes;

        protected JunctionTreeManager() {
            visitedTrees = new HashSet<>();
            activeNodes = new ArrayList<>(5);
        }

        protected JunctionTreeManager(JunctionTreeManager toCopy) {
            this.visitedTrees = new HashSet<>(toCopy.visitedTrees);
            this.activeNodes = new ArrayList<>(toCopy.activeNodes);
        }

        // Add a junction tree, ensuring that there is a valid tree in order to check.
        // NOTE: this method filters out empty trees or trees that have already been visited on this path
        // Return true if the tree was informative and was actually added
        public boolean addJunctionTree(final JunctionTreeLinkedDeBruijnGraph.ThreadingTree junctionTreeForNode) {
            if (!visitedTrees.contains(junctionTreeForNode) && !junctionTreeForNode.getRootNode().hasNoEvidence()) {
                visitedTrees.add(junctionTreeForNode);
                activeNodes.add(junctionTreeForNode.getRootNode());
                return true;
            }
            return false;
        }

        // method to handle incrementing all of the nodes in the tree simultaneously
        public void traverseEdgeForAllTrees(E edgeTaken) {
            activeNodes = activeNodes.stream()
                    .filter(node -> node.getChildrenNodes().containsKey(edgeTaken))
                    .map(node -> node.getChildrenNodes().get(edgeTaken))
                    .filter(node -> !node.hasNoEvidence())
                    .collect(Collectors.toList());
        }

        // Returns an iterable list of nodes that have sufficient data in the tree (performs pruning of empty nodes)
        public Iterable<JunctionTreeLinkedDeBruijnGraph.ThreadingNode> removeEmptyNodesAndReturnIterator() {
            activeNodes = activeNodes.stream().filter(node -> getTotalOutForBranch(node) > 0).collect(Collectors.toList());
            return activeNodes;
        }

        private JunctionTreeLinkedDeBruijnGraph.ThreadingNode get(int i) {
            return activeNodes.get(i);
        }

        private int size() {
            return activeNodes == null ? 0 : activeNodes.size();
        }

        private void removeOldestTree() {
            activeNodes.remove(0);
        }

        public boolean hasJunctionTreeEvidence() {
            return !activeNodes.isEmpty();
        }
    }
}
