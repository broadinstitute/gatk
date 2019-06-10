package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ExperimentalReadThreadingGraph;
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
public class JTBestHaplotype<T extends BaseVertex, E extends BaseEdge> extends KBestHaplotype<T, E> {
    // TODO these nodes probably shouldn't be exposed and held externally methinks
    private List<ExperimentalReadThreadingGraph.ThreadingNode> activeNodes;

    public JTBestHaplotype(final JTBestHaplotype p, final List<E> edgesToExtend, final double edgePenalty) {
        super(p, edgesToExtend, edgePenalty);
        //this.threadingTreesToConsult = p.threadingTreesToConsult.cop;
        activeNodes = new ArrayList<ExperimentalReadThreadingGraph.ThreadingNode>(p.activeNodes);
    }

    // Constructor to be used for internal calls from {@link #getApplicableNextEdgesBasedOnJunctionTrees()}
    public JTBestHaplotype(final JTBestHaplotype p, final List<E> chain, final int edgeMultiplicity, final int totalOutgoingMultiplicity) {
        super(p, chain, computeLogPenaltyScore( edgeMultiplicity, totalOutgoingMultiplicity));
        activeNodes = new ArrayList<ExperimentalReadThreadingGraph.ThreadingNode>(p.activeNodes);
        // Ensure that the relevant edge has been traversed
        takeEdge(chain.get(chain.size() - 1));
    }

    public JTBestHaplotype(final T initialVertex, final BaseGraph<T,E> graph) {
        super(initialVertex, graph);
        activeNodes = new ArrayList<>();
    }

    //TODO this needs to be the same logic as the blow method, this is temporary
    // returns true if there is a symbolic edge pointing to the reference end or if there is insufficient node data
    public boolean hasStoppingEvidence(final int weightThreshold) {
        ExperimentalReadThreadingGraph.ThreadingNode eldestTree = activeNodes.isEmpty() ? null : activeNodes.get(0);
        int totalOut = getTotalOutForBranch(eldestTree);

        // Keep removing trees until we find one under our threshold TODO this should be in a helper method
        while (eldestTree != null && totalOut < weightThreshold) {
            activeNodes.remove(0);
            eldestTree = activeNodes.isEmpty() ? null : activeNodes.get(0);
            totalOut = getTotalOutForBranch(eldestTree);
        }

        if (eldestTree != null) {
            for ( ExperimentalReadThreadingGraph.ThreadingNode node : eldestTree.getChildrenNodes().values()) {
                if (node.isSymbolicEnd()) {
                    return true;
                }
            }
            return false;
        }
        return true;
    }


    private int getTotalOutForBranch(ExperimentalReadThreadingGraph.ThreadingNode eldestTree) {
        int totalOut = 0;
        if (eldestTree != null) {
            for (ExperimentalReadThreadingGraph.ThreadingNode node : eldestTree.getChildrenNodes().values()) {
                totalOut += node.getCount();
            }
        }
        return totalOut;
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
    public List<JTBestHaplotype<T, E>> getApplicableNextEdgesBasedOnJunctionTrees(final List<E> chain, final int weightThreshold) {
        List<JTBestHaplotype<T, E>> output = new ArrayList<>();
        ExperimentalReadThreadingGraph.ThreadingNode eldestTree = activeNodes.isEmpty() ? null : activeNodes.get(0);
        while (eldestTree != null) {
            //TODO this can be better, need to create a tree "view" object that tracks the current node more sanely
            int totalOut = getTotalOutForBranch(eldestTree);
            // This right here is what handles dealing with weight thresholds
            // TODO for the reviewer - Perhaps the approach here is to check younger trees as well but don't throw out informative older trees
            if (totalOut >= weightThreshold) {
                //TODO add SOME sanity check to ensure that the vertex we stand on and the edges we are polling line up
                for (Map.Entry<MultiSampleEdge, ExperimentalReadThreadingGraph.ThreadingNode> childNode : eldestTree.getChildrenNodes().entrySet()) {
                    // Don't add edges to the symbolic end vertex here at all, thats handled elsewhere
                    if (!childNode.getValue().isSymbolicEnd()) {
                        ExperimentalReadThreadingGraph.ThreadingNode child = childNode.getValue();
                        List<E> chainCopy = new ArrayList<>(chain);
                        chainCopy.add((E) childNode.getKey());
                        output.add(new JTBestHaplotype<>(this, chainCopy, child.getCount(), totalOut));
                    }
                }
                return output;

            // If there aren't enough outgoing nodes we just remove
            } else {
                activeNodes.remove(0);
                eldestTree = activeNodes.isEmpty() ? null : activeNodes.get(0);
            }
        }
        return output;
    }

    // method to handle incrementing all of the nodes in the tree
    private void takeEdge(final E edgeTaken) {
        activeNodes = activeNodes.stream().map(node -> {
            if (!node.getChildrenNodes().containsKey(edgeTaken)) {
                return null;
            }
            return node.getChildrenNodes().get(edgeTaken);
        }).filter(Objects::nonNull).collect(Collectors.toList());
    }

    public void addJunctionTree(final ExperimentalReadThreadingGraph.ThreadingTree junctionTreeForNode) {
        activeNodes.add(junctionTreeForNode.getRootNode());
    }
}
