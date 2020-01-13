package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.AbstractReadThreadingGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.JunctionTreeLinkedDeBruinGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;

public class JunctionTreeKBestHaplotypeFinder<V extends BaseVertex, E extends BaseEdge> extends KBestHaplotypeFinder<V, E> {
    public static final int DEFAULT_MAX_UPSTREAM_REFERENCE_JUMP_TO_ALLOW = 40;
    public static final int DEFAULT_NUM_UPSTREAM_REFERENCE_EDGES_TO_TOLERATE = 10;
    public static final int DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE = 3;
    public static final int DEFAULT_MINIMUM_WEIGHT_FOR_JT_BRANCH_TO_NOT_BE_PRUNED = 1;
    public static final int DEFAULT_MAX_ACCEPTABLE_DECISION_EDGES_WITHOUT_JT_GUIDANCE = 5;
    public static final int DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH = 1;
    // Workarounds for relatively rare complex sites that loop pathologically and do not generate ending paths
    public static final int DEFAULT_MAX_PATHS_TO_CONSIDER_WITHOUT_RESULT = 1000;
    public static final int DEFAULT_MAX_PATHS_TO_EVER_CONSIDER = 10000;
    private int weightThreshold = DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE;

    // List for mapping vertexes that start chains of kmers that do not diverge, used to cut down on repeated graph traversal
    Map<V, List<E>> graphKmerChainCache = new HashMap<>();

    // Graph to be operated on, in this case cast as an JunctionTreeLinkedDeBruinGraph
    JunctionTreeLinkedDeBruinGraph junctionTreeLinkedDeBruinGraph;
    final boolean experimentalEndRecoveryMode;

    public JunctionTreeKBestHaplotypeFinder(final BaseGraph<V, E> graph, final Set<V> sources, Set<V> sinks, final int branchWeightThreshold, final boolean experimentalEndRecoveryMode) {
        super(sinks, sources, graph);
        Utils.validate(graph instanceof JunctionTreeLinkedDeBruinGraph, "JunctionTreeKBestHaplotypeFinder requires an JunctionTreeLinkedDeBruinGraph be provided");
        junctionTreeLinkedDeBruinGraph = (JunctionTreeLinkedDeBruinGraph) graph;
        this.experimentalEndRecoveryMode = experimentalEndRecoveryMode;
        ParamUtils.isPositive(weightThreshold, "Pruning Weight Threshold must be a positive number greater than 0");
    }

    /**
     * Constructor for the special case of a single source and sink
     */
    public JunctionTreeKBestHaplotypeFinder(final BaseGraph<V, E> graph, final V source, final V sink) {
        this(graph, Collections.singleton(source), Collections.singleton(sink), DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE, false);
    }

    /**
     * Constructor for the special case of a single source and sink
     */
    public JunctionTreeKBestHaplotypeFinder(final BaseGraph<V, E> graph, final V source, final V sink, final int branchWeightThreshold, final boolean endRecoveryMode) {
        this(graph, Collections.singleton(source), Collections.singleton(sink), branchWeightThreshold, endRecoveryMode);
    }

    /**
     * Constructor for the default case of all sources and sinks
     */
    public JunctionTreeKBestHaplotypeFinder(final BaseGraph<V, E> graph) {
        this(graph, graph.getReferenceSourceVertex(), graph.getReferenceSinkVertex());
    }

    @Override
    public boolean keepCycles() {
        return true;
    }

    @VisibleForTesting
    public JunctionTreeKBestHaplotypeFinder<V, E> setWeightThreshold(final int outgoingWeight) {
        Utils.validate(weightThreshold > 0, "Pruning Weight Threshold must be a positive number greater than 0");
        weightThreshold = outgoingWeight;
        return this;
    }

    /**
     * The primary engine for finding haplotypes based on junction trees.
     * <p>
     * Best haplotypes are discovered by a modified Djikstra's algorithm. Paths are genearated from the reference start
     * kmer and proceed forward until they encounter one of the following stop conditions and are treated accordingly:
     *  (1) Encounter a kmer with a JunctionTree - In this case the junction tree is added to a JunctionTree queue maintained
     *                                             by the KBestHaplotypePath.
     *  (2) Encounter a kmer with out degree > 1 - In this case we consult the oldest JunctionTree in the experimental path,
     *                                             and add new KBestHaplotypePaths to the Queue for each branch on the tree
     *                                             for that node with weights calculated based on edge weight in the junction tree.
     *                                             If the eldest tree has insufficient data then it is popped off the queue and a younger
     *                                             tree is consulted, if no younger tree is consulted then the raw graph edge
     *                                             weights for the fork are used to generate paths instead.
     *  (3) Encounter a ReferenceSink - If a reference sink is encountered, close out the path and add it to the results. TODO this will be subjected to change
     *
     * Note: A cache is maintained of vertexes to contiguous sequences of paths in order to cut down on repeated traversal
     *       of the same segments of the graph for each path.
     *
     * @param maxNumberOfHaplotypes maximum number of haplotypes to discover.
     * @return A list of the best scoring haplotypes for hte GraphBasedKBestHaplotypeFinder
     */
    @Override
    @SuppressWarnings({"unchecked"})
    public List<KBestHaplotype<V, E>> findBestHaplotypes(final int maxNumberOfHaplotypes) {
        //pre-process step: find pivotal edges so they can be marked off as visited
        final LinkedHashSet<E> unvisitedPivotalEdges = createMapOfPivotalEdgesInTopologicalOrder();

        final List<JTBestHaplotype<V, E>> result = new ArrayList<>();
        final PriorityQueue<JTBestHaplotype<V, E>> queue = new PriorityQueue<>(Comparator.comparingDouble(KBestHaplotype<V, E>::score).reversed());
        sources.forEach(source -> queue.add(new JTBestHaplotype<>(source, graph)));

        // Iterate over paths in the queue, unless we are out of paths of maxHaplotypes to find
        while (result.size() < maxNumberOfHaplotypes && (!queue.isEmpty() || !unvisitedPivotalEdges.isEmpty())) {
            // TODO this may change at some point
            // stopgap to handle edge cases (save ourselves the risk of infinite looping)
            if (result.isEmpty() ?
                    queue.size() > DEFAULT_MAX_PATHS_TO_CONSIDER_WITHOUT_RESULT : // restrict the number of branching paths examined
                    queue.size() > DEFAULT_MAX_PATHS_TO_EVER_CONSIDER) {
                break;
            }

            // breakout condition, pop a new path onto the tree from unvisited pivotal edges if
            if ( queue.isEmpty() ) {
                final E firstEdge = unvisitedPivotalEdges.stream().findFirst().get();
                unvisitedPivotalEdges.remove(firstEdge);
                // TODO this may change when the logic changes... but it seems like a no-brainer check for now
                // Check at this stage that are not starting a path that can never succeed
                final V pivotalVerex = graph.getEdgeSource(firstEdge);
                Optional<JTBestHaplotype<V, E>> bestMatchingHaplotype = result.stream().filter(path -> path.containsVertex(pivotalVerex)).max(Comparator.comparingDouble(JTBestHaplotype::score));
                if (bestMatchingHaplotype.isPresent()) {
                    JTBestHaplotype<V, E> pathToAdd = constructArtificialHaplotypePath(firstEdge, pivotalVerex, bestMatchingHaplotype.get());
                    if (pathToAdd != null) {
                        queue.add(pathToAdd);
                    }
                }
                continue;
            }

            final JTBestHaplotype<V, E> pathToExtend = queue.poll();

            // This safeguards against infinite loops and degenerate excessively long paths, only allow 4 decisions without junction tree guidance
            if (pathToExtend.getDecisionEdgesTakenSinceLastJunctionTreeEvidence() > DEFAULT_MAX_ACCEPTABLE_DECISION_EDGES_WITHOUT_JT_GUIDANCE) {
                continue;
            }
            ////////////////////////////////////////////////////////////
            // code to discover where the next interesting node is
            ////////////////////////////////////////////////////////////
            V vertexToExtend = pathToExtend.getLastVertex();
            Set<E> outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);

            // Check if we have cached the current vertex in a contiguous sequence
            final List<E> chain = graphKmerChainCache.computeIfAbsent(vertexToExtend, k -> new ArrayList<>());
            // if not, step forward until a vertex meets one of conditions (1), (2), or (3) are met
            if (chain.isEmpty()) {
                // Keep going until we reach a fork, reference sink, or fork
                while (outgoingEdges.size() == 1 && // Case (2)
                        !junctionTreeLinkedDeBruinGraph.getJunctionTreeForNode((MultiDeBruijnVertex) vertexToExtend).isPresent() && // Case (1)
                        !sinks.contains(vertexToExtend))// Case (3)
                {
                    final E edge = outgoingEdges.iterator().next();
                    // Defensive check for looping edges, don't extend the chain in such a way as to create loops, let other code handle that
                    if (chain.contains(edge)) {
                        break;
                    }
                    chain.add(edge);
                    vertexToExtend = graph.getEdgeTarget(edge);
                    outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);
                }
                // Cache the chain result
                graphKmerChainCache.put(pathToExtend.getLastVertex(), chain);

            } else {
                // We have already expanded this part of the graph, use it.
                vertexToExtend = graph.getEdgeTarget(chain.get(chain.size() - 1));
                outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);
            }
            // vertexToExtend and outgoingEdges are necessarily at the next "interesting" point


            ////////////////////////////////////////////////////////////
            // code to decide what to do at that interesting node
            ////////////////////////////////////////////////////////////
            // In the event we have a junction tree on top of a vertex with outDegree > 1, we add this first before we traverse paths
            junctionTreeLinkedDeBruinGraph.getJunctionTreeForNode((MultiDeBruijnVertex) vertexToExtend).ifPresent(pathToExtend::addJunctionTree);

            //TODO this can probabaly be 100% consumed by getApplicableNextEdgesBasedOnJunctionTrees() as a check... that would simplify things somewhat
            // If we are at a reference end then we close out the path
            if (sinks.contains(vertexToExtend) && pathToExtend.hasStoppingEvidence(weightThreshold)) {
                //TODO this will probably be resolved using a junction tree on that node and treating it as an edge to extend
                //todo the proposal here would be to check if there is an active tree left for us at this point and if so keep going
                JTBestHaplotype<V, E> newPath =
                        reconcilePathMissingReferenceStartPositions(
                                chain.isEmpty() ? pathToExtend : new JTBestHaplotype<>(pathToExtend, chain, 0),
                                result,
                                graph);
                annotatePathBasedOnGraph(newPath, junctionTreeLinkedDeBruinGraph);
                result.add(newPath);
                pathToExtend.getEdges().forEach(e -> unvisitedPivotalEdges.remove(e));
            }
            // NOTE: even if we are at the reference stop and there is evidence in the junction trees of a stop we still want to explore other edges potentially

            // We must be at a point where the path diverges, use junction trees to resolve if possible
            if (outgoingEdges.size() > 1) {
                List<JTBestHaplotype<V, E>> jTPaths = pathToExtend.getApplicableNextEdgesBasedOnJunctionTrees(chain, outgoingEdges, weightThreshold);
                if (jTPaths.isEmpty() && !sinks.contains(vertexToExtend)) {
//                    throw new GATKException("Found no path based on the junction trees or exisiting paths, this should not have happened");
                    System.out.println("Found nothing Queue has this many: " + queue.size() + "\nPath that failed to extend was junction tree: " + pathToExtend.getVertices());
                }
                // Filter out paths that involve the same kmer too many times (if we were directed by a junction tree or have trees to follow then don't worry about repeated kmers)
                List<JTBestHaplotype<V, E>> filteredPaths = jTPaths.stream()
                        .filter(path -> path.hasJunctionTreeEvidence() || path.wasLastEdgeFollowedBasedOnJTEvidence() ||
                                        // Count the number of occurances of the laset vertex, if there are more than DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH throw away the path
                                        path.getVertices().stream().filter(v -> v.equals(path.getLastVertex())).count() <= DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH)
                        .collect(Collectors.toList());
                if (jTPaths.isEmpty() && !sinks.contains(vertexToExtend)) {
//                    throw new GATKException("Found no path based on the junction trees or exisiting paths, this should not have happened");
                    System.out.println("A path was filtered because it was looping without junction tree support");
                }

                queue.addAll(filteredPaths);

                // Otherwise just take the next node forward
            } else {
                // If there are no outgoing edges from this node, then just kill this branch from the queue
                // NOTE: this branch in the future might be responsible for dangling end merging
                if (outgoingEdges.size() > 0) {
                    //TODO evaluate the expense of asking this quesion, there are ways to mitigate the cost. This particuar case is almost always triggered
                    // Defensive check, if we see the same vertex
                    final V finalVertexToExtend = vertexToExtend;
                    if (!pathToExtend.hasJunctionTreeEvidence() &&
                            pathToExtend.getVertices()
                                    .stream()
                                    .filter(v -> v.equals(finalVertexToExtend))
                                    .count() > DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH) {
                        // do nothing
                    } else {
                        // otherwie add the path
                        List<E> chainCopy = new ArrayList<>(chain);
                        chainCopy.add(outgoingEdges.iterator().next());
                        queue.add(new JTBestHaplotype<>(pathToExtend, chainCopy, 0));
                    }
                }
            }
        }

        return result.stream().map(n -> (KBestHaplotype<V, E>) n).collect(Collectors.toList());
    }

    //TODO thiss might be best computed in the paths as they are being expanded
    private void annotatePathBasedOnGraph(final JTBestHaplotype<V, E> newPath, final JunctionTreeLinkedDeBruinGraph graph) {
        int farthestReferenceEdgeReached = 0;
        int numUpstreamRefEdgesEncountered = 0;
        int lastReferenceEdgeVisited = 0;
        for (E edge : newPath.getEdges()) {
            List<Integer> refOccurances = ((MultiSampleEdge)edge).getReferencePathIndexes();

            // if we are not on a reference edge don't worry
            if (!refOccurances.isEmpty()) {
                // find the next lowest ref occurance that is incrementally higher than our current one (assumes sorted ref occurrences list)
                int refIndex = 0;
                for (Integer index : refOccurances) {
                    refIndex = index;
                    if (refIndex > lastReferenceEdgeVisited) {
                        break;
                    }
                }

                // check if we are too far upstream
                if (farthestReferenceEdgeReached > refIndex + DEFAULT_MAX_UPSTREAM_REFERENCE_JUMP_TO_ALLOW) {
                    numUpstreamRefEdgesEncountered++;
                } else if (refIndex > farthestReferenceEdgeReached) {
                    farthestReferenceEdgeReached = refIndex;
                }
                lastReferenceEdgeVisited = refIndex;
            }
        }

        // if we saw too many out of sequence reference edges then we report it as a potentially bad haplotype
        if (numUpstreamRefEdgesEncountered > DEFAULT_NUM_UPSTREAM_REFERENCE_EDGES_TO_TOLERATE) {
            newPath.setIsWonky(true);
        }
    }

    /**
     * Method that performs the stitching of another paths begining onto the front of a new path that is being recovered after not being discovered.
     *
     * @param firstEdge
     * @param pivotalVerex
     * @param bestMatchingHaplotype
     * @return
     */
    private JTBestHaplotype<V, E> constructArtificialHaplotypePath(final E firstEdge, final V pivotalVerex, final JTBestHaplotype<V, E> bestMatchingHaplotype) {
        // Now we try to construct a reference covering haplotype from the one we just discovered
        List<E> bestMatchingHaplotypeEdges = bestMatchingHaplotype.getEdges();
        List<E> edgesIncomingToSplitPoint = bestMatchingHaplotypeEdges.stream().filter(edge -> graph.getEdgeTarget(edge).equals(pivotalVerex)).collect(Collectors.toList());
        //todo it is either an error state to find nothing or could mean we accidentally did the search over the source vertex, either way shoudl be a bug
        if (edgesIncomingToSplitPoint.isEmpty()) {
            return null;
        }

        // TODO maybe this will matter some day, simply select the last edge
        List<E> edgesBeforeSplit = new ArrayList<>(bestMatchingHaplotypeEdges.subList(0, bestMatchingHaplotypeEdges.lastIndexOf(edgesIncomingToSplitPoint.get(edgesIncomingToSplitPoint.size() - 1)) + 1));
        edgesBeforeSplit.add(firstEdge);

        // create a new path with the beginging of the best edge stapled to the front
        JTBestHaplotype<V, E> pathToAdd = new JTBestHaplotype<>(new JTBestHaplotype<>(bestMatchingHaplotype.getFirstVertex(), graph), edgesBeforeSplit, bestMatchingHaplotype.score());
        List<JunctionTreeLinkedDeBruinGraph.ThreadingTree> treesPassed = pathToAdd.getVertices().stream()
                .map(v -> junctionTreeLinkedDeBruinGraph.getJunctionTreeForNode((MultiDeBruijnVertex) v))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());
        pathToAdd.markTreesAsVisited(treesPassed);
        return pathToAdd;
    }


    /**
     * Helper method that takes BestHaplotypePaths that may or may not start at a valid reference start position
     *
     * @return true if a valid reference-starting path was able to be constructed.
     */
    // TODO maybe if this fails to find we should include the reference path explicitly
    private JTBestHaplotype<V, E> reconcilePathMissingReferenceStartPositions(JTBestHaplotype<V, E> pathToReconcile, List<JTBestHaplotype<V, E>> validReturnPaths, BaseGraph<V, E> graph) {
        // check that the path is valid, if so don't modify it
        if (  sources.contains(pathToReconcile.getVertices().get(0))) {
            return pathToReconcile;
        }

        throw new RuntimeException("e");
//
//        V pivotalVerex = pathToReconcile.getFirstVertex();
//        //TODO this can be MUCH faster than a simple contains search here
//        List<JTBestHaplotype<V, E>> candidatePaths = validReturnPaths.stream().filter(path -> path.containsVertex(pivotalVerex)).collect(Collectors.toList());
//
//        // todo, perhaps something more drastic can be done here, this arises from either uncovered reference path or from branches that lead to loops that are unresolvable... perhaps try to capture the first one
//        if (candidatePaths.isEmpty()) {
//            // this is a failure state for now
//            return null;
//        }
//
//        //TODO this is totally simple for right now, will choose a better approach soon.
//        JTBestHaplotype<V, E> bestMatchingHaplotype = candidatePaths.stream().max(Comparator.comparingDouble(JTBestHaplotype::score)).get();
//
//        // Now we try to construct a reference covering haplotype from the one we just discovered
//        List<E> bestMatchingHaplotypeEdges = bestMatchingHaplotype.getEdges();
//        List<E> edgesIncomingToSplitPoint = bestMatchingHaplotypeEdges.stream().filter(edge -> graph.getEdgeTarget(edge).equals(pivotalVerex)).collect(Collectors.toList());
//        //todo it is either an error state to find nothing or could mean we accidentally did the search over the source vertex, either way shoudl be a bug
//        if (edgesIncomingToSplitPoint.isEmpty()) {
//            return null;
//        }
//
//        // TODO maybe this will matter some day, simply select the last edge
//        List<E> outputEdges = new ArrayList<>(bestMatchingHaplotypeEdges.subList(0, bestMatchingHaplotypeEdges.lastIndexOf(edgesIncomingToSplitPoint.get(edgesIncomingToSplitPoint.size() - 1)) + 1));
//        outputEdges.addAll(pathToReconcile.getEdges());
//
//        return new JTBestHaplotype<V,E>(new JTBestHaplotype<V,E>(bestMatchingHaplotype.getFirstVertex(), graph), outputEdges, bestMatchingHaplotype.score() + pathToReconcile.score());
    }


    /**
     * Method responsible for handling the missed edges from the junction tree graph.
     * <p>
     * New behavior is as follows, this creates a map of "pivotal edges" based on topographical order, thowing away those that have been visited
     */
    private void initializeTogpographicalMap() {

    }

    private boolean isPivotalEdge(E edge) {
        return graph.outgoingEdgesOf(graph.getEdgeSource(edge)).size() > 1;
    }


    /**
     * This method is used as a pre-processing step in order to find all of the pivotal edges in the graph in a sorted fashion.
     *
     * In order to accomplish this, the entire graph is traversed starting at the reference inputs and edges that might
     * or might not be considered by the KBestHaplotypeFinder.
     *
     * Note: Currently "Topological order is defined as the number of edges taken from the start, which is accomplished with
     *       a breadth first search.
     *
     * @return returns a LinkedHashSet object where all of the pivotal edges in the graph will be iterated in
     */
    //TODO for this first implementation I have chosen (for simplicity sake) to base the event on the
    //TODO this can be optimized, don't go down this rabbit hole too eagerly
    @VisibleForTesting
    LinkedHashSet<E> createMapOfPivotalEdgesInTopologicalOrder() {
        final Set<E> visitedEdges = new HashSet<>(); // used to save ourselves the trouble of excessive graph traversal (and repetative)
        final PriorityQueue<tinyEdgeHelper> edgesToVisit = new PriorityQueue<>(Comparator.comparingDouble(tinyEdgeHelper::score));
        final LinkedHashSet<E> outputEdgesInOrder = new LinkedHashSet<>();

        // zips through the entire graph, searching for pivotal edges and adding them based on the number of edges since the referecne
        // TODO decide on true topography or distance from reference. To add confusion, how do you guarintee the topographical order is visited correctly?

        // Initialize the graph with the start
        edgesToVisit.addAll(sources.stream()
                .flatMap(s -> graph.outgoingEdgesOf(s).stream())
                .map(e -> new tinyEdgeHelper(e, 0))
                .collect(Collectors.toList()));

        while (!edgesToVisit.isEmpty()) {
            tinyEdgeHelper nextEdge = edgesToVisit.poll();
            List<E> outgoingEdges = new ArrayList<>(graph.outgoingEdgesOf(graph.getEdgeTarget(nextEdge.edge)));

            // If there are multiple interesting outgoing edges mark them as pivotal (this is currently non-deterministic and maybe could be improved?)
            if (outgoingEdges.size() > 1) {
                // not that uncovered reference edges do not get added to the output (because they were discounted from the graph for a reason)
                outputEdgesInOrder.addAll(outgoingEdges.stream().filter(e -> !(e.isRef() && e.getMultiplicity() == 1) && !visitedEdges.contains(e)).collect(Collectors.toList()));
            }

            // visit the unvisited edges and add any edges not already found to the output linked hash map
            outgoingEdges.stream()
                    .filter(e -> !visitedEdges.contains(e))
                    .forEach(e -> {
                        edgesToVisit.add(new tinyEdgeHelper(e, nextEdge.score + 1));
                        visitedEdges.add(e);
                    });
        }

        return outputEdgesInOrder;
    }


//        while (!edgesToVisit.isEmpty()) {
//            tinyEdgeHelper startEdge = edgesToVisit.poll();
//            E currentEdge = startEdge.edge;
//            int count = startEdge.score;
//
//            // follow
//            while (currentEdge != null) {
//                // don't visit the same edge twice
//                if (visitedEdges.contains(currentEdge)) {
//                    continue;
//                }
//                visitedEdges.add(currentEdge);
//                count++;
//                V currentVertex = graph.getEdgeTarget(currentEdge);
//                Set<E> outgoingEdges = graph.outgoingEdgesOf(currentVertex);
//                // Follow the single track while we can
//                if (outgoingEdges.size() == 1) {
//                    currentEdge = outgoingEdges.stream().findFirst().get();
//                } else if (outgoingEdges.size() > 0) {
//                    Iterator<E> edges = outgoingEdges.iterator();
//                    currentEdge = edges.next();
//                    if (!paintedEdges.contains(currentEdge)) {
//                        outputQueue.add()
//                    }
//                    while (edges.hasNext()) {
//                        outputQueue.add(new tinyEdgeHelper(edges.next(), count));
//                    }
//                } else {
//                    // Then we have reached the end
//                }
//            }
//        }

    //todo figure out if this is necessary
    // TODO also this could be greatly optimized if we just remember where the suspect edges are during the first pass...
    private class tinyEdgeHelper {
        E edge;
        int score;

        public tinyEdgeHelper(E e, int i) {
            edge = e;
            score = score();
        }

        public int score() {
            return score;
        }

    }
}
