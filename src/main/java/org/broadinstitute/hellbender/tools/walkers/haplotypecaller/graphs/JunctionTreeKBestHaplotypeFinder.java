package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.JunctionTreeLinkedDeBruijnGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;

public class JunctionTreeKBestHaplotypeFinder<V extends BaseVertex, E extends BaseEdge> extends KBestHaplotypeFinder<V, E> {
    private Logger logger = LogManager.getLogger(getClass());

    public static final int DEFAULT_MAX_UPSTREAM_REFERENCE_JUMP_TO_ALLOW = 40;
    public static final int DEFAULT_NUM_UPSTREAM_REFERENCE_EDGES_TO_TOLERATE = 10;
    public static final int DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE = 3;
    public static final int DEFAULT_MINIMUM_WEIGHT_FOR_JT_BRANCH_TO_NOT_BE_PRUNED = 1;
    public static final int DEFAULT_MAX_ACCEPTABLE_DECISION_EDGES_WITHOUT_JT_GUIDANCE = 5;
    public static final int DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH = 1;
    // Workarounds for relatively rare complex sites that loop pathologically and do not generate ending paths
    public static final int DEFAULT_MAX_PATHS_TO_CONSIDER_WITHOUT_RESULT = 1000;
    public static final int DEFAULT_MAX_PATHS_TO_EVER_CONSIDER = 10000;
    private int junctionTreeEvidenceWeightThreshold = DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE;

    // List for mapping vertexes that start chains of kmers that do not diverge, used to cut down on repeated graph traversal
    private Map<V, List<E>> graphKmerChainCache = new HashMap<>();

    // Graph to be operated on, in this case cast as an JunctionTreeLinkedDeBruijnGraph
    private JunctionTreeLinkedDeBruijnGraph junctionTreeLinkedDeBruijnGraph;
    private final boolean experimentalEndRecoveryMode;

    public JunctionTreeKBestHaplotypeFinder(final BaseGraph<V, E> graph, final Set<V> sources, Set<V> sinks, final int branchWeightThreshold, final boolean experimentalEndRecoveryMode) {
        super(sinks, sources, graph);
        Utils.validate(graph instanceof JunctionTreeLinkedDeBruijnGraph, "JunctionTreeKBestHaplotypeFinder requires an JunctionTreeLinkedDeBruijnGraph be provided");
        this.junctionTreeLinkedDeBruijnGraph = (JunctionTreeLinkedDeBruijnGraph) graph;
        this.experimentalEndRecoveryMode = experimentalEndRecoveryMode;
        ParamUtils.isPositive(junctionTreeEvidenceWeightThreshold, "Pruning Weight Threshold must be a positive number greater than 0");
    }

    /**
     * Constructor for the special case of a single source and sink
     */
    public JunctionTreeKBestHaplotypeFinder(final BaseGraph<V, E> graph, final V source, final V sink) {
        this(graph, Collections.singleton(source), Collections.singleton(sink), DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE, true);
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
    public JunctionTreeKBestHaplotypeFinder<V, E> setJunctionTreeEvidenceWeightThreshold(final int outgoingWeight) {
        Utils.validate(junctionTreeEvidenceWeightThreshold > 0, "Pruning Weight Threshold must be a positive number greater than 0");
        junctionTreeEvidenceWeightThreshold = outgoingWeight;
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
     *  (3) Encounter a ReferenceSink - If a reference sink is encountered, close out the path and add it to the results. TODO: this will be subjected to change
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
        //pre-process step: find pivotal edges so they can be marked off as visited (if we want to recover edges uncovered in the graph).
        final LinkedHashSet<E> unvisitedPivotalEdges = experimentalEndRecoveryMode ? createMapOfPivotalEdgesInTopologicalOrder() : new LinkedHashSet<>();

        final List<JTBestHaplotype<V, E>> result = new ArrayList<>();
        final PriorityQueue<JTBestHaplotype<V, E>> queue = new PriorityQueue<>(Comparator.comparingDouble(KBestHaplotype<V, E>::score).reversed());
        sources.forEach(source -> queue.add(new JTBestHaplotype<>(source, graph)));

        // Iterate over paths in the queue, unless we are out of paths of maxHaplotypes to find
        while (result.size() < maxNumberOfHaplotypes && (!queue.isEmpty() || !unvisitedPivotalEdges.isEmpty())) {
            // check that we aren't caught in a hopelessly complicated graph for which we can't hope to recover
            if (queue.size() > (result.isEmpty() ? DEFAULT_MAX_PATHS_TO_CONSIDER_WITHOUT_RESULT : DEFAULT_MAX_PATHS_TO_EVER_CONSIDER)) {
                break;
            }

            // breakout condition, pop a new path onto the tree from unvisited pivotal edges if
            if ( queue.isEmpty() ) {
                enqueueNextPivotalEdge(unvisitedPivotalEdges, result, queue);
                continue;
            }

            final JTBestHaplotype<V, E> pathToExtend = queue.poll();

            // This safeguards against infinite loops and degenerate excessively long paths, only allow 4 decisions without junction tree guidance
            if (pathToExtend.getDecisionEdgesTakenSinceLastJunctionTreeEvidence() > DEFAULT_MAX_ACCEPTABLE_DECISION_EDGES_WITHOUT_JT_GUIDANCE) {
                continue;
            }
            ////////////////////////////////////////////////////////////
            // code to discover where the next interesting node is (or use the cached result)
            ////////////////////////////////////////////////////////////
            V vertexToExtend = pathToExtend.getLastVertex();
            Set<E> outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);

            // Check if we have cached the current vertex in a contiguous sequence
            final List<E> chain = graphKmerChainCache.computeIfAbsent(vertexToExtend, k -> new ArrayList<>());
            // if not, step forward until a vertex meets one of conditions (1), (2), or (3) are met
            if (chain.isEmpty()) {
                // Keep going until we reach a fork, reference sink, or fork
                while (outgoingEdges.size() == 1 && // Case (2)
                        !junctionTreeLinkedDeBruijnGraph.getJunctionTreeForNode((MultiDeBruijnVertex) vertexToExtend).isPresent() && // Case (1)
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
            junctionTreeLinkedDeBruijnGraph.getJunctionTreeForNode((MultiDeBruijnVertex) vertexToExtend).ifPresent(pathToExtend::addJunctionTree);

            //TODO this can probabaly be 100% consumed by getApplicableNextEdgesBasedOnJunctionTrees() as a check... that would simplify things somewhat
            // If we are at a reference end then we close out the path
            if (sinks.contains(vertexToExtend) && pathToExtend.hasStoppingEvidence(junctionTreeEvidenceWeightThreshold)) {
                //TODO this will probably be resolved using a junction tree on that node and treating it as an edge to extend
                //todo the proposal here would be to check if there is an active tree left for us at this point and if so keep going
                JTBestHaplotype<V, E> newPath = reconcilePathMissingReferenceStartPositions(chain.isEmpty() ?
                        pathToExtend :  new JTBestHaplotype<>(pathToExtend, chain, 0));
                // check that we were able to recover the missing path
                if (newPath != null) {
                    //TODO this code corresponds to where we check how well the path matches with the reference path, its not currently enabled but left in as we will do further evaluations
                    //annotatePathBasedOnGraph(newPath, junctionTreeLinkedDeBruijnGraph);
                    result.add(newPath);
                }
                pathToExtend.getEdges().forEach(unvisitedPivotalEdges::remove);
            }
            // NOTE: even if we are at the reference stop and there is evidence in the junction trees of a stop we still want to explore other edges potentially

            // We must be at a point where the path diverges, use junction trees to resolve if possible
            if (outgoingEdges.size() > 1) {
                List<JTBestHaplotype<V, E>> jTPaths = pathToExtend.getApplicableNextEdgesBasedOnJunctionTrees(chain, outgoingEdges, junctionTreeEvidenceWeightThreshold);
                if (jTPaths.isEmpty() && !sinks.contains(vertexToExtend)) {
                    logger.debug("Found nothing Queue has this many: " + queue.size() + "\nPath that failed to extend was junction tree: " + pathToExtend.getVertices());
                }
                // Filter out paths that involve the same kmer too many times without permission from a junction tree
                List<JTBestHaplotype<V, E>> filteredPaths = jTPaths.stream()
                        .filter(path -> path.hasJunctionTreeEvidence() || path.wasLastEdgeFollowedBasedOnJTEvidence() ||
                                // Count the number of occurrences of the latest vertex, if there are more than DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH throw away the path
                                path.getVertices().stream().filter(v -> v.equals(path.getLastVertex())).count() <= DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH)
                        .collect(Collectors.toList());
                if (jTPaths.isEmpty() && !sinks.contains(vertexToExtend)) {
                    logger.debug("A path was filtered because it was looping without junction tree support");
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

    /**
     * Helper method that controls the logic for pivotal (ie. edge where a decision is made) edges.
     *
     * The logic for pivotal edges currently is this:
     *  1: pop the next pivotal edge in our tree
     *  2: search our results set for any paths that cross the uncovered vertex, choose the one with the highest score
     *  3: if one is found, search for the last occurrence of the vertex for the pivotal edge in the path.
     *  4: construct an "artificial" haplotype that consists of all the edges of the chosen path up to the pivotal vertex
     *     with the chosen pivotal edge appended to the end. Add this to the provided queue.
     *
     *
     * NOTE: there is a limitation to this approach, while appending the paths at the front is faster and simpler as it allows
     *       the new "artificial" haplotype to remember what edges it has visited, it does not necessarily mean that the chosen
     *       path is the closest match to the resulting haplotype... See {@link JunctiontreeKbesthaplotypeFinderUnitTest.testRecoveryOfDroppedPathChoosingMostLikePathDespiteThatPathHavingAWorseScore()}
     *       for an illustration of this problem.
     *
     * @param unvisitedPivotalEdges ordered list of edges to try connecting
     * @param result                completed paths in the graph to use for construction
     * @param queue                 path priority queue to deposit new edges into
     */
    private void enqueueNextPivotalEdge(LinkedHashSet<E> unvisitedPivotalEdges, List<JTBestHaplotype<V, E>> result, PriorityQueue<JTBestHaplotype<V, E>> queue) {
        final E firstEdge = unvisitedPivotalEdges.stream().findFirst().get();
        final V pivotalVerex = graph.getEdgeSource(firstEdge);
        unvisitedPivotalEdges.remove(firstEdge);

        // Check at this stage that are not starting a path that can never succeed //TODO this might cost a lot of runtime, check in profiler
        Optional<JTBestHaplotype<V, E>> bestMatchingHaplotype = result.stream().filter(path -> path.containsVertex(pivotalVerex)).max(Comparator.comparingDouble(JTBestHaplotype::score));
        if (bestMatchingHaplotype.isPresent()) {
            // Now we try to construct a reference covering haplotype from the one we just discovered
            final List<E> bestMatchingHaplotypeEdges = bestMatchingHaplotype.get().getEdges();
            final List<E> edgesIncomingToSplitPoint = bestMatchingHaplotypeEdges.stream().filter(edge -> graph.getEdgeTarget(edge).equals(pivotalVerex)).collect(Collectors.toList());
            //todo it is either an error state to find nothing or could mean we accidentally did the search over the source vertex, either way shoudl be a bug
            if (edgesIncomingToSplitPoint.isEmpty()) {
                return;
            }

            // From that best haplotype we choose the last occurrence of the pivotal branch vertex as representative
            // TODO maybe this will matter some day, simply select the last edge
            List<E> edgesBeforeSplit = new ArrayList<>(bestMatchingHaplotypeEdges.subList(0, bestMatchingHaplotypeEdges.lastIndexOf(edgesIncomingToSplitPoint.get(edgesIncomingToSplitPoint.size() - 1)) + 1));
            edgesBeforeSplit.add(firstEdge);

            // create a new path with the beginging of the best edge stapled to the front
            JTBestHaplotype<V, E> pathToAdd = new JTBestHaplotype<>(new JTBestHaplotype<>(bestMatchingHaplotype.get().getFirstVertex(), graph), edgesBeforeSplit, bestMatchingHaplotype.get().score());
            List<JunctionTreeLinkedDeBruijnGraph.ThreadingTree> treesPassed = pathToAdd.getVertices().stream()
                    .map(v -> junctionTreeLinkedDeBruijnGraph.getJunctionTreeForNode((MultiDeBruijnVertex) v))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.toList());
            pathToAdd.markTreesAsVisited(treesPassed);
            queue.add(pathToAdd);
        }
    }

    //TODO this probably needs to be more testing...
    //TODO thiss might be best computed in the paths as they are being expanded
    private void annotatePathBasedOnGraph(final JTBestHaplotype<V, E> newPath, final JunctionTreeLinkedDeBruijnGraph graph) {
        int farthestReferenceEdgeReached = 0;
        int numUpstreamRefEdgesEncountered = 0;
        int lastReferenceEdgeVisited = 0;
        for (E edge : newPath.getEdges()) {
            List<Integer> refOccurances = ((MultiSampleEdge)edge).getReferencePathIndexes();

            // if we are not on a reference edge don't worry
            if (!refOccurances.isEmpty()) {
                // find the next lowest ref occurrence that is incrementally higher than our current one (assumes sorted ref occurrences list)
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
            newPath.setWasPoorlyRecovered(true);
        }
    }


    /**
     * Helper method that takes BestHaplotypePaths that may or may not start at a valid reference start position
     *
     * @return true if a valid reference-starting path was able to be constructed.
     */
    // TODO maybe if this fails to find we should include the reference path explicitly
    private JTBestHaplotype<V, E> reconcilePathMissingReferenceStartPositions(final JTBestHaplotype<V, E> pathToReconcile) {
        // check that the path is valid, if so don't modify it
        if (  sources.contains(pathToReconcile.getVertices().get(0))) {
            return pathToReconcile;
        }

        // TODO this might change if we want to construct better paths in the future
        throw new GATKException.ShouldNeverReachHereException("It looks like we have tried to merge a haplotype to the output that does not start on the reference source. This should not have happened");
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
        final PriorityQueue<TinyEdgeHelper> edgesToVisit = new PriorityQueue<>(Comparator.comparingDouble(TinyEdgeHelper::score));
        final LinkedHashSet<E> outputEdgesInOrder = new LinkedHashSet<>();

        // zips through the entire graph, searching for pivotal edges and adding them based on the number of edges since the referecne
        // Initialize the graph with the start
        edgesToVisit.addAll(sources.stream()
                .flatMap(s -> graph.outgoingEdgesOf(s).stream())
                .map(e -> new TinyEdgeHelper(e, 0))
                .collect(Collectors.toList()));

        while (!edgesToVisit.isEmpty()) {
            TinyEdgeHelper nextEdge = edgesToVisit.poll();
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
                        edgesToVisit.add(new TinyEdgeHelper(e, nextEdge.score + 1));
                        visitedEdges.add(e);
                    });
        }

        return outputEdgesInOrder;
    }


    private class TinyEdgeHelper {
        E edge;
        int score;

        private TinyEdgeHelper(E e, int i) {
            edge = e;
            score = score();
        }

        private int score() {
            return score;
        }
    }
}
