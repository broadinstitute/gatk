package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.JunctionTreeLinkedDeBruinGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;

public class JunctionTreeKBestHaplotypeFinder<V extends BaseVertex, E extends BaseEdge> extends KBestHaplotypeFinder<V, E> {
    public static final int DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE = 3;
    public static final int DEFAULT_MINIMUM_WEIGHT_FOR_JT_BRANCH_TO_NOT_BE_PRUNED = 2;
    public static final int DEFAULT_MAX_ACCEPTABLE_DECISION_EDGES_WITHOUT_JT_GUIDANCE = 7;
    public static final int DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH = 10;
    private int weightThreshold = DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE;

    // List for mapping vertexes that start chains of kmers that do not diverge, used to cut down on repeated graph traversal
    Map<V, List<E>> graphKmerChainCache = new HashMap<>();

    // Graph to be operated on, in this case cast as an JunctionTreeLinkedDeBruinGraph
    JunctionTreeLinkedDeBruinGraph junctionTreeLinkedDeBruinGraph;

    public JunctionTreeKBestHaplotypeFinder(final BaseGraph<V, E> graph, final Set<V> sources, Set<V> sinks, final int branchWeightThreshold) {
        super(sinks, sources, graph);
        Utils.validate(graph instanceof JunctionTreeLinkedDeBruinGraph, "JunctionTreeKBestHaplotypeFinder requires an JunctionTreeLinkedDeBruinGraph be provided");
        junctionTreeLinkedDeBruinGraph = (JunctionTreeLinkedDeBruinGraph) graph;
        ParamUtils.isPositive(weightThreshold, "Pruning Weight Threshold must be a positive number greater than 0");
    }

    /**
     * Constructor for the special case of a single source and sink
     */
    public JunctionTreeKBestHaplotypeFinder(final BaseGraph<V, E> graph, final V source, final V sink) {
        this(graph, Collections.singleton(source), Collections.singleton(sink), DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE);
    }

    /**
     * Constructor for the special case of a single source and sink
     */
    public JunctionTreeKBestHaplotypeFinder(final BaseGraph<V, E> graph, final V source, final V sink, final int branchWeightThreshold) {
        this(graph, Collections.singleton(source), Collections.singleton(sink), branchWeightThreshold);
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
     *
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
        final List<JTBestHaplotype<V, E>> result = new ArrayList<>();
        final PriorityQueue<JTBestHaplotype<V, E>> queue = new PriorityQueue<>(Comparator.comparingDouble(KBestHaplotype<V, E>::score).reversed());
        sources.forEach(source -> queue.add(new JTBestHaplotype<>(source, graph)));

        // Iterate over paths in the queue, unless we are out of paths of maxHaplotypes to find
        while (!queue.isEmpty() && result.size() < maxNumberOfHaplotypes) {
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
            if (chain.isEmpty()){
                // Keep going until we reach a fork, reference sink, or fork
                while ( outgoingEdges.size() == 1 && // Case (2)
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
                if (chain.isEmpty()) {
                    result.add(pathToExtend);
                } else {
                    result.add(new JTBestHaplotype<>(pathToExtend, chain, 0));
                }
            }
            // NOTE: even if we are at the reference stop and there is evidence in the junction trees of a stop we still want to explore other edges potentially

            // We must be at a point where the path diverges, use junction trees to resolve if possible
            if (outgoingEdges.size() > 1) {
                List<JTBestHaplotype<V, E>> jTPaths = pathToExtend.getApplicableNextEdgesBasedOnJunctionTrees(chain, outgoingEdges, weightThreshold);
                if (jTPaths.isEmpty() && !sinks.contains(vertexToExtend)) {
//                    throw new GATKException("Found no path based on the junction trees or exisiting paths, this should not have happened");
                    System.out.println("Found nothing Queue has this many: "+queue.size()+"\nPath that failed to extend was junction tree: "+pathToExtend.getVertices());
                }
                queue.addAll(jTPaths);

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

}
