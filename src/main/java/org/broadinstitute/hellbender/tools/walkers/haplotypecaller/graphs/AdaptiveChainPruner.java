package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;

public class AdaptiveChainPruner<V extends BaseVertex, E extends BaseEdge> extends ChainPruner<V,E> {
    private final double initialErrorProbability;
    private final double logOddsThreshold;
    private final double seedingLogOddsThreshold;   // threshold for seeding subgraph of good vertices
    private final int maxUnprunedVariants;

    public AdaptiveChainPruner(final double initialErrorProbability, final double logOddsThreshold, final double seedingLogOddsThreshold, final int maxUnprunedVariants) {
        ParamUtils.isPositive(initialErrorProbability, "Must have positive error probability");
        this.initialErrorProbability = initialErrorProbability;
        this.logOddsThreshold = logOddsThreshold;
        this.seedingLogOddsThreshold = seedingLogOddsThreshold;
        this.maxUnprunedVariants = maxUnprunedVariants;
    }

    @Override
    protected Collection<Path<V,E>> chainsToRemove(final List<Path<V, E>> chains) {
        if (chains.isEmpty()) {
            return Collections.emptyList();
        }

        final BaseGraph<V,E> graph = chains.get(0).getGraph();



        Collection<Path<V,E>> probableErrorChains = likelyErrorChains(chains, graph, initialErrorProbability);
        final int errorCount = probableErrorChains.stream().mapToInt(c -> c.getLastEdge().getMultiplicity()).sum();
        final int totalBases = chains.stream().mapToInt(c -> c.getEdges().stream().mapToInt(E::getMultiplicity).sum()).sum();
        final double errorRate = (double) errorCount / totalBases;

        return likelyErrorChains(chains, graph, errorRate).stream().filter(c -> !c.getEdges().stream().anyMatch(BaseEdge::isRef)).collect(Collectors.toList());
    }

    private Collection<Path<V,E>> likelyErrorChains(final List<Path<V, E>> chains, final BaseGraph<V,E> graph, final double errorRate) {
        // pre-compute the left and right log odds of each chain
        final Map<Path<V,E>, Pair<Double, Double>> chainLogOdds = chains.stream()
                .collect(Collectors.toMap(c -> c, c-> chainLogOdds(c, graph, errorRate)));

        // compute correspondence of vertices to incident chains with log odds above the seeding and extending thresholds
        final Multimap<V, Path<V,E>> vertexToSeedableChains = ArrayListMultimap.create();
        final Multimap<V, Path<V,E>> vertexToGoodIncomingChains = ArrayListMultimap.create();
        final Multimap<V, Path<V,E>> vertexToGoodOutgoingChains = ArrayListMultimap.create();

        for (final Path<V,E> chain : chains) {
            if (chainLogOdds.get(chain).getRight() >= logOddsThreshold) {
                vertexToGoodIncomingChains.put(chain.getLastVertex(), chain);
            }

            if (chainLogOdds.get(chain).getLeft() >= logOddsThreshold) {
                vertexToGoodOutgoingChains.put(chain.getFirstVertex(), chain);
            }

            // seed-worthy chains must pass the more stringent seeding log odds threshold on both sides
            // in addition to that, we only seed from vertices with multiple such chains incoming or outgoing (see below)
            if (chainLogOdds.get(chain).getRight() >= seedingLogOddsThreshold && chainLogOdds.get(chain).getLeft() >= seedingLogOddsThreshold) {
                vertexToSeedableChains.put(chain.getFirstVertex(), chain);
                vertexToSeedableChains.put(chain.getLastVertex(), chain);
            }
        }

        // find a subset of good vertices from which to grow the subgraph of good chains.
        final Queue<V> verticesToProcess = new ArrayDeque<>();
        final Path<V,E> maxWeightChain = getMaxWeightChain(chains);
        verticesToProcess.add(maxWeightChain.getFirstVertex());
        verticesToProcess.add(maxWeightChain.getLastVertex());

        // look for vertices with two incoming or two outgoing chains (plus one outgoing or incoming for a total of 3 or more) with good log odds to seed the subgraph of good vertices
        // the logic here is that a high-multiplicity error chain A that branches into a second error chain B and a continuation-of-the-original-error chain A'
        // may have a high log odds for A'.  However, only in the case of true variation will  multiple branches leaving the same vertex have good log odds.
        vertexToSeedableChains.keySet().stream()
                .filter(v -> vertexToSeedableChains.get(v).size() > 2)
                .forEach(verticesToProcess::add);

        final Set<V> processedVertices = new LinkedHashSet<>();
        final Set<Path<V,E>> goodChains = new LinkedHashSet<>();

        // starting from the high-confidence seed vertices, grow the "good" subgraph along chains with above-threshold log odds,
        // discovering good chains as we go.
        while (!verticesToProcess.isEmpty()) {
            final V vertex = verticesToProcess.poll();
            processedVertices.add(vertex);
            for (final Path<V,E> outgoingChain : vertexToGoodOutgoingChains.get(vertex)) {
                goodChains.add(outgoingChain);
                if(!processedVertices.contains(outgoingChain.getLastVertex())) {
                    verticesToProcess.add(outgoingChain.getLastVertex());
                }
            }

            for (final Path<V,E> incomingChain : vertexToGoodIncomingChains.get(vertex)) {
                goodChains.add(incomingChain);
                if(!processedVertices.contains(incomingChain.getFirstVertex())) {
                    verticesToProcess.add(incomingChain.getFirstVertex());
                }
            }
        }

        final Set<Path<V,E>> errorChains = chains.stream().filter(c -> !goodChains.contains(c)).collect(Collectors.toSet());

        // A vertex with N > 0 outgoing good chains corresponds to N - 1 variants
        int numberOfVariantsInGraph = vertexToGoodOutgoingChains.keySet().stream()
                .mapToInt(v -> Math.max(vertexToGoodOutgoingChains.get(v).size() - 1, 0)).sum();

        if (numberOfVariantsInGraph > maxUnprunedVariants) {
            errorChains.addAll(goodChainsToPruneToSatisfyMaxVariants(chainLogOdds, goodChains));
        }

        return errorChains;
    }

    /**
     * If there are too many apparent variants in the graph we do extra pruning of the worst chains until we have few enough variants.
     *
     * @param chainLogOdds  two-sided log odds for each chain in the graph
     * @param goodChains    chains that were initially considered good.
     */
    private Set<Path<V, E>> goodChainsToPruneToSatisfyMaxVariants(Map<Path<V, E>, Pair<Double, Double>> chainLogOdds, Set<Path<V, E>> goodChains) {
        final Set<Path<V, E>> newlyPrunedChains = new HashSet<>();

        // define the subgraph of good chains -- we're going to shrink this
        final Multimap<V, Path<V, E>> incoming = ArrayListMultimap.create();
        final Multimap<V, Path<V, E>> outgoing = ArrayListMultimap.create();
        for (final Path<V,E> goodChain : goodChains) {
            outgoing.put(goodChain.getFirstVertex(), goodChain);
            incoming.put(goodChain.getLastVertex(), goodChain);
        }

        // find "subgraph chains", that is, maximal non-branching paths within the subgraph of good chains.
        // These may comprise multiple good chains that were split by error chains.
        final List<List<Path<V,E>>> subgraphChains = new ArrayList<>();
        for (final Path<V,E> startChain : goodChains) {
            // find chains that start subgraph chains
            if (incoming.get(startChain.getFirstVertex()).size() == 1 && outgoing.get(startChain.getFirstVertex()).size() == 1) {
                continue;
            }

            final List<Path<V,E>> subgraphChain = new ArrayList<>();
            V leftVertex = startChain.getFirstVertex();

            do {
                Path<V, E> chain = outgoing.get(leftVertex).iterator().next();
                subgraphChain.add(chain);
                leftVertex = chain.getLastVertex();
            } while (incoming.get(leftVertex).size() == 1 && outgoing.get(leftVertex).size() == 1);

            subgraphChains.add(subgraphChain);
        }

        // done up to here

        int numberOfVariantsInGraph = outgoing.keySet().stream().mapToInt(v -> Math.max(outgoing.get(v).size() - 1, 0)).sum();

        // start with the worst good variants
        final PriorityQueue<Path<V,E>> pruningQueue = new PriorityQueue<>(
                Comparator.comparingDouble((Path<V,E> c) -> Math.min(chainLogOdds.get(c).getLeft(), chainLogOdds.get(c).getRight()))
                        .thenComparing((Path<V,E> c) -> c.getFirstVertex().getSequence(), BaseUtils.BASES_COMPARATOR));

        pruningQueue.addAll(goodChains);

        while( !pruningQueue.isEmpty() && numberOfVariantsInGraph > maxUnprunedVariants) {

            final Path<V,E> worstGoodChain = pruningQueue.poll();
            if (newlyPrunedChains.contains(worstGoodChain)) {   // it's possible we already deleted this chain when extending an earlier chain
                continue;
            }

            // This chain may be part of a larger non-branching path within the good subgraph, all of which we delete.

            // first, move to the start of the subgraph chain
            V leftVertex = worstGoodChain.getFirstVertex();
            while (incoming.get(leftVertex).size() == 1 && outgoing.get(leftVertex).size() == 1) {
                leftVertex = incoming.get(leftVertex).iterator().next().getFirstVertex();
            }

            // now delete every chain in the subgraph chain
            do {
                Path<V, E> chainToRemove = outgoing.get(leftVertex).iterator().next();
                final V rightVertex = chainToRemove.getLastVertex();
                newlyPrunedChains.add(chainToRemove);
                outgoing.remove(leftVertex, chainToRemove);
                incoming.remove(rightVertex, chainToRemove);
                leftVertex = rightVertex;
            } while (incoming.get(leftVertex).isEmpty() && outgoing.get(leftVertex).size() == 1);
            numberOfVariantsInGraph--;
        }
        return newlyPrunedChains;
    }

    // find the chain containing the edge of greatest weight, taking care to break ties deterministically
    private Path<V, E> getMaxWeightChain(final Collection<Path<V, E>> chains) {
        return chains.stream()
                .max(Comparator.comparingInt((Path<V, E> chain) -> chain.getEdges().stream().mapToInt(BaseEdge::getMultiplicity).max().orElse(0))
                        .thenComparingInt(Path::length)
                        .thenComparing((Path<V,E> c) -> c.getFirstVertex().getSequence(), BaseUtils.BASES_COMPARATOR))
                .get();
    }

    // left and right chain log odds
    private Pair<Double, Double> chainLogOdds(final Path<V,E> chain, final BaseGraph<V,E> graph, final double errorRate) {
        final int leftTotalMultiplicity = MathUtils.sumIntFunction(graph.outgoingEdgesOf(chain.getFirstVertex()), E::getMultiplicity);
        final int rightTotalMultiplicity = MathUtils.sumIntFunction(graph.incomingEdgesOf(chain.getLastVertex()), E::getMultiplicity);

        final int leftMultiplicity = chain.getEdges().get(0).getMultiplicity();
        final int rightMultiplicity = chain.getLastEdge().getMultiplicity();

        final double leftLogOdds = graph.isSource(chain.getFirstVertex()) ? 0.0 :
                Mutect2Engine.logLikelihoodRatio(leftTotalMultiplicity - leftMultiplicity, leftMultiplicity, errorRate);
        final double rightLogOdds = graph.isSink(chain.getLastVertex()) ? 0.0 :
                Mutect2Engine.logLikelihoodRatio(rightTotalMultiplicity - rightMultiplicity, rightMultiplicity, errorRate);

        return ImmutablePair.of(leftLogOdds, rightLogOdds);
    }

}
