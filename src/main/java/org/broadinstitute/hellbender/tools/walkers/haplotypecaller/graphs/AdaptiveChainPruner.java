package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
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
        /////// NEW

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

        final Set<V> processedVertices = new HashSet<>();
        final Set<Path<V,E>> goodChains = new HashSet<>();

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

        // add non-error chains to error chains if maximum number of variants has been exceeded
        chains.stream()
                .filter(c -> !errorChains.contains(c))
                .filter(c -> isChainPossibleVariant(c, graph))
                .sorted(Comparator.comparingDouble(c -> Math.min(chainLogOdds.get(c).getLeft(), chainLogOdds.get(c).getLeft())).reversed())
                .skip(maxUnprunedVariants)
                .forEach(errorChains::add);

        return errorChains;

    }

    // find the chain containing the edge of greatest weight, taking care to break ties deterministically
    private Path<V, E> getMaxWeightChain(final Collection<Path<V, E>> chains) {
        return chains.stream()
                .max(Comparator.comparingInt((Path<V, E> chain) -> chain.getEdges().stream().mapToInt(BaseEdge::getMultiplicity).max().orElse(0))
                        .thenComparingInt(Path::length)
                        .thenComparingInt(Path::hashCode)).get();
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
    
    private boolean isChainPossibleVariant(final Path<V,E> chain, final BaseGraph<V,E> graph) {
        final int leftTotalMultiplicity = MathUtils.sumIntFunction(graph.outgoingEdgesOf(chain.getFirstVertex()), E::getMultiplicity);
        final int rightTotalMultiplicity = MathUtils.sumIntFunction(graph.incomingEdgesOf(chain.getLastVertex()), E::getMultiplicity);

        final int leftMultiplicity = chain.getEdges().get(0).getMultiplicity();
        final int rightMultiplicity = chain.getLastEdge().getMultiplicity();

        return leftMultiplicity <= leftTotalMultiplicity / 2 || rightMultiplicity <= rightTotalMultiplicity / 2;
    }
}
