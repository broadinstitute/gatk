package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class AdaptiveChainPruner<V extends BaseVertex, E extends BaseEdge> extends ChainPruner<V,E> {
    private final double initialErrorProbability;
    private final double logOddsThreshold;
    private final int maxUnprunedVariants;

    public AdaptiveChainPruner(final double initialErrorProbability, final double logOddsThreshold, final int maxUnprunedVariants) {
        ParamUtils.isPositive(initialErrorProbability, "Must have positive error probability");
        this.initialErrorProbability = initialErrorProbability;
        this.logOddsThreshold = logOddsThreshold;
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

        return likelyErrorChains(chains, graph, errorRate);
    }

    private Collection<Path<V,E>> likelyErrorChains(final List<Path<V, E>> chains, final BaseGraph<V,E> graph, final double errorRate) {
        final Map<Path<V,E>, Double> chainLogOdds = chains.stream()
                .collect(Collectors.toMap(c -> c, c-> chainLogOdds(c, graph, errorRate)));

        final Set<Path<V,E>> result = new HashSet<>(chains.size());

        chainLogOdds.forEach((chain, lod) -> {
            if (lod < logOddsThreshold) {
                result.add(chain);
            }
        });

        chains.stream().filter(c -> isChainPossibleVariant(c, graph))
                .sorted(Comparator.comparingDouble(chainLogOdds::get).reversed())
                .skip(maxUnprunedVariants)
                .forEach(result::add);

        return result;

    }

    private double chainLogOdds(final Path<V,E> chain, final BaseGraph<V,E> graph, final double errorRate) {
        if (chain.getEdges().stream().anyMatch(E::isRef)) {
            return Double.POSITIVE_INFINITY;
        }

        final int leftTotalMultiplicity = MathUtils.sumIntFunction(graph.outgoingEdgesOf(chain.getFirstVertex()), E::getMultiplicity);
        final int rightTotalMultiplicity = MathUtils.sumIntFunction(graph.incomingEdgesOf(chain.getLastVertex()), E::getMultiplicity);

        final int leftMultiplicity = chain.getEdges().get(0).getMultiplicity();
        final int rightMultiplicity = chain.getLastEdge().getMultiplicity();

        final double leftLogOdds = graph.isSource(chain.getFirstVertex()) ? 0.0 :
                Mutect2Engine.lnLikelihoodRatio(leftTotalMultiplicity - leftMultiplicity, leftMultiplicity, errorRate);
        final double rightLogOdds = graph.isSink(chain.getLastVertex()) ? 0.0 :
                Mutect2Engine.lnLikelihoodRatio(rightTotalMultiplicity - rightMultiplicity, rightMultiplicity, errorRate);

        return FastMath.max(leftLogOdds, rightLogOdds);
    }

    // is the chain
    private boolean isChainPossibleVariant(final Path<V,E> chain, final BaseGraph<V,E> graph) {
        final int leftTotalMultiplicity = MathUtils.sumIntFunction(graph.outgoingEdgesOf(chain.getFirstVertex()), E::getMultiplicity);
        final int rightTotalMultiplicity = MathUtils.sumIntFunction(graph.incomingEdgesOf(chain.getLastVertex()), E::getMultiplicity);

        final int leftMultiplicity = chain.getEdges().get(0).getMultiplicity();
        final int rightMultiplicity = chain.getLastEdge().getMultiplicity();

        return leftMultiplicity <= leftTotalMultiplicity / 2 || rightMultiplicity <= rightTotalMultiplicity / 2;
    }
}
