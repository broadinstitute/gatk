package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Prune all chains from this graph where all edges in the path have multiplicity < pruneFactor
 *
 *
 * For A -[1]> B -[1]> C -[1]> D would be removed with pruneFactor 2
 * but A -[1]> B -[2]> C -[1]> D would not be because the linear chain includes an edge with weight >= 2
 *
 */
public final class LowWeightChainPruner<V extends BaseVertex, E extends BaseEdge> extends ChainPruner<V,E> {
    private final int pruneFactor;

    public LowWeightChainPruner(final int pruneFactor) {
        Utils.validateArg( pruneFactor >= 0, "pruneFactor must be >= 0 but got " + pruneFactor);
        this.pruneFactor = pruneFactor;
    }

    @Override
    protected Collection<Path<V,E>> chainsToRemove(final List<Path<V, E>> chains) {
        return chains.stream().filter(this::needsPruning).collect(Collectors.toList());
    }

    private boolean needsPruning(final Path<V,E> chain) {
        return chain.getEdges().stream().allMatch(e -> e.getPruningMultiplicity() < pruneFactor && !e.isRef());
    }
}
