package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.util.List;

/**
 * Represents a result from a K-best haplotype search.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class KBestHaplotype<V extends BaseVertex, E extends BaseEdge> extends Path<V, E>{
    private double score;
    private boolean isReference;

    public double score() { return score; }
    public boolean isReference() { return isReference; }

    public KBestHaplotype(final V initialVertex, final BaseGraph<V,E> graph) {
        super(initialVertex, graph);
        score = 0;
    }

    public KBestHaplotype(final KBestHaplotype<V, E> p, final E edge, final int totalOutgoingMultiplicity) {
        super(p, edge);
        score = p.score + computeLogPenaltyScore( edge.getMultiplicity(), totalOutgoingMultiplicity);
        isReference &= edge.isRef();
    }

    public static double computeLogPenaltyScore(int edgeMultiplicity, int totalOutgoingMultiplicity) {
        return Math.log10(edgeMultiplicity) - Math.log10(totalOutgoingMultiplicity);
    }

    public KBestHaplotype(final KBestHaplotype<V, E> p, final List<E> edgesToExtend, final double edgePenalty) {
        super(p, edgesToExtend);
        score = p.score() + edgePenalty;
        isReference &= edgesToExtend.get(edgesToExtend.size() - 1).isRef();
    }

    public final Haplotype haplotype() {
        final Haplotype haplotype = new Haplotype(getBases(),isReference());
        haplotype.setScore(score());
        return haplotype;
    }
}
