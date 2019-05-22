package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

/**
 * Represents a result from a K-best haplotype search.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class KBestHaplotype<V extends BaseVertex, E extends BaseEdge> extends Path<V, E>{
    private double score;
    private boolean isReference;

    public double score() { return score; }
    public boolean isReference() { return isReference; }

    public KBestHaplotype(final V initialVertex, final BaseGraph<V,E> graph) {
        super(initialVertex, graph);
        score = 0;
    }

    public KBestHaplotype(final KBestHaplotype p, final E edge, final int totalOutgoingMultiplicity) {
        super(p, edge);
        score = p.score() + MathUtils.log10(edge.getMultiplicity()) - MathUtils.log10(totalOutgoingMultiplicity);
        isReference &= edge.isRef();
    }

    public final Haplotype haplotype() {
        final Haplotype haplotype = new Haplotype(getBases(),isReference());
        haplotype.setScore(score());
        return haplotype;
    }
}
