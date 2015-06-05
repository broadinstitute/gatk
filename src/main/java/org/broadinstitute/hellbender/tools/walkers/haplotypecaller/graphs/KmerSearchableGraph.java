package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;

/**
 * Common interface for those graphs that implement vertex by kmer look-up.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public interface KmerSearchableGraph<V extends BaseVertex, E extends BaseEdge> {

    /**
     * Returns the vertex that represents or contains the last base of a given kmer.
     * @param k the query kmer.
     *
     * @throws NullPointerException if {@code k} is {@code null}.
     * @return {@code null} if there is no such a kmer in the graph or it is not unique.
     */
    V findKmer(Kmer k);

    /**
     * The kmer-size of indexed kmers.
     *
     * @return greater than 0.
     */
    int getKmerSize();

}
