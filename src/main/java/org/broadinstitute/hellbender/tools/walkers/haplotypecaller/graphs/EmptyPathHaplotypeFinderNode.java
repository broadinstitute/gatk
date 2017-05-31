package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Collections;
import java.util.Set;

/**
 * Trivial k-best sub-haplotype finder where the source and sink vertex are the same one.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
final class EmptyPathHaplotypeFinderNode implements KBestSubHaplotypeFinder {

    /**
     * Caches the only solution returned by this finder.
     */
    private final KBestHaplotype singleHaplotypePath;

    /**
     * Constructs a new empty k-best haplotype finder.
     *
     * @param graph the search graph.
     * @param vertex the source and sink vertex of the only solution returned by this finder.
     */
    EmptyPathHaplotypeFinderNode(final SeqGraph graph, final SeqVertex vertex) {
        singleHaplotypePath = new MyBestHaplotypePath(graph,vertex);
    }

    @Override
    public String id() {
        return "v" + Integer.valueOf(singleHaplotypePath.head().getId());
    }

    @Override
    public String label() {
        return singleHaplotypePath.head().getSequenceString();
    }

    @Override
    public Set<Pair<? extends KBestSubHaplotypeFinder, String>> subFinderLabels() {
        return Collections.emptySet();
    }

    @Override
    public int getCount() {
        return 1;
    }

    @Override
    public KBestHaplotype getKBest(final int k) {
        ParamUtils.isPositiveOrZero(k, "k cannot be negative");
        Utils.validateArg(k == 0, "k cannot greater than the possible haplotype count");
        return singleHaplotypePath;
    }

    @Override
    public boolean isReference() {
        return singleHaplotypePath.isReference();
    }

    @Override
    public double score(final byte[] bases, final int offset, final int length) {
        Utils.nonNull(bases, "bases cannot be null");
        ParamUtils.isPositiveOrZero(offset, "the offset cannot be negative");
        ParamUtils.isPositiveOrZero(length, "the length cannot be negative");
        Utils.validateArg(offset + length <= bases.length, "the offset and length go beyond the array size");
        final byte[] vertexBases = singleHaplotypePath.head().getSequence();
        if (length != vertexBases.length) {
            return Double.NaN;
        } else {
            return Utils.equalRange(bases, offset, vertexBases, 0, length) ? 0.0 : Double.NaN;
        }
    }

    /**
     * Custom extension of {@link KBestHaplotype} that implements the single solution behaviour.
     */
    private static final class MyBestHaplotypePath extends KBestHaplotype {

        /**
         * The solution's only vertex.
         */
        private final SeqVertex vertex;

        /**
         * The search graph.
         */
        private final SeqGraph graph;

        /**
         * Whether the vertex is a reference vertex.
         */
        private final boolean isReference;

        /**
         * Constructs a new empty k-best haplotype solution.
         *
         * @param graph the search graph.
         * @param vertex the source and sink vertex of the only solution returned by the outer finder.
         */
        MyBestHaplotypePath(final SeqGraph graph, final SeqVertex vertex) {
            this.vertex = vertex;
            this.graph = graph;
            isReference = graph.isReferenceNode(vertex);
        }

        @Override
        public SeqGraph graph() {
            return graph;
        }

        @Override
        public double score() {
            return 0;
        }

        @Override
        public int rank() {
            return 0;
        }

        @Override
        protected SeqVertex head() {
            return vertex;
        }

        @Override
        protected KBestHaplotype tail() {
            return null;
        }

        @Override
        public boolean isReference() {
            return isReference;
        }
    }
}