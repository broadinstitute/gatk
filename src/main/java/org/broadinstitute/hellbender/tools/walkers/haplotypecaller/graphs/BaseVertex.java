package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;

/**
 * A graph vertex that holds some sequence information
 */
public class BaseVertex {
    /** placeholder to store additional information for debugging purposes */
    private String additionalInfo = "";
    protected final byte[] sequence;
    private final int cachedHashCode;

    /**
     * Create a new sequence vertex with sequence
     *
     * This code doesn't copy sequence for efficiency reasons, so sequence must absolutely not be modified
     * in any way after passing this sequence to the BaseVertex
     *
     * @param sequence a non-null sequence of bases contained in this vertex
     */
    public BaseVertex(final byte[] sequence) {
        Utils.nonNull(sequence, "Sequence cannot be null");
        this.sequence = sequence;
        cachedHashCode = Arrays.hashCode(sequence);

    }

    /**
     * For testing purposes only -- low performance
     * @param sequence the sequence as a string
     */
    @VisibleForTesting
    BaseVertex(final String sequence) {
        this(sequence.getBytes());
    }


    /**
     * Does this vertex have an empty sequence?
     *
     * That is, is it a dummy node that's only present for structural reasons but doesn't actually
     * contribute to the sequence of the graph?
     *
     * @return true if sequence is empty, false otherwise
     */
    public final boolean isEmpty() {
        return length() == 0;
    }

    /**
     * Get the length of this sequence
     * @return a positive integer >= 1
     */
    public final int length() {
        return sequence.length;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final BaseVertex that = (BaseVertex) o;
        if (hashCode() != that.hashCode()){
            return false;
        }

        return seqEquals(that);
    }

    /**
     * necessary to override here so that graph.containsVertex() works the same way as vertex.equals() as one might expect
     * @return
     */
    @Override
    public int hashCode() {
        return cachedHashCode;
    }

    /**
     * Are b and this equal according to their base sequences?
     *
     * @param b the vertex to compare ourselves to
     * @return true if b and this have the same sequence, regardless of other attributes that might differentiate them
     */
    public final boolean seqEquals(final BaseVertex b) {
        return Arrays.equals(getSequence(), b.getSequence());
    }

    @Override
    public String toString() {
        return getSequenceString();
    }

    /**
     * Get the sequence of bases contained in this vertex
     *
     * Do not modify these bytes in any way!
     *
     * @return a non-null pointer to the bases contained in this vertex
     */
    public final byte[] getSequence() {
        return sequence;
    }

    /**
     * Get a string representation of the bases in this vertex
     * @return a non-null String
     */
    public final String getSequenceString() {
        return new String(sequence);
    }

    /**
     * Get the sequence unique to this vertex
     *
     * This function may not return the entire sequence stored in the vertex, as kmer graphs
     * really only provide 1 base of additional sequence (the last base of the kmer).
     *
     * The base implementation simply returns the sequence.
     *
     * @param source is this vertex a source vertex (i.e., no in nodes) in the graph
     * @return a byte[] of the sequence added by this vertex to the overall sequence
     */
    public byte[] getAdditionalSequence(final boolean source) {
        return getSequence();
    }

    /**
     * Set additional debugging information for this vertex
     * @param info the new info value.
     */
    public final void setAdditionalInfo(final String info) {
        Utils.nonNull(info, "info cannot be null");
        additionalInfo = info;
    }

    /**
     * @return the additional information for display about this vertex
     */
    public String getAdditionalInfo() { return additionalInfo; }

    /**
     * Checks whether the vertex sequence is ambiguous or not.
     *
     * <p>
     *     Ambiguity may come about as a result of either:
     *     <ul>
     *        <li>by construction as the generating sequence (read or haplotype) had ambiguous bases</li>
     *        <li>or because this vertex is the result of merging two or more vertices with some variation upstream
     *        no more than kmerSize bases away.</li>
     *     </ul>
     * </p>
     *
     * @return {@code true} iff so.
     */
    public final boolean hasAmbiguousSequence() {
        for (final byte base : sequence) {
            switch (Character.toUpperCase(base)) {
                case 'A':
                case 'T':
                case 'G':
                case 'C':
                    continue;
                default:
                    return true;
            }
        }
        return false;
    }
}
