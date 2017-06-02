package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import java.util.Arrays;

/**
 * A graph vertex containing a sequence of bases and a unique ID that
 * allows multiple distinct nodes in the graph to have the same sequence.
 *
 * This is essential when thinking about representing the actual sequence of a haplotype
 * in a graph.  There can be many parts of the sequence that have the same sequence, but
 * are distinct elements in the graph because they have a different position in the graph.  For example:
 *
 * A -> C -> G -> A -> T
 *
 * The two As are not the same, because they occur with different connections.  In a kmer graph equals()
 * is based on the sequence itself, as each distinct kmer can only be represented once.  But the transformation
 * of the kmer graph into a graph of base sequences, without their kmer prefixes, means that nodes that
 * where once unique including their prefix can become equal after shedding the prefix.  So we need to
 * use some mechanism -- here a unique ID per node -- to separate nodes that have the same sequence
 * but are distinct elements of the graph.
 *
 */
public final class SeqVertex extends BaseVertex {
    /**
     * Create a new SeqVertex with sequence and the next available id
     * @param sequence our base sequence
     */
    public SeqVertex(final byte[] sequence) {
        super(sequence);
    }

    /**
     * Create a new SeqVertex having bases of sequence.getBytes()
     * @param sequence the string representation of our bases
     */
    public SeqVertex(final String sequence) {
        super(sequence);
    }

    /**
     * Get the unique ID for this SeqVertex
     * @return a positive integer >= 0
     */
    public int getId() {
        return hashCode();
    }

    @Override
    public String toString() {
        return "SeqVertex_id_" + hashCode() + "_seq_" + getSequenceString();
    }

    /**
     * Two SeqVertex are equal only if their ids are equal
     * @param o
     * @return
     */
    @Override
    public boolean equals(final Object o) { return o == this; }

    @Override
    public int hashCode() {
        return System.identityHashCode(this);
    }

    /**
     * Return a new SeqVertex derived from this one but not including the suffix bases
     *
     * @param suffix the suffix bases to remove from this vertex
     * @return a newly allocated SeqVertex with appropriate prefix, or null if suffix removes all bases from this node
     */
    public SeqVertex withoutSuffix(final byte[] suffix) {
        final int prefixSize = sequence.length - suffix.length;
        return prefixSize > 0 ? new SeqVertex(Arrays.copyOf(sequence, prefixSize)) : null;
    }

    /**
     * Return a new SeqVertex derived from this one but not including prefix or suffix bases
     *
     * @param prefix the previx bases to remove
     * @param suffix the suffix bases to remove from this vertex
     * @return a newly allocated SeqVertex
     */
    public SeqVertex withoutPrefixAndSuffix(final byte[] prefix, final byte[] suffix) {
        final int start = prefix.length;
        final int length = sequence.length - suffix.length - prefix.length;
        final int stop = start + length;
        return length > 0 ? new SeqVertex(Arrays.copyOfRange(sequence, start, stop)) : null;
    }
}
