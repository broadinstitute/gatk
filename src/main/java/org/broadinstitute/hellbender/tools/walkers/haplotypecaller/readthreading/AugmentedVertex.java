package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseVertex;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.LinkedList;
import java.util.List;

/**
 * A DeBruijnVertex that supports multiple copies of the same kmer
 */
public final class AugmentedVertex extends BaseVertex {
    private static final byte[][] sufficesAsByteArray = new byte[Byte.MAX_VALUE - Byte.MIN_VALUE + 1][];
    static {
        for ( int i = 0; i < sufficesAsByteArray.length; i++ ) {
            sufficesAsByteArray[i] = new byte[]{(byte) (i & 0xFF)};
        }
    }

    private final int hashCode;
    private final int position;


    public AugmentedVertex(final byte[] sequence, final int position) {
        super(sequence);
        this.position = position;

        //cache hashcode because its computation shows on profiler
        hashCode = 7 * super.hashCode() + 13 * position;
    }

    @Override
    public boolean equals(final Object o) {
        return super.equals(o) && ((AugmentedVertex) o).position == this.position;
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    @Override
    public String toString() {
        return "AugmentedVertex" + "_seq_" + getSequenceString() + "_pos_" + position;
    }

    /**
     * Get the kmer size for this DeBruijnVertex
     * @return integer >= 1
     */
    public int getKmerSize() {
        return sequence.length;
    }

    /**
     * Get the suffix byte of this DeBruijnVertex
     *
     * The suffix byte is simply the last byte of the kmer sequence, so if this is holding sequence ACT
     * getSuffix would return T
     *
     * @return a byte
     */
    public byte getSuffix() {
        return sequence[getKmerSize() - 1];
    }

    /**
     * Optimized version that returns a byte[] for the single byte suffix of this graph without allocating memory.
     *
     * Should not be modified
     *
     * @return a byte[] that contains 1 byte == getSuffix()
     */
    private byte[] getSuffixAsArray() {
        return sufficesAsByteArray[getSuffix()];
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public byte[] getAdditionalSequence(final boolean source) {
        return source ? super.getAdditionalSequence(source) : getSuffixAsArray();
    }
}
