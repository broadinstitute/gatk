package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseVertex;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.LinkedList;
import java.util.List;

/**
 * A DeBruijnVertex that supports multiple copies of the same kmer
 */
public final class MultiDeBruijnVertex extends BaseVertex {
    private static final byte[][] sufficesAsByteArray = new byte[Byte.MAX_VALUE - Byte.MIN_VALUE + 1][];
    static {
        for ( int i = 0; i < sufficesAsByteArray.length; i++ ) {
            sufficesAsByteArray[i] = new byte[]{(byte) (i & 0xFF)};
        }
    }

    private static final boolean KEEP_TRACK_OF_READS = false;

    private final List<String> reads = new LinkedList<>();
    private final boolean mergeIdenticalNodes;
    private final int hashCode;

    /**
     * Create a new MultiDeBruijnVertex with kmer sequence
     * @param mergeIdenticalNodes should nodes with the same sequence be treated as equal?
     * @param sequence the kmer sequence
     */
    public MultiDeBruijnVertex(final byte[] sequence, final boolean mergeIdenticalNodes) {
        super(sequence);
        this.mergeIdenticalNodes = mergeIdenticalNodes;

        //cache hashcode because its computation shows on profiler
        this.hashCode = mergeIdenticalNodes ? super.hashCode() : System.identityHashCode(this);
    }

    /**
     * Create a new MultiDeBruijnVertex with kmer sequence
     * @param sequence the kmer sequence
     */
    public MultiDeBruijnVertex(final byte[] sequence) {
        this(sequence, false);
    }

    @Override
    public boolean equals(final Object o) {
        if (mergeIdenticalNodes) {
            return super.equals(o);
        } else {
            return o == this;
        }
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    @Override
    public String toString() {
        return "MultiDeBruijnVertex_id_" + hashCode() + "_seq_" + getSequenceString();
    }

    /**
     * Add name information to this vertex for debugging
     *
     * This information will be captured as a list of strings, and displayed in DOT if this
     * graph is written out to disk
     *
     * This functionality is only enabled when KEEP_TRACK_OF_READS is true
     *
     * @param name a non-null string
     */
    public void addRead(final String name) {
        Utils.nonNull(name, "name cannot be null");
        if ( KEEP_TRACK_OF_READS ) {
            reads.add(name);
        }
    }

    @Override
    public String getAdditionalInfo() {
        if (reads.contains("ref")) {
            return super.getAdditionalInfo();
        } else {
            return super.getAdditionalInfo() + (KEEP_TRACK_OF_READS ? "__" + Utils.join(",", reads) : "");
        }
    }

    /**
     * Get the kmer size for this DeBruijnVertex
     * @return integer >= 1
     */
    public int getKmerSize() {
        return sequence.length;
    }

    /**
     * Get the string representation of the suffix of this DeBruijnVertex
     * @return a non-null non-empty string
     */
    public String getSuffixString() {
        return new String(getSuffixAsArray());
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
