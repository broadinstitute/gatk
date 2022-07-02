package org.broadinstitute.hellbender.tools.longreads.graph;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;

import java.util.Arrays;
import java.util.Objects;

public class AlignedBaseVertex extends SeqVertex implements Comparable<AlignedBaseVertex> {

    private static final long serialVersionUID = 0x1337L;

    private String readName;
    private GenomicAndInsertionPosition pos;

    public AlignedBaseVertex(final byte[] bases,
                             final String contig,
                             final int genomePos,
                             final int insertionOffset,
                             final String readName) {
        super(bases);
        this.pos = new GenomicAndInsertionPosition(contig, genomePos, insertionOffset);
        this.readName = readName;
    }

    public AlignedBaseVertex(final byte[] bases,
                             final GenomicAndInsertionPosition pos,
                             final String readName) {
        super(bases);
        this.pos = pos;
        this.readName = readName;
    }


    public String getReadName() { return readName; }
    public GenomicAndInsertionPosition getPos() { return pos; }
    public boolean isSingleBase() { return getSequence().length == 1; }

    /**
     * Check to see if the given sequence and the sequence of {@code this} {@link AlignedBaseVertex} are the same.
     * @param seq The {@code byte[]} containing a sequence to check for equality to the sequence in {@code this} {@link AlignedBaseVertex}.
     * @return {@code true} iff the given sequence and the sequence of {@code this} {@link AlignedBaseVertex} are the same.
     */
    public boolean isSequenceEqual(final byte[] seq ) {
        return Arrays.equals(getSequence(), seq);
    }

    /**
     * Check to see if the sequences underneath {@code this} {@link AlignedBaseVertex} and {@code that} {@link AlignedBaseVertex} are the same.
     * @param that The {@link AlignedBaseVertex} to check for sequence equality.
     * @return {@code true} iff the sequence of  {@code this} {@link AlignedBaseVertex} and {@code that} {@link AlignedBaseVertex} are the same.
     */
    public boolean isSequenceEqual(final AlignedBaseVertex that) {
        return isSequenceEqual(that.getSequence());
    }

    /**
     * Checks for adjacency between two {@link AlignedBaseVertex}s.
     * @param that The other {@link AlignedBaseVertex} to check for adjacency.
     * @return {@code true} iff that {@link AlignedBaseVertex} is adjacent to this {@link AlignedBaseVertex}.
     */
    public boolean isAdjacentTo(final AlignedBaseVertex that) {
        return this.getPos().isAdjacentTo(that.getPos());
    }

    /**
     * Get the unique ID for this SeqVertex
     * @return a positive integer >= 0
     */
    public int getId() {
        return System.identityHashCode(this);
    }

    @Override
    public String getGexfAttributesString() {
        return "<attvalue for=\"0\" value=\"" + getSequenceString().length() + "\" />" +
               "<attvalue for=\"1\" value=\"" + getReadName() + "\" />" +
               "<attvalue for=\"2\" value=\"" + getPos().getStart() + "." + getPos().getInsertionOffset() + "\" />";
    }

    @Override
    public String toString() {
        return "AlignedBaseVertex_id_" + System.identityHashCode(this) + "_seq_" + getSequenceString() + "_read_" + getReadName();
    }

    @Override
    public int compareTo(final AlignedBaseVertex o) {

        // 1 - position
        // 2 - insertion offset
        // 3 - read name

        // Check contig / pos / insertion pos:
        int result = this.getPos().compareTo(o.getPos());
        if ( result == 0 ) {
            // compare read name:
            result = this.getReadName().compareTo(o.getReadName());
        }

        return result;
    }

    @Override
    public boolean equals(final Object o) {
        if ( this == o ) return true;
        if ( o == null || getClass() != o.getClass() ) return false;
        if ( !super.equals(o) ) return false;

        final AlignedBaseVertex vertex = (AlignedBaseVertex) o;

        if ( !Objects.equals(readName, vertex.readName) ) return false;
        return Objects.equals(pos, vertex.pos);
    }

    // NOTE: You cannot overwrite the hash code function!
    //       It breaks the graph code.
}
