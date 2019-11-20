package org.broadinstitute.hellbender.tools.longreads.graph;

import java.util.Objects;

public class GenomicAndInsertionPosition implements Comparable<GenomicAndInsertionPosition> {

    private String contig;
    private int start;
    private int insertionOffset;

    public GenomicAndInsertionPosition(final String contig, final int start, final int insertionOffset) {
        this.contig = contig;
        this.start = start;
        this.insertionOffset = insertionOffset;
    }

    public String getContig() {
        return contig;
    }

    public int getStart() {
        return start;
    }

    public int getInsertionOffset() {
        return insertionOffset;
    }

    /**
     * Checks for adjacency between two {@link GenomicAndInsertionPosition}s.
     * Adjacent positions are those that are next to eachother in genomic coordinates, including those with equal
     * start positions, but insertion offsets that differ by 1.
     * A {@link GenomicAndInsertionPosition} is NEVER adjacent to itself.
     * @param that The other {@link GenomicAndInsertionPosition} to check for adjacency.
     * @return {@code true} iff that {@link GenomicAndInsertionPosition} is adjacent to this {@link GenomicAndInsertionPosition}.
     */
    public boolean isAdjacentTo( final GenomicAndInsertionPosition that ) {
        boolean isAdjacent = false;

        if ( !this.equals(that) ) {
            if ( getContig().equals(that.getContig()) ) {
                final int posOffset = Math.abs( this.getStart() - that.getStart() ) +
                        Math.abs( this.getInsertionOffset() - that.getInsertionOffset() );

                isAdjacent = posOffset == 1;
            }
        }

        return isAdjacent;
    }

    @Override
    public String toString() {
        if ( insertionOffset != 0 ) {
            return String.format("GenomicAndInsertionPosition(%s:%d+%d)", contig, start, insertionOffset);
        }
        else {
            return String.format("GenomicAndInsertionPosition(%s:%d)", contig, start);
        }
    }

    @Override
    public int compareTo(final GenomicAndInsertionPosition o) {
        int result = this.getContig().compareTo(o.getContig());
        if ( result == 0 ) {
            // Check start pos:
            result = Integer.compare(this.getStart(), o.getStart());
            if (result == 0 ) {
                // Check insertion offset:
                result = Integer.compare(this.getInsertionOffset(), o.getInsertionOffset());
            }
        }
        return result;
    }

    @Override
    public boolean equals(final Object o) {
        if ( this == o ) return true;
        if ( o == null || getClass() != o.getClass() ) return false;

        final GenomicAndInsertionPosition that = (GenomicAndInsertionPosition) o;

        if ( start != that.start ) return false;
        if ( insertionOffset != that.insertionOffset ) return false;
        return Objects.equals(contig, that.contig);
    }

    @Override
    public int hashCode() {
        int result = contig != null ? contig.hashCode() : 0;
        result = 31 * result + start;
        result = 31 * result + insertionOffset;
        return result;
    }
}
