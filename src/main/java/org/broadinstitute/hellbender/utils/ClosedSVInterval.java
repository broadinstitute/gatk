package org.broadinstitute.hellbender.utils;

public final class ClosedSVInterval extends SVInterval {
    private final int contig;
    private final int start;
    private final int end;
    public ClosedSVInterval(int contig, int start, int end) {
        super(contig, start, end);
        this.contig = contig;
        this.start = start;
        this.end = end;
    }

    @Override
    public SimpleInterval toSimpleInterval(String[] contigIDToName) {
        Integer contigID = this.getContig();  // non-negative
        if (contigID >= contigIDToName.length) {
            throw new ArrayIndexOutOfBoundsException("Contig ID " + contigID + " out of bounds of provided contig ID to name map");
        }
        return new SimpleInterval(contigIDToName[contigID], this.getStart(), this.getEnd());
    }

    @Override
    public int getLength() { return end - start + 1; }

    // The following functions are not technically overrides because comparison input "that" is different class
    public boolean overlaps( final ClosedSVInterval that ) {
        return this.contig == that.contig && this.start <= that.end && that.start <= this.end;
    }

    public boolean isDisjointFrom( final ClosedSVInterval that ) {
        return !overlaps(that);
    }

    /**
     * Assumes {@link #isUpstreamOf(ClosedSVInterval)}.
     */
    public int gapLen( final ClosedSVInterval that ) {
        if ( this.contig != that.contig ) return Integer.MAX_VALUE;
        return that.start - this.end - 1;
    }

    public int overlapLen( final ClosedSVInterval that ) {
        if ( this.isDisjointFrom(that)) return 0;
        return Math.max(0, Math.min(this.end, that.end) - Math.max(this.start, that.start)) + 1;
    }

    public boolean isUpstreamOf( final ClosedSVInterval that ) {
        return this.contig < that.contig || (this.contig == that.contig && this.end < that.start);
    }
}
