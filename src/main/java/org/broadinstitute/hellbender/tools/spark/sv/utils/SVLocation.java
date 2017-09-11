package org.broadinstitute.hellbender.tools.spark.sv.utils;

public final class SVLocation implements Comparable<SVLocation> {
    private int contig;
    private int position;

    public SVLocation( final int contig, final int position ) {
        this.contig = contig;
        this.position = position;
    }

    public int getContig() { return contig; }
    public int getPosition() { return position; }

    @Override
    public boolean equals( final Object obj ) {
        if ( this == obj ) return true;
        if ( obj == null || obj.getClass() != SVLocation.class ) return false;
        SVLocation that = (SVLocation)obj;
        return contig == that.contig && position == that.position;
    }

    @Override
    public int hashCode() {
        final int mult = 16777619;
        return mult*(mult*(0x811c9dc5 ^ contig) ^ position);
    }

    @Override
    public int compareTo( final SVLocation that ) {
        int result = Integer.compare(contig, that.contig);
        if ( result == 0 ) result = Integer.compare(position, that.position);
        return result;
    }

    @Override
    public String toString() {
        return contig + ":" + position;
    }
}
