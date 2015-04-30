package org.broadinstitute.hellbender.engine.dataflow.datasources;

import htsjdk.samtools.util.Locatable;

/**
 * Created by davidada on 5/15/15.
 */
public final class ReferenceShard {
    private int shardNumber;
    private String contig;

    public ReferenceShard(int shardNumber, String contig) {
        this.shardNumber = shardNumber;
        this.contig = contig;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ReferenceShard that = (ReferenceShard) o;

        if (getShardNumber() != that.getShardNumber()) return false;
        return getContig().equals(that.getContig());

    }

    @Override
    public int hashCode() {
        int result = getShardNumber();
        result = 31 * result + getContig().hashCode();
        return result;
    }

    public String getContig() {
        return contig;
    }

    public int getShardNumber() {

        return shardNumber;
    }

    static public ReferenceShard getShardNumberFromInterval(final Locatable location) {
        final int referenceShardSize = 100000;
        return new ReferenceShard(location.getStart()/referenceShardSize, location.getContig());
    }

}
