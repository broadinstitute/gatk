package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;

/**
 * ReferenceShard is section of the reference genome that's used for sharding work for pairing things with
 * the reference. The primary use case is pairing reads with the reference.
 * This class is designed to be a simple data-storage class with related static utilities.
 */
public final class ReferenceShard implements Serializable {
    private static final long serialVersionUID = 1L;

    private final int shardNumber; // shardNumber is zero-based.
    private final String contig;

    public static final int REFERENCE_SHARD_SIZE = 10000; // This value is subject to change (by humans).

    public ReferenceShard(int shardNumber, String contig) {
        this.shardNumber = shardNumber;
        this.contig = Utils.nonNull(contig);
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
        //using a 4 digit prime because we found empirically that 31 ended up with badly dispersed values in spark
        result = 5779 * result + getContig().hashCode();
        return result;
    }

    public String getContig() {
        return contig;
    }

    public int getShardNumber() {
        return shardNumber;
    }

    @Override
    public String toString() {
        return "ReferenceShard{" +
                "shardNumber=" + shardNumber +
                ", contig='" + contig + '\'' +
                '}';
    }

    /**
     * getShardNumberFromInterval returns the ReferenceShard that overlap the read's start position.
     * @param location, the start of which is used to determine the shard
     * @return the shard (contig + id)
     */
    public static ReferenceShard getShardNumberFromInterval(final Locatable location) {
        return new ReferenceShard(location.getStart()/REFERENCE_SHARD_SIZE, location.getContig());
    }

}
