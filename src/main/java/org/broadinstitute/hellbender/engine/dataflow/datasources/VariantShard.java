package org.broadinstitute.hellbender.engine.dataflow.datasources;

import htsjdk.samtools.util.Locatable;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by davidada on 5/15/15.
 */
public final class VariantShard {
    private int shardNumber;
    private String contig;

    public VariantShard(int shardNumber, String contig) {
        this.shardNumber = shardNumber;
        this.contig = contig;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        VariantShard that = (VariantShard) o;

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

    /*
      * May have some bugs...
      */
    static public List<VariantShard> getVariantShardsFromInterval(final Locatable location) {
        final int variantShardSize = 100000;
        List<VariantShard> intervalList = new ArrayList<>();
        // Get all of the shard numbers that span the start and end of the interval.
        int startShard = location.getStart()/variantShardSize;
        int endShard = location.getEnd()/variantShardSize;
        for (int i = startShard; i <= endShard; ++i) {
            intervalList.add(new VariantShard(i, location.getContig()));
        }
        return intervalList;
    }
}
