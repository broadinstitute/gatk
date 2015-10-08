package org.broadinstitute.hellbender.engine;

import com.google.cloud.dataflow.sdk.coders.DelegateCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.coders.StringUtf8Coder;
import com.google.cloud.dataflow.sdk.coders.VarIntCoder;
import com.google.cloud.dataflow.sdk.values.KV;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;

/**
 * VariantShard is section of the genome that's used for sharding work for pairing things with
 * variants. The primary use case is pairing reads with overlapping variants.
 * This is class designed to be simple a data-storage class with related static utilities.
 */
public final class VariantShard {
    private final int shardNumber;
    private final String contig;
    public static final int VARIANT_SHARDSIZE = 1000; // This value is subject to change (by humans)

    public VariantShard(int shardNumber, String contig) {
        this.shardNumber = shardNumber;
        this.contig = Utils.nonNull(contig);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        VariantShard that = (VariantShard) o;

        return getShardNumber() == that.getShardNumber() && getContig().equals(that.getContig());

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

    /**
     * getVariantShardsFromInterval calculates *all* shards overlapping location.
     * @param location the range of sites determines which shards are overlapping
     * @return All overlapping VariantShards
     */
    public static List<VariantShard> getVariantShardsFromInterval(final Locatable location) {
        if (location.getContig()==null) {
            // don't feed me unmapped reads!
            throw new GATKException("getVariantShardsFromInterval requires locations to be mapped");
        }
        List<VariantShard> shardList = new ArrayList<>();
        // Get all of the shard numbers that span the start and end of the interval.
        int startShard = location.getStart()/ VARIANT_SHARDSIZE;
        int endShard = location.getEnd()/ VARIANT_SHARDSIZE;
        for (int i = startShard; i <= endShard; ++i) {
            shardList.add(new VariantShard(i, location.getContig()));
        }
        return shardList;
    }

    @Override
    public String toString() {
        return "VariantShard{" +
                "shardNumber=" + shardNumber +
                ", contig='" + contig + '\'' +
                '}';
    }

    public static final DelegateCoder<VariantShard, KV<Integer, String>> CODER =
            DelegateCoder.of(
                    KvCoder.of(VarIntCoder.of(), StringUtf8Coder.of()),
                    new DelegateCoder.CodingFunction<VariantShard, KV<Integer, String>>() {
                        private static final long serialVersionUID = 1L;
                        @Override
                        public KV<Integer, String> apply(VariantShard ref) throws Exception {
                            return KV.of(ref.getShardNumber(), ref.getContig());
                        }
                    },
                    new DelegateCoder.CodingFunction<KV<Integer, String>, VariantShard>() {
                        private static final long serialVersionUID = 1L;
                        @Override
                        public VariantShard apply(KV<Integer, String> kv) throws Exception {
                            return new VariantShard(kv.getKey(), kv.getValue());
                        }
                    }
            );
}
