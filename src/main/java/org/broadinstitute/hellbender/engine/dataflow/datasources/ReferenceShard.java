package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.coders.*;
import com.google.cloud.dataflow.sdk.values.KV;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;

/**
 * ReferenceShard is section of the reference genome that's used for sharding work for pairing things with
 * the reference. The primary use case is pairing reads with the reference.
 * This class designed to be simple a data-storage class with related static utilities.
 */
public class ReferenceShard implements Serializable {
    private static final long serialVersionUID = 1L;

    private final int shardNumber; // shardNumber is zero-based.
    private final String contig;

    public static int REFERENCE_SHARD_SIZE = 100000; // This value is likely to be tweaked.

    public ReferenceShard(int shardNumber, String contig) {
        this.shardNumber = shardNumber;
        Utils.nonNull(contig);
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

    @Override
    public String toString() {
        return "ReferenceShard{" +
                "shardNumber=" + shardNumber +
                ", contig='" + contig + '\'' +
                '}';
    }

    /**
     * getShardNumberFromInterval determines which shard the location maps to (using the start position)
     * @param location, the start of which is used to determine the shard
     * @return the shard (contig + id)
     */
    static public ReferenceShard getShardNumberFromInterval(final Locatable location) {
        return new ReferenceShard(location.getStart()/REFERENCE_SHARD_SIZE, location.getContig());
    }

    public static final DelegateCoder<ReferenceShard, KV<Integer, String>> CODER =
            DelegateCoder.of(
                    KvCoder.of(VarIntCoder.of(), StringUtf8Coder.of()),
                    new DelegateCoder.CodingFunction<ReferenceShard, KV<Integer, String>>() {
                        private static final long serialVersionUID = 1L;
                        @Override
                        public KV<Integer, String> apply(ReferenceShard ref) throws Exception {
                            return KV.of(ref.getShardNumber(), ref.getContig());
                        }
                    },
                    new DelegateCoder.CodingFunction<KV<Integer, String>, ReferenceShard>() {
                        private static final long serialVersionUID = 1L;
                        @Override
                        public ReferenceShard apply(KV<Integer, String> kv) throws Exception {
                            return new ReferenceShard(kv.getKey(), kv.getValue());
                        }
                    }
            );
}
