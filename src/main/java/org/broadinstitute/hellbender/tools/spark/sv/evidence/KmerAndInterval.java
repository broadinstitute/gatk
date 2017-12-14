package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerLong;

import java.util.Map;

/**
 * A <Kmer,IntervalId> pair.
 */
@DefaultSerializer(KmerAndInterval.Serializer.class)
public final class KmerAndInterval extends SVKmerLong implements Map.Entry<SVKmer, Integer> {
    private final int intervalId;

    public KmerAndInterval( final SVKmer kmer, final int intervalId ) {
        super(kmer);
        this.intervalId = intervalId;
    }

    private KmerAndInterval( final Kryo kryo, final Input input ) {
        super(kryo, input);
        intervalId = input.readInt();
    }

    @Override
    protected void serialize( final Kryo kryo, final Output output ) {
        super.serialize(kryo, output);
        output.writeInt(intervalId);
    }

    @Override
    public SVKmer getKey() { return new SVKmerLong(this); }
    @Override
    public Integer getValue() { return intervalId; }
    @Override
    public Integer setValue( final Integer value ) {
        throw new UnsupportedOperationException("Can't set KmerAndInterval.intervalId");
    }

    @Override
    public boolean equals( final Object obj ) {
        return obj instanceof KmerAndInterval && equals((KmerAndInterval)obj);
    }

    public boolean equals( final KmerAndInterval that ) {
        return super.equals(that) && this.intervalId == that.intervalId;
    }

    public int getIntervalId() { return intervalId; }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<KmerAndInterval> {
        @Override
        public void write( final Kryo kryo, final Output output, final KmerAndInterval kmerAndInterval) {
            kmerAndInterval.serialize(kryo, output);
        }

        @Override
        public KmerAndInterval read( final Kryo kryo, final Input input,
                                     final Class<KmerAndInterval> klass ) {
            return new KmerAndInterval(kryo, input);
        }
    }
}
