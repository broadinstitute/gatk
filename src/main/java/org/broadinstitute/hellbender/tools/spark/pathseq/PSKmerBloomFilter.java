package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.utils.LongBloomFilter;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Kmer Bloom Filter class that encapsulates the filter, kmer size, and kmer mask
 */
@DefaultSerializer(PSKmerBloomFilter.Serializer.class)
public final class PSKmerBloomFilter extends PSKmerCollection {

    private final LongBloomFilter kmerSet;
    private final int kmerSize;
    private final SVKmerShort kmerMask;
    private final double falsePositiveProbability;

    /**
     * Note values in the Bloom filter should have been run through PSKmerCollection's canonicalizeAndMask()
     */
    public PSKmerBloomFilter(final LongBloomFilter maskedKmerBloomFilter, final int kmerSize, final SVKmerShort kmerMask,
                             final long numElements) {
        Utils.nonNull(maskedKmerBloomFilter);
        Utils.nonNull(kmerMask);
        this.kmerSet = maskedKmerBloomFilter;
        this.kmerSize = kmerSize;
        this.kmerMask = kmerMask;
        this.falsePositiveProbability = kmerSet.getTheoreticalFPP(numElements);
    }

    private PSKmerBloomFilter(final Kryo kryo, final Input input) {
        this.kmerSize = input.readInt();
        this.kmerMask = new SVKmerShort(input.readLong());
        this.kmerSet = kryo.readObject(input, LongBloomFilter.class);
        this.falsePositiveProbability = input.readDouble();
    }

    /**
     * Input should not be canonicalized/masked
     */
    @Override
    public boolean contains(final SVKmerShort rawKmer) {
        return kmerSet.contains(canonicalizeAndMask(rawKmer, kmerSize, kmerMask));
    }

    @Override
    public int kmerSize() {
        return kmerSize;
    }

    @Override
    public SVKmerShort getMask() {
        return kmerMask;
    }

    @Override
    public double getFalsePositiveProbability() {
        return falsePositiveProbability;
    }

    private void serialize(final Kryo kryo, final Output output) {
        output.writeInt(kmerSize);
        output.writeLong(kmerMask.getLong());
        kryo.writeObject(output, kmerSet);
        output.writeDouble(falsePositiveProbability);
        output.close();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PSKmerBloomFilter)) return false;

        final PSKmerBloomFilter psKmerSet = (PSKmerBloomFilter) o;

        if (kmerSize != psKmerSet.kmerSize) return false;
        if (!kmerSet.equals(psKmerSet.kmerSet)) return false;
        return kmerMask.equals(psKmerSet.kmerMask);
    }

    @Override
    public int hashCode() {
        int result = kmerSet.hashCode();
        result = 31 * result + kmerSize;
        result = 31 * result + kmerMask.hashCode();
        return result;
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<PSKmerBloomFilter> {
        @Override
        public void write(final Kryo kryo, final Output output, final PSKmerBloomFilter kmerSet) {
            kmerSet.serialize(kryo, output);
        }

        @Override
        public PSKmerBloomFilter read(final Kryo kryo, final Input input, final Class<PSKmerBloomFilter> klass) {
            return new PSKmerBloomFilter(kryo, input);
        }
    }

}
