package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Kmer Hopscotch set class that encapsulates the filter, kmer size, and kmer mask
 */
@DefaultSerializer(PSKmerSet.Serializer.class)
public final class PSKmerSet extends PSKmerCollection {

    private final LargeLongHopscotchSet kmerSet;
    private final int kmerSize;
    private final SVKmerShort kmerMask;

    /**
     * Note values in the input set should have been run through PSKmerCollection's canonicalizeAndMask()
     */
    public PSKmerSet(final LargeLongHopscotchSet maskedKmerSet, final int kmerSize, final SVKmerShort kmerMask) {
        Utils.nonNull(maskedKmerSet);
        Utils.nonNull(kmerMask);
        this.kmerSet = maskedKmerSet;
        this.kmerSize = kmerSize;
        this.kmerMask = kmerMask;
    }

    private PSKmerSet(final Kryo kryo, final Input input) {
        this.kmerSize = input.readInt();
        this.kmerMask = new SVKmerShort(input.readLong());
        this.kmerSet = kryo.readObject(input, LargeLongHopscotchSet.class);
    }

    /**
     * Input should not be canonicalized/masked
     */
    @Override
    public boolean contains(final SVKmerShort rawKmer) {
        return kmerSet.contains(canonicalizeAndMask(rawKmer, kmerSize, kmerMask));
    }

    public LongIterator iterator() {
        return kmerSet.iterator();
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
        return 0;
    }

    public long setSize() {
        return kmerSet.size();
    }

    private void serialize(final Kryo kryo, final Output output) {
        output.writeInt(kmerSize);
        output.writeLong(kmerMask.getLong());
        kryo.writeObject(output, kmerSet);
        output.close();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PSKmerSet)) return false;

        final PSKmerSet psKmerSet = (PSKmerSet) o;

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

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<PSKmerSet> {
        @Override
        public void write(final Kryo kryo, final Output output, final PSKmerSet kmerSet) {
            kmerSet.serialize(kryo, output);
        }

        @Override
        public PSKmerSet read(final Kryo kryo, final Input input, final Class<PSKmerSet> klass) {
            return new PSKmerSet(kryo, input);
        }
    }

}
