package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * An immutable SVKmerShort. This class is the same as SVKmerLong but uses a 1 long instead of 2 longs to store the kmer.
 * K must be between 1 and 31
 * Canonicalization is unimplemented for even K.
 */
@DefaultSerializer(SVKmerShort.Serializer.class)
public class SVKmerShort extends SVKmer implements Comparable<SVKmerShort> {
    // these are treated as K-bit unsigned integers
    private final long valLow; // Kmer bits

    /**
     * Makes an empty SVKmerShort.  If you call toString on it, it'll look like poly-A.
     */
    public SVKmerShort() {
        valLow = 0;
    }

    /**
     *  Makes an empty SVKmerShort.  If you call toString on it, it'll look like poly-A.
     */
    public SVKmerShort( final int kSize ) {
        Utils.validateArg(kSize >= 1 && kSize < 32, "Kmer length must be between 1 and 31.");
        valLow = 0;
    }

    public SVKmerShort( final SVKmerShort that ) { this.valLow = that.valLow; }

    public SVKmerShort( final long valLow ) { this.valLow = valLow; }

    private SVKmerShort(final Kryo kryo, final Input input) {
        valLow = input.readLong();
    }

    private void serialize(final Kryo kryo, final Output output) {
        output.writeLong(valLow);
    }

    public long getLong() {
        return valLow;
    }

    /**
     * Returns a new SVKmerShort that's like this one, but with its leading base discarded and a new one added to the end.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.successor(SVKmerShort.Base.C,5).toString(5) is "CTGAC".
     *
     * @param base must be 0, 1, 2, or 3, corresponding to A, C, G, or T.
     */
    public final SVKmerShort successor( final Base base, final int kSize ) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        // note we multiply kSize by two in SVKmerShort because we no longer divide the bits into two longs
        final long mask = (1L << kSize * 2) - 1L;
        // move all the bits up two places, OR in the pair of successor bits at the bottom, and mask to kSize bits
        final long newV2 = ((valLow << 2) | (base.value & 3L)) & mask;
        return new SVKmerShort(newV2);
    }

    /**
     * Returns a new SVKmerShort that's like this one, but with its trailing base discarded and a new one added to the start.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.predecessor(SVKmerShort.Base.T,5).toString(5) is "TACTG".
     *
     * @param base must be 0, 1, 2, or 3, corresponding to A, C, G, or T.
     */
    public final SVKmerShort predecessor( final Base base, final int kSize ) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        // note we multiply kSize by two in SVKmerShort because we no longer divide the bits into two longs
        final long mask = (1L << kSize * 2) - 1L;
        // move all the bits down two places and mask to kSize bits
        final long newV2 = ((valLow >> 2) | (base.value << (kSize * 2 - 2))) & mask;
        return new SVKmerShort(newV2);
    }

    /**
     * Returns a new SVKmerShort that's the reverse-complement of this one.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.rc(5).toString(5) is "TCAGT".
     */
    public final SVKmerShort reverseComplement(final int kSize) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        // note we multiply kSize by two in SVKmerShort because we no longer divide the bits into two longs
        final long mask = (1L << kSize * 2) - 1L;
        // number of unused bits at the top
        final int compK = 64 - kSize * 2;
        // move the significant bits up to the top, reverse complement, and mask to kSize bits.
        return new SVKmerShort(reverseComplement(valLow << compK) & mask);
    }

    /**
     * Returns a SVKmerShort that is a canonical representation of this one.
     * An odd-K SVKmerShort is in canonical form if its middle base is A or C.
     * The reverse-complement of a non-canonical SVKmerShort is a canonical SVKmerShort, and vice versa.  (Think about it.)
     * Canonical form is not defined for even-K Kmers (too expensive to compute routinely).
     */
    public SVKmerShort canonical( final int kSize ) {
        Utils.validateArg( (kSize & 1) != 0, "Kmer length must be odd to canonicalize.");
        // test middle base's value by ANDing with 1.  if it's zero the middle base is A or C and we're good to go.
        if (((valLow >> kSize) & 1L) == 0) return this;
        // middle base is G or T.  reverse complement.
        return reverseComplement(kSize);
    }

    public final Base firstBase( final int kSize ) { return Base.values()[(int)(valLow >> (kSize*2-2))]; }
    public final Base lastBase() { return Base.values()[(int)(valLow & 3)]; }
    public final int firstTrimer(final int kSize ) { return (int)(valLow >>> (kSize*2-6)); }
    public final int lastTrimer() { return (int)valLow & 0x3F; }

    @Override
    public boolean equals(final Object obj) {
        return obj instanceof SVKmerShort && equals((SVKmerShort) obj);
    }

    public final boolean equals(final SVKmerShort that) {
        return this.valLow == that.valLow;
    }

    @Override
    public final int hashCode() {
        // 32-bit FNV-1a algorithm
        return (int) SVUtils.fnvLong64(2166136261L, valLow);
    }

    /**
     * SVKmerShort comparison is consistent with equals.
     * It's also the same as the lexicographic ordering you'd get using toString on the Kmers.
     */
    @Override
    public final int compareTo(final SVKmerShort that) {
        return Long.compare(this.valLow, that.valLow);
    }

    //Creates kmer mask given an array of base 0-based positions and the kmer size
    public static SVKmerShort getMask(final byte[] positions, final int kSize) {
        long mask = 0;
        for (final byte pos : positions) {
            mask |= 3L << 2*(kSize - pos - 1);
        }
        return new SVKmerShort(~mask);
    }

    public SVKmerShort mask(final SVKmerShort mask) {
        return new SVKmerShort(valLow & mask.valLow);
    }

    /**
     * Not an override.  An SVKmerShort doesn't know what K is, so it has to be supplied.
     */
    public final String toString(final int kSize) {
        final StringBuilder sb = new StringBuilder(kSize);

        // we'll produce the string in reverse order and reverse it at the end
        long val = valLow;
        for (int idx = 0; idx != kSize; ++idx) {
            // grab the two least significant bits to index into the BASE_CHARS array
            sb.append(BaseUtils.BASE_CHARS[(int) val & 3]);
            // roll the whole mess down two bits
            val >>= 2;
        }
        // we built the string in least-significant to most-significant bit order.  reverse it now.
        return sb.reverse().toString();
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SVKmerShort> {
        @Override
        public void write( final Kryo kryo, final Output output, final SVKmerShort svKmer ) {
            svKmer.serialize(kryo, output);
        }

        @Override
        public SVKmerShort read( final Kryo kryo, final Input input, final Class<SVKmerShort> klass ) {
            return new SVKmerShort(kryo, input);
        }
    }

}
