package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * An immutable SVKmerLong.
 * K must be between 1 and 63 (but it's silly to use this class for K < 33).
 * Canonicalization is unimplemented for even K.
 */
@DefaultSerializer(SVKmerLong.Serializer.class)
public class SVKmerLong extends SVKmer implements Comparable<SVKmerLong>  {
    // these are treated as K-bit unsigned integers
    private final long valHigh; // most significant K bits
    private final long valLow; // least significant K bits

    public SVKmerLong() {
        valHigh = valLow = 0;
    }

    /**
     *  Makes an empty SVKmerLong.  If you call toString on it, it'll look like poly-A.
     */
    public SVKmerLong( final int kSize ) {
        Utils.validateArg(kSize >= 1 && kSize < 64, "Kmer length must be between 1 and 63.");
        valHigh = valLow = 0;
    }

    public SVKmerLong( final SVKmerLong that ) { this.valHigh = that.valHigh; this.valLow = that.valLow; }

    public SVKmerLong( final SVKmer that ) {
        SVKmerLong thatLong = (SVKmerLong)that;
        this.valHigh = thatLong.valHigh;
        this.valLow = thatLong.valLow;
    }

    private SVKmerLong( final long valHigh, final long valLow ) { this.valHigh = valHigh; this.valLow = valLow; }

    protected SVKmerLong( final Kryo kryo, final Input input ) {
        valHigh = input.readLong();
        valLow = input.readLong();
    }

    protected void serialize( final Kryo kryo, final Output output ) {
        output.writeLong(valHigh);
        output.writeLong(valLow);
    }

    /**
     * Returns a new SVKmerLong that's like this one, but with its leading base discarded and a new one added to the end.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.successor(SVKmerLong.Base.C,5).toString(5) is "CTGAC".
     * @param base must be 0, 1, 2, or 3, corresponding to A, C, G, or T.
     */
    public final SVKmerLong successor( final Base base, final int kSize ) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        final long mask = (1L << kSize) - 1L;
        // move all the bits up two places, OR in the top two bits from valLow at the bottom, and mask to kSize bits
        final long newV1 = ((valHigh << 2) | (valLow >> (kSize-2))) & mask;
        // move all the bits up two places, OR in the pair of successor bits at the bottom, and mask to kSize bits
        final long newV2 = ((valLow << 2) | (base.value & 3L)) & mask;
        return new SVKmerLong(newV1, newV2);
    }

    /**
     * Returns a new SVKmerLong that's like this one, but with its trailing base discarded and a new one added to the start.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.predecessor(SVKmerLong.Base.T,5).toString(5) is "TACTG".
     * @param base must be 0, 1, 2, or 3, corresponding to A, C, G, or T.
     */
    public final SVKmerLong predecessor( final Base base, final int kSize ) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        final long mask = (1L << kSize) - 1L;
        // move all the bits down two places, OR in the successor bits at the top, and mask to kSize bits
        final long newV1 = ((valHigh >> 2) | (base.value << (kSize-2))) & mask;
        // move all the bits down two places, OR in the bottom two bits from valHigh at the top, and mask to kSize bits
        final long newV2 = ((valLow >> 2) | (valHigh << (kSize-2))) & mask;
        return new SVKmerLong(newV1, newV2);
    }

    /**
     * Returns a new SVKmerLong that's the reverse-complement of this one.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.rc(5).toString(5) is "TCAGT".
     */
    public final SVKmerLong reverseComplement( final int kSize ) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        final long mask = (1L << kSize) - 1L;
        // number of unused bits at the top
        final int compK = 64 - kSize;
        // if kSize is even
        if ( (kSize&1) == 0 ) {
            // for each val, move the significant bits up to the top, reverse complement, and mask to kSize bits.
            // make a new kmer where the reverse complemented lowVal becomes the highVal and vice versa.
            return new SVKmerLong(reverseComplement(valLow << compK) & mask, reverseComplement(valHigh << compK) & mask);
        }
        // this is complicated: the middle base's bits straddle the two values (top bit of valLow, bottom of valHigh).
        // it's like the kSize even operation, but we flip the "missing" bit that's in the other value using an XOR,
        // when necessary.
        // i think you're just going to have to trust the unit tests.
        final long newV1 = ((reverseComplement(valLow << (compK+1)) << 1) & mask) ^ ((valHigh & 1L) ^ 1L);
        final long newV2 = (reverseComplement(valHigh << compK) & mask) ^ (valLow & (1L << (kSize-1)));
        return new SVKmerLong(newV1, newV2);
    }

    /**
     * Returns a SVKmerLong that is a canonical representation of this one.
     * An odd-K SVKmerLong is in canonical form if its middle base is A or C.
     * The reverse-complement of a non-canonical SVKmerLong is a canonical SVKmerLong, and vice versa.  (Think about it.)
     * Canonical form is not defined for even-K Kmers (too expensive to compute routinely).
     */
    public SVKmerLong canonical( final int kSize ) {
        Utils.validateArg( (kSize & 1) != 0, "Kmer length must be odd to canonicalize.");
        // for odd-size kmers, the high bit of the middle base is in least significant position in valHigh.
        // test its value by ANDing with 1.  if it's zero the middle base is A or C and we're good to go.
        if ( (valHigh & 1L) == 0 ) return this;
        // middle base is G or T.  reverse complement.
        return reverseComplement(kSize);
    }

    public final Base firstBase( final int kSize ) { return Base.values()[(int)(valHigh >> (kSize-2))]; }
    public final Base lastBase() { return Base.values()[(int)(valLow & 3)]; }
    public final int firstTrimer(final int kSize ) { return (int)(valHigh >>> (kSize-6)); }
    public final int lastTrimer() { return (int)valLow & 0x3F; }
    @Override
    public boolean equals( final Object obj ) {
        return obj instanceof SVKmerLong && equals((SVKmerLong)obj);
    }

    public final boolean equals( final SVKmerLong that ) {
        return this.valHigh == that.valHigh && this.valLow == that.valLow;
    }

    @Override
    public final int hashCode() {
        return (int) SVUtils.fnvLong64(SVUtils.fnvLong64(valHigh), valLow);
    }

    /**
     * SVKmerLong comparison is consistent with equals.
     * It's also the same as the lexicographic ordering you'd get using toString on the Kmers.
     */
    @Override
    public final int compareTo( final SVKmerLong that ) {
        int result = Long.compare(this.valHigh, that.valHigh);
        if ( result == 0 ) result = Long.compare(this.valLow, that.valLow);
        return result;
    }

    /**
     * Not an override.  An SVKmerLong doesn't know what K is, so it has to be supplied.
     */
    public final String toString( final int kSize ) {
        final StringBuilder sb = new StringBuilder(kSize);

        // we'll produce the string in reverse order and reverse it at the end
        long val = valLow;
        for ( int nChars = kSize/2; nChars > 0; --nChars ) {
            // grab the two least significant bits to index into the BASE_CHARS array
            sb.append(BaseUtils.BASE_CHARS[(int)val & 3]);
            // roll the whole mess down two bits
            val >>= 2;
        }
        // if kSize is even
        if ( (kSize & 1) == 0 ) {
            val = valHigh; // we've used all the bits from valLow -- just move on to valHigh
        } else { // kSize is odd
            val |= valHigh << 1; // there's one leftover bit from valLow that need to be accounted for
        }
        // this for loop will have one more iteration than the previous one if kSize is odd
        for ( int nChars = (kSize+1)/2; nChars > 0; --nChars ) {
            // grab two lowest bits to index into array
            sb.append(BaseUtils.BASE_CHARS[(int)val & 3]);
            // move 'em down
            val >>= 2;
        }
        // we built the string in least-significant to most-significant bit order.  reverse it now.
        return sb.reverse().toString();
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SVKmerLong> {
        @Override
        public void write( final Kryo kryo, final Output output, final SVKmerLong svKmer ) {
            svKmer.serialize(kryo, output);
        }

        @Override
        public SVKmerLong read( final Kryo kryo, final Input input, final Class<SVKmerLong> klass ) {
            return new SVKmerLong(kryo, input);
        }
    }

}
