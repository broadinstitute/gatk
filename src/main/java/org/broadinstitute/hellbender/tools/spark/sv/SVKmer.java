package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.BaseUtils;

import java.io.Serializable;

/**
 * An immutable SVKmer.
 * K must be between 1 and 63 (but it's silly to use this class for K < 33).
 * Canonicalization is unimplemented for even K.
 */
public final class SVKmer implements Comparable<SVKmer>, Serializable {
    private static final long serialVersionUID = 1L;

    // these are treated as K-bit unsigned integers
    private final long valHigh; // most significant K bits
    private final long valLow; // least significant K bits

    // Lookup table for reverse-complementing each possible byte value.
    // Each pair of bits represents a base, so you have to reverse bits pairwise and then invert all bits.
    // This is most quickly and easily done with a lookup table.
    private static final long BYTEWISE_REVERSE_COMPLEMENT[];
    static {
        BYTEWISE_REVERSE_COMPLEMENT = new long[256];
        for ( int idx = 0; idx != 256; ++idx ) {
            BYTEWISE_REVERSE_COMPLEMENT[idx] = reverseComplementByteValueAsLong(idx);
        }
    }

    public enum Base {
        A(0L),
        C(1L),
        G(2L),
        T(3L);

        public long value;

        Base(final long value) { this.value = value; }
    }

    /**
     *  Makes an empty SVKmer.  If you call toString on it, it'll look like poly-A.
     */
    public SVKmer( final int kSize ) {
        if ( kSize < 1 || kSize > 63 ) throw new IllegalArgumentException("K must be between 1 and 63.");
        valHigh = valLow = 0;
    }

    private SVKmer( final long valHigh, final long valLow ) { this.valHigh = valHigh; this.valLow = valLow; }

    /**
     * Returns a new SVKmer that's like this one, but with its leading base discarded and a new one added to the end.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.successor(SVKmer.Base.C,5).toString(5) is "CTGAC".
     * @param base must be 0, 1, 2, or 3, corresponding to A, C, G, or T.
     */
    public SVKmer successor( final Base base, final int kSize ) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        final long mask = (1L << kSize) - 1L;
        // move all the bits up two places, OR in the top two bits from valLow at the bottom, and mask to kSize bits
        final long newV1 = ((valHigh << 2) | (valLow >> (kSize-2))) & mask;
        // move all the bits up two places, OR in the pair of successor bits at the bottom, and mask to kSize bits
        final long newV2 = ((valLow << 2) | (base.value & 3L)) & mask;
        return new SVKmer(newV1, newV2);
    }

    /**
     * Returns a new SVKmer that's like this one, but with its trailing base discarded and a new one added to the start.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.predecessor(SVKmer.Base.T,5).toString(5) is "TACTG".
     * @param base must be 0, 1, 2, or 3, corresponding to A, C, G, or T.
     */
    public SVKmer predecessor( final Base base, final int kSize ) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        final long mask = (1L << kSize) - 1L;
        // move all the bits down two places, OR in the successor bits at the top, and mask to kSize bits
        final long newV1 = ((valHigh >> 2) | (base.value << (kSize-2))) & mask;
        // move all the bits down two places, OR in the bottom two bits from valHigh at the top, and mask to kSize bits
        final long newV2 = ((valLow >> 2) | (valHigh << (kSize-2))) & mask;
        return new SVKmer(newV1, newV2);
    }

    /**
     * Returns a new SVKmer that's the reverse-complement of this one.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.rc(5).toString(5) is "TCAGT".
     */
    public SVKmer reverseComplement( final int kSize ) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        final long mask = (1L << kSize) - 1L;
        // number of unused bits at the top
        final int compK = 64 - kSize;
        // if kSize is even
        if ( (kSize&1) == 0 ) {
            // for each val, move the significant bits up to the top, reverse complement, and mask to kSize bits.
            // make a new kmer where the reverse complemented lowVal becomes the highVal and vice versa.
            return new SVKmer(reverseComplement(valLow << compK) & mask, reverseComplement(valHigh << compK) & mask);
        }
        // this is complicated: the middle base's bits straddle the two values (top bit of valLow, bottom of valHigh).
        // it's like the kSize even operation, but we flip the "missing" bit that's in the other value using an XOR,
        // when necessary.
        // i think you're just going to have to trust the unit tests.
        final long newV1 = ((reverseComplement(valLow << (compK+1)) << 1) & mask) ^ ((valHigh & 1L) ^ 1L);
        final long newV2 = (reverseComplement(valHigh << compK) & mask) ^ (valLow & (1L << (kSize-1)));
        return new SVKmer(newV1, newV2);
    }

    /**
     * Returns a SVKmer that is a canonical representation of this one.
     * An odd-K SVKmer is in canonical form if its middle base is A or C.
     * The reverse-complement of a non-canonical SVKmer is a canonical SVKmer, and vice versa.  (Think about it.)
     * Canonical form is not defined for even-K Kmers (too expensive to compute routinely).
     */
    public SVKmer canonical( final int kSize ) {
        if ( (kSize & 1) == 0 ) throw new IllegalArgumentException("K must be odd to canonicalize.");
        // for odd-size kmers, the high bit of the middle base is in least significant position in valHigh.
        // test its value by ANDing with 1.  if it's zero the middle base is A or C and we're good to go.
        if ( (valHigh & 1L) == 0 ) return this;
        // middle base is G or T.  reverse complement.
        return reverseComplement(kSize);
    }

    @Override
    public boolean equals( final Object obj ) {
        return obj instanceof SVKmer && equals((SVKmer)obj);
    }

    public boolean equals( final SVKmer that ) {
        return this.valHigh == that.valHigh && this.valLow == that.valLow;
    }

    @Override
    public int hashCode() { return (int)(valHigh ^ (valHigh >> 32) ^ valLow ^ (valLow >> 32)); }

    /**
     * SVKmer comparison is consistent with equals.
     * It's also the same as the lexicographic ordering you'd get using toString on the Kmers.
     */
    @Override
    public int compareTo( final SVKmer that ) {
        int result = Long.compare(this.valHigh, that.valHigh);
        if ( result == 0 ) result = Long.compare(this.valLow, that.valLow);
        return result;
    }

    /**
     * Not an override.  An SVKmer doesn't know what K is, so it has to be supplied.
     */
    public String toString( final int kSize ) {
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

    private static long reverseComplementByteValueAsLong( final int bIn ) {
        // this turns the 8 bits [b1 b2 b3 b4 b5 b6 b7 b8] into [~b7 ~b8 ~b5 ~b6 ~b3 ~b4 ~b1 ~b2]
        return ~(((bIn & 3) << 6) | (((bIn >> 2) & 3) << 4) | (((bIn >> 4) & 3) << 2) | ((bIn >> 6) & 3)) & 0xffL;
    }

    // Reverse-complement a long by taking the reverse-complement of each of its bytes in reverse order.
    private static long reverseComplement( long val ) {
        // process val one byte at a time
        long result = BYTEWISE_REVERSE_COMPLEMENT[(int)val & 0xFF]; // handle the low-order byte
        int nBytes = 8;
        while ( --nBytes != 0 ) { // pre-decrementing:  we'll go through the loop 7 times
            // rotate down by a byte
            val >>= 8;
            // rotate up by a byte and OR in the reverse complement of the next byte
            result = (result << 8) | BYTEWISE_REVERSE_COMPLEMENT[(int)val & 0xFF];
        }
        return result;
    }
}
