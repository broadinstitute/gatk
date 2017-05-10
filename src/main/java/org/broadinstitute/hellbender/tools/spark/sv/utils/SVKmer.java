package org.broadinstitute.hellbender.tools.spark.sv.utils;

//Superclass for SVKmerLong and SVKmerShort

public abstract class SVKmer {

    public enum Base {
        A(0L),
        C(1L),
        G(2L),
        T(3L);

        public final long value;

        Base(final long value) { this.value = value; }
    }

    // Lookup table for reverse-complementing each possible byte value.
    // Each pair of bits represents a base, so you have to reverse bits pairwise and then invert all bits.
    // This is most quickly and easily done with a lookup table.
    protected static final long[] BYTEWISE_REVERSE_COMPLEMENT;
    static {
        BYTEWISE_REVERSE_COMPLEMENT = new long[256];
        for ( int idx = 0; idx != 256; ++idx ) {
            BYTEWISE_REVERSE_COMPLEMENT[idx] = reverseComplementByteValueAsLong(idx);
        }
    }

    public SVKmer() {}
    public abstract SVKmer successor( final Base base, final int kSize );
    public abstract SVKmer predecessor( final Base base, final int kSize );
    public abstract SVKmer reverseComplement( final int kSize );
    public abstract SVKmer canonical( final int kSize );
    public abstract Base firstBase( final int kSize );
    public abstract Base lastBase();
    public abstract int firstTrimer(final int kSize );
    public abstract int lastTrimer();
    public abstract String toString( final int kSize );

    // Reverse-complement a long by taking the reverse-complement of each of its bytes in reverse order.
    protected static long reverseComplement( long val ) {
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

    protected static long reverseComplementByteValueAsLong( final int bIn ) {
        // this turns the 8 bits [b1 b2 b3 b4 b5 b6 b7 b8] into [~b7 ~b8 ~b5 ~b6 ~b3 ~b4 ~b1 ~b2]
        return ~(((bIn & 3) << 6) | (((bIn >> 2) & 3) << 4) | (((bIn >> 4) & 3) << 2) | ((bIn >> 6) & 3)) & 0xffL;
    }

}
