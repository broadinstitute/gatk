package org.broadinstitute.hellbender.tools.spark.sv;

import java.io.Serializable;

/**
 * Created by tsharpe on 12/10/15.
 *
 * An immutable Kmer.
 * K must be between 1 and 63 (but it's silly to use for K < 33).
 * Canonicalization is unimplemented for even K.
 */
public final class Kmer implements Comparable<Kmer>, Serializable
{
    private static final long serialVersionUID = 1L;

    public Kmer( final int K )
    {
        if ( K < 1 || K > 63 ) throw new IllegalArgumentException("K must be between 1 and 63.");
        v1 = v2 = 0;
    }

    public Kmer successor( final long base, final int K )
    {
        final long mask = (1L << K) - 1L;
        final long newV1 = ((v1 << 2) | (v2 >> (K-2))) & mask;
        final long newV2 = ((v2 << 2) | (base & 3L)) & mask;
        return new Kmer(newV1,newV2);
    }

    public Kmer predecessor( final long base, final int K )
    {
        final long mask = (1L << K) - 1L;
        final long newV1 = ((v1 >> 2) | (base << (K-2))) & mask;
        final long newV2 = ((v2 >> 2) | (v1 << (K-2))) & mask;
        return new Kmer(newV1,newV2);
    }

    public Kmer rc( final int K )
    {
        final long mask = (1L << K) - 1L;
        int compK = 64 - K;
        if ( (K&1) == 0 )
            return new Kmer(rcLong(v2 << compK)&mask,rcLong(v1 << compK)&mask);
        long newV1 = ((rcLong(v2 << (compK+1)) << 1) & mask) ^ ((v1 & 1L) ^ 1L);
        long newV2 = (rcLong(v1 << compK) & mask) ^ (v2 & (1L << (K-1)));
        return new Kmer(newV1, newV2);
    }

    // an odd-K kmer is in canonical form if its middle base is A or C
    public Kmer canonical( final int K )
    {
        if ( (K & 1) == 0 ) throw new IllegalArgumentException("K must be odd to canonicalize.");
        if ( (v1 & 1L) == 0 ) return this;
        return rc(K);
    }

    @Override
    public boolean equals( final Object obj )
    {
        if ( !(obj instanceof Kmer) ) return false;
        final Kmer that = (Kmer)obj;
        return this.v1 == that.v1 && this.v2 == that.v2;
    }

    public boolean equals( final Kmer that )
    {
        return this.v1 == that.v1 && this.v2 == that.v2;
    }

    @Override
    public int hashCode() { return (int)(v1 ^ (v1 >> 32) ^ v2 ^ (v2 >> 32)); }

    @Override
    public int compareTo( final Kmer that )
    {
        int result = Long.compare(this.v1,that.v1);
        if ( result == 0 ) result = Long.compare(this.v2,that.v2);
        return result;
    }

    public String toString( final int K )
    {
        final StringBuilder sb = new StringBuilder(K);
        long val = v2;
        int nnn = K >> 1;
        while ( nnn-- != 0 )
        {
            sb.append(baseChars[(int)val & 3]);
            val >>= 2;
        }
        if ( (K & 1) == 0 )
            val = v1;
        else
            val |= v1 << 1;
        nnn = (K + 1) >> 1;
        while ( nnn-- != 0 )
        {
            sb.append(baseChars[(int)val & 3]);
            val >>= 2;
        }
        return sb.reverse().toString();
    }


    private Kmer( final long aV1, final long aV2 ) { v1 = aV1; v2 = aV2; }

    // these are treated as K-bit unsigned integers
    private final long v1; // most significant K bits
    private final long v2; // least significant K bits

    private static long rcLong( long val )
    {
        long result = rcByte[(int)val & 0xFF];
        int nBytes = 8;
        while ( --nBytes != 0 )
        {
            val >>= 8;
            result = (result << 8) | rcByte[(int)val & 0xFF];
        }
        return result;
    }

    private static final char baseChars[] = { 'A', 'C', 'G', 'T' };
    private static final long rcByte[] =
    {
        0xff, 0xbf, 0x7f, 0x3f, 0xef, 0xaf, 0x6f, 0x2f,
        0xdf, 0x9f, 0x5f, 0x1f, 0xcf, 0x8f, 0x4f, 0x0f,
        0xfb, 0xbb, 0x7b, 0x3b, 0xeb, 0xab, 0x6b, 0x2b,
        0xdb, 0x9b, 0x5b, 0x1b, 0xcb, 0x8b, 0x4b, 0x0b,
        0xf7, 0xb7, 0x77, 0x37, 0xe7, 0xa7, 0x67, 0x27,
        0xd7, 0x97, 0x57, 0x17, 0xc7, 0x87, 0x47, 0x07,
        0xf3, 0xb3, 0x73, 0x33, 0xe3, 0xa3, 0x63, 0x23,
        0xd3, 0x93, 0x53, 0x13, 0xc3, 0x83, 0x43, 0x03,
        0xfe, 0xbe, 0x7e, 0x3e, 0xee, 0xae, 0x6e, 0x2e,
        0xde, 0x9e, 0x5e, 0x1e, 0xce, 0x8e, 0x4e, 0x0e,
        0xfa, 0xba, 0x7a, 0x3a, 0xea, 0xaa, 0x6a, 0x2a,
        0xda, 0x9a, 0x5a, 0x1a, 0xca, 0x8a, 0x4a, 0x0a,
        0xf6, 0xb6, 0x76, 0x36, 0xe6, 0xa6, 0x66, 0x26,
        0xd6, 0x96, 0x56, 0x16, 0xc6, 0x86, 0x46, 0x06,
        0xf2, 0xb2, 0x72, 0x32, 0xe2, 0xa2, 0x62, 0x22,
        0xd2, 0x92, 0x52, 0x12, 0xc2, 0x82, 0x42, 0x02,
        0xfd, 0xbd, 0x7d, 0x3d, 0xed, 0xad, 0x6d, 0x2d,
        0xdd, 0x9d, 0x5d, 0x1d, 0xcd, 0x8d, 0x4d, 0x0d,
        0xf9, 0xb9, 0x79, 0x39, 0xe9, 0xa9, 0x69, 0x29,
        0xd9, 0x99, 0x59, 0x19, 0xc9, 0x89, 0x49, 0x09,
        0xf5, 0xb5, 0x75, 0x35, 0xe5, 0xa5, 0x65, 0x25,
        0xd5, 0x95, 0x55, 0x15, 0xc5, 0x85, 0x45, 0x05,
        0xf1, 0xb1, 0x71, 0x31, 0xe1, 0xa1, 0x61, 0x21,
        0xd1, 0x91, 0x51, 0x11, 0xc1, 0x81, 0x41, 0x01,
        0xfc, 0xbc, 0x7c, 0x3c, 0xec, 0xac, 0x6c, 0x2c,
        0xdc, 0x9c, 0x5c, 0x1c, 0xcc, 0x8c, 0x4c, 0x0c,
        0xf8, 0xb8, 0x78, 0x38, 0xe8, 0xa8, 0x68, 0x28,
        0xd8, 0x98, 0x58, 0x18, 0xc8, 0x88, 0x48, 0x08,
        0xf4, 0xb4, 0x74, 0x34, 0xe4, 0xa4, 0x64, 0x24,
        0xd4, 0x94, 0x54, 0x14, 0xc4, 0x84, 0x44, 0x04,
        0xf0, 0xb0, 0x70, 0x30, 0xe0, 0xa0, 0x60, 0x20,
        0xd0, 0x90, 0x50, 0x10, 0xc0, 0x80, 0x40, 0x00
    };
}
