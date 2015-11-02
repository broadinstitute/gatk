package org.broadinstitute.hellbender.tools.spark.sv;

import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Iterator over successive Kmers from a sequence of characters.
 * Silently skips over parts of the sequence that has characters other than A, C, G, or T.
 */
public final class SVKmerizer implements Iterator<SVKmer> {
    private final CharSequence seq;
    private final int kSize;
    private int idx = 0;
    private SVKmer nextKmer;

    public SVKmerizer( final CharSequence seq, final int kSize ) {
        this.seq = seq;
        this.kSize = kSize;
        this.nextKmer = nextKmer(new SVKmer(kSize), 0);
    }

    public SVKmerizer( final byte[] seq, final int kSize ) {
        this(new ASCIICharSequence(seq), kSize);
    }

    @Override
    public boolean hasNext() {
        return nextKmer != null;
    }

    @Override
    public SVKmer next() {
        if ( nextKmer == null ) throw new NoSuchElementException("Kmerization sequence exhausted.");
        final SVKmer result = nextKmer;
        nextKmer = nextKmer(nextKmer, kSize-1);
        return result;
    }

    public static SVKmer toKmer( final CharSequence seq ) {
        final SVKmerizer sk = new SVKmerizer(seq, seq.length());
        if ( !sk.hasNext() ) throw new IllegalArgumentException("Can't make a SVKmer from '"+seq+"'");
        return sk.next();
    }

    public static SVKmer toKmer( final byte[] seq ) {
        return toKmer(new ASCIICharSequence(seq));
    }

    public static Stream<SVKmer> stream( final CharSequence seq, final int kSize ) {
        return StreamSupport.stream(((Iterable<SVKmer>)() -> new SVKmerizer(seq, kSize)).spliterator(), false);
    }

    public static Stream<SVKmer> stream( final byte[] seq, final int kSize ) {
        return stream(new ASCIICharSequence(seq),kSize);
    }

    private SVKmer nextKmer( SVKmer tmpKmer, int validBaseCount ) {
        final int len = seq.length();
        while ( idx < len ) {
            switch ( seq.charAt(idx) ) {
                case 'a': case 'A': tmpKmer = tmpKmer.successor(SVKmer.Base.A, kSize); break;
                case 'c': case 'C': tmpKmer = tmpKmer.successor(SVKmer.Base.C, kSize); break;
                case 'g': case 'G': tmpKmer = tmpKmer.successor(SVKmer.Base.G, kSize); break;
                case 't': case 'T': tmpKmer = tmpKmer.successor(SVKmer.Base.T, kSize); break;
                default: validBaseCount = -1;
            }
            idx += 1;

            if ( ++validBaseCount == kSize ) return tmpKmer;
        }
        return null;
    }

    // a shim to turn a byte array into a character sequence by treating the bytes as ASCII characters
    public static final class ASCIICharSequence implements CharSequence {
        ASCIICharSequence( final byte[] bytes ) { this.bytes = bytes; }

        @Override public int length() { return bytes.length; }

        @Override public char charAt( final int index ) { return (char)(bytes[index] & 0xff); }

        @Override public CharSequence subSequence( final int start, final int end ) {
            return new ASCIICharSequence(Arrays.copyOfRange(bytes, start, end));
        }

        @Override public String toString() { return new StringBuilder(this).toString(); }

        final byte[] bytes;
    }
}
