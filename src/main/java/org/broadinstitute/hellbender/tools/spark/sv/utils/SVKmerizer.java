package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Iterator over successive Kmers from a sequence of characters.
 * Silently skips over parts of the sequence that has characters other than A, C, G, or T.
 */
public class SVKmerizer<KmerType extends SVKmer> implements Iterator<KmerType> {
    protected final CharSequence seq;
    protected final int kSize, kAdvance;
    protected int idx = 0;
    protected KmerType nextKmer;

    //3-arg constructor creates a Kmerizer that returns consecutive kmers in the sequence
    public SVKmerizer( final byte[] seq, final int kSize, final KmerType kmer ) {
        this(seq, kSize, 1, kmer );
    }

    //Overloaded constructor with kSpace, the number of bases to advance between kmers
    public SVKmerizer( final byte[] seq, final int kSize, final int kSpace, final KmerType kmer ) {
        this(new ASCIICharSequence(seq), kSize, kSpace, kmer );
    }

    public SVKmerizer( final CharSequence seq, final int kSize, final KmerType kmer ) {
        this(seq, kSize, 1, kmer);
    }

    public SVKmerizer( final CharSequence seq, final int kSize, final int kSpace, final KmerType kmer ) {
        this.seq = seq;
        this.kSize = kSize;
        this.kAdvance = kSize - kSpace;
        this.nextKmer = nextKmer(kmer, 0);
    }

    protected SVKmerizer( final int kSize, final CharSequence seq ) {
        this(kSize, 1, seq);
    }

    protected SVKmerizer( final int kSize, final int kSpace, final CharSequence seq ) {
        this.seq = seq;
        this.kSize = kSize;
        this.kAdvance = kSize - kSpace;
    }

    @Override
    public boolean hasNext() {
        return nextKmer != null;
    }

    @Override
    public KmerType next() {
        if ( nextKmer == null ) throw new NoSuchElementException("Kmerization sequence exhausted.");
        final KmerType result = nextKmer;
        nextKmer = nextKmer(nextKmer, kAdvance);
        return result;
    }

    /** returns distance from beginning of sequence to beginning of current kmer */
    public int getOffset() { return idx - kSize; }

    public static <KmerType extends SVKmer> KmerType toKmer( final CharSequence seq, final KmerType kmer ) {
        final SVKmerizer<KmerType> sk = new SVKmerizer<>(seq, seq.length(), 1, kmer);
        Utils.validateArg(sk.hasNext(), () -> "Can't make a SVKmer from '"+seq+"'");
        return sk.next();
    }

    public static <KmerType extends SVKmer> KmerType toKmer( final byte[] seq, final KmerType kmer ) {
        return toKmer(new ASCIICharSequence(seq),kmer);
    }

    public static <KmerType extends SVKmer> Stream<KmerType> stream( final CharSequence seq,
                                                                     final int kSize,
                                                                     final int kSpace,
                                                                     final KmerType kmer ) {
        return StreamSupport.stream(((Iterable<KmerType>)() ->
                new SVKmerizer<>(seq, kSize, kSpace, kmer)).spliterator(), false);
    }

    public static <KmerType extends SVKmer> Stream<KmerType> stream( final CharSequence seq,
                                                                     final int kSize,
                                                                     final KmerType kmer ) {
        return stream(seq, kSize, 1, kmer);
    }

    public static <KmerType extends SVKmer> Stream<KmerType> stream( final byte[] seq,
                                                                     final int kSize,
                                                                     final int kSpace,
                                                                     final KmerType kmer ) {
        return stream(new ASCIICharSequence(seq), kSize, kSpace, kmer);
    }

    public static <KmerType extends SVKmer> Stream<KmerType> stream( final byte[] seq,
                                                                     final int kSize,
                                                                     final KmerType kmer ) {
        return stream(seq, kSize, 1, kmer);
    }

    @SuppressWarnings("unchecked")
    public static <KmerType extends SVKmer> Stream<KmerType> canonicalStream( final byte[] seq,
                                                                              final int kSize,
                                                                              final KmerType kmer ) {
        return stream(seq, kSize, 1, kmer).map(kkk -> (KmerType)kkk.canonical(kSize));
    }

    @SuppressWarnings("unchecked")
    public static <KmerType extends SVKmer> Stream<KmerType> canonicalStream( final CharSequence seq,
                                                                              final int kSize,
                                                                              final KmerType kmer ) {
        return stream(seq, kSize, 1, kmer).map(kkk -> (KmerType)kkk.canonical(kSize));
    }

    @SuppressWarnings("unchecked")
    protected KmerType nextKmer( final KmerType tmpKmerArg, final int validBaseCountArg ) {
        KmerType tmpKmer = tmpKmerArg;
        int validBaseCount = validBaseCountArg;
        final int len = seq.length();
        while ( idx < len ) {
            switch ( seq.charAt(idx) ) {
                case 'a': case 'A': tmpKmer = (KmerType)tmpKmer.successor(SVKmer.Base.A, kSize); break;
                case 'c': case 'C': tmpKmer = (KmerType)tmpKmer.successor(SVKmer.Base.C, kSize); break;
                case 'g': case 'G': tmpKmer = (KmerType)tmpKmer.successor(SVKmer.Base.G, kSize); break;
                case 't': case 'T': tmpKmer = (KmerType)tmpKmer.successor(SVKmer.Base.T, kSize); break;
                default: validBaseCount = -1;
            }
            idx += 1;

            if ( ++validBaseCount == kSize ) return tmpKmer;
        }
        return null;
    }

    // a shim to turn a byte array into a character sequence by treating the bytes as ASCII characters
    public static final class ASCIICharSequence implements CharSequence {
        public ASCIICharSequence( final byte[] bytes ) { this.bytes = bytes; }

        @Override public int length() { return bytes.length; }

        @Override public char charAt( final int index ) { return (char)(bytes[index] & 0xff); }

        @Override public CharSequence subSequence( final int start, final int end ) {
            return new ASCIICharSequence(Arrays.copyOfRange(bytes, start, end));
        }

        @Override public String toString() { return new StringBuilder(this).toString(); }

        final byte[] bytes;
    }
}
