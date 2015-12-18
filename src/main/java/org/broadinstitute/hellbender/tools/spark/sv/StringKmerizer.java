package org.broadinstitute.hellbender.tools.spark.sv;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Iterator over successive Kmers from a String.
 * Silently skips over parts of the sequence that has characters other than A, C, G, or T.
 *
 * Created by tsharpe on 12/10/15.
 */
public final class StringKmerizer implements Iterator<Kmer>
{
    private final String seq;
    private final int kSize;
    private int idx;
    private Kmer kmer;

    public StringKmerizer( final String someSeq, final int someKSize )
    {
        seq = someSeq;
        kSize = someKSize;
        idx = 0;
        kmer = new Kmer(kSize);
        fill(0);
    }

    @Override
    public boolean hasNext()
    {
        return kmer != null;
    }

    @Override
    public Kmer next()
    {
        if ( kmer == null ) throw new NoSuchElementException("Kmerization sequence exhausted.");
        final Kmer result = kmer;
        fill(kSize -1);
        return result;
    }

    public static Kmer toKmer( final String seq )
    {
        final StringKmerizer sk = new StringKmerizer(seq,seq.length());
        if ( !sk.hasNext() )
            throw new IllegalArgumentException("Can't make a Kmer from '"+seq+"'");
        return sk.next();
    }

    private void fill( int validBaseCount )
    {
        Kmer tmpKmer = kmer;
        kmer = null;

        final int len = seq.length();
        while ( idx < len )
        {
            switch ( seq.charAt(idx) )
            {
                case 'A': tmpKmer = tmpKmer.successor(0, kSize); break;
                case 'C': tmpKmer = tmpKmer.successor(1, kSize); break;
                case 'G': tmpKmer = tmpKmer.successor(2, kSize); break;
                case 'T': tmpKmer = tmpKmer.successor(3, kSize); break;
                default: validBaseCount = -1;
            }
            idx += 1;

            if ( ++validBaseCount == kSize )
            {
                kmer = tmpKmer;
                break;
            }
        }
    }
}
