package org.broadinstitute.hellbender.tools.spark.sv;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Created by tsharpe on 12/10/15.
 *
 * Iterator over Kmers from a String.
 * Silently skips over parts of the sequence that has characters other than A, C, G, or T.
 */
public final class Kmerizer implements Iterator<Kmer>
{
    public Kmerizer( String someSeq, int someK )
    {
        seq = someSeq;
        K = someK;
        idx = 0;
        kmer = new Kmer(K);
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
        Kmer result = kmer;
        fill(K-1);
        return result;
    }

    private void fill( int validBaseCount )
    {
        Kmer kkk = kmer;
        kmer = null;

        final int len = seq.length();
        while ( idx < len )
        {
            switch ( seq.charAt(idx) )
            {
                case 'A': kkk = kkk.successor(0,K); break;
                case 'C': kkk = kkk.successor(1,K); break;
                case 'G': kkk = kkk.successor(2,K); break;
                case 'T': kkk = kkk.successor(3,K); break;
                default: validBaseCount = -1;
            }
            idx += 1;

            if ( ++validBaseCount == K )
            {
                kmer = kkk;
                break;
            }
        }
    }

    private final String seq;
    private final int K;
    private int idx;
    private Kmer kmer;
}
