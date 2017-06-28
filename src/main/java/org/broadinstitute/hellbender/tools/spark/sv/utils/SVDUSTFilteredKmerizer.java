package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.stream.Stream;

/** An iterator over kmers with a specified maximum DUST-style, low-complexity score.
 *
 *  DUST scores a sequence by counting the number of occurrences of each of the 64 possible trimers in that sequence,
 *  and summing these counts according to:  sumOverTrimerCounts( count[t]*(count[t]-1)/2 ) / (sequence.length - 3)
 *  (Where count[t] is the number of occurrences of the trimer having index t in the sequence under consideration,
 *   and sumOverTrimerCounts just a for-loop accumulator that lets the index t run from 0 to 63).
 *  See "DUST for masking low-complexity DNA sequences", Morgulies et al., 2006 [DOI: 10.1089/cmb.2006.13.1028].
 *
 *  Since we're dealing with kmers of a fixed sequence length, we ignore the division by sequence length.
 *  So our DUST-like score is (kSize-3) too big, but proportional to a real DUST score for the kmer,
 *  which is good enough and much faster to calculate.
 */
public class SVDUSTFilteredKmerizer extends SVKmerizer {
    private final int maxDUSTScore;
    private final int[] trimerCounts;
    private int curDUSTScore;

    public SVDUSTFilteredKmerizer(final byte[] seq, final int kSize, final int maxDUSTScore, final SVKmer kmer ) {
        this(new ASCIICharSequence(seq), kSize, maxDUSTScore, kmer);
    }

    public SVDUSTFilteredKmerizer(final CharSequence seq, final int kSize, final int maxDUSTScore, final SVKmer kmer ) {
        super(kSize, seq);
        // kmers of length < 3 don't have any trimers.
        // kmers of length 3 always have a score of 0, and are uninformative.
        // kSize probably ought to be considerably larger than 4, but we know that numbers less than 4 are plain silly.
        if ( kSize < 4 ) {
            throw new GATKException("kmer size must be at least 4 for this filter to work properly.");
        }
        this.maxDUSTScore = maxDUSTScore;
        trimerCounts = new int[64]; // there are 64 different trimers
        final int polyACount = kSize - 2; // the current kmer is poly-A, so trimer 0 (AAA) has kSize-2 counts
        trimerCounts[0] = polyACount;
        curDUSTScore = polyACount*(polyACount - 1)/2; // see formula above
        nextKmer = nextKmer(kmer, 0);
    }

    public static Stream<SVKmer> stream( final CharSequence seq, final int kSize, final int maxDUSTScore, final SVKmer kmer ) {
        return Utils.stream(new SVDUSTFilteredKmerizer(seq, kSize, maxDUSTScore, kmer));
    }

    public static Stream<SVKmer> stream( final byte[] seq, final int kSize, final int maxDUSTScore, final SVKmer kmer ) {
        return Utils.stream(new SVDUSTFilteredKmerizer(seq, kSize, maxDUSTScore, kmer));
    }

    @Override
    protected SVKmer nextKmer( final SVKmer initialKmer, int validBaseCount ) {
        final int len = seq.length();
        SVKmer result = initialKmer;
        while ( idx < len ) {
            // adjust score for the contribution of the lead trimer, which is disappearing
            // note that count*(count-1)/2 is equal to sum(0, 1, 2, ..., (count-1))
            curDUSTScore -= --trimerCounts[result.firstTrimer(kSize)];
            switch ( seq.charAt(idx) ) {
                case 'a': case 'A':
                    result = result.successor(SVKmer.Base.A, kSize);
                    break;
                case 'c': case 'C':
                    result = result.successor(SVKmer.Base.C, kSize);
                    break;
                case 'g': case 'G':
                    result = result.successor(SVKmer.Base.G, kSize);
                    break;
                case 't': case 'T':
                    result = result.successor(SVKmer.Base.T, kSize);
                    break;
                default:
                    result = result.successor(SVKmer.Base.A, kSize);
                    validBaseCount = -1;
                    break;
            }
            // adjust score for the new trailing trimer that just got rolled in
            curDUSTScore += trimerCounts[result.lastTrimer()]++;
            idx += 1;

            if ( ++validBaseCount >= kSize && curDUSTScore <= maxDUSTScore ) return result;
        }
        return null;
    }
}
