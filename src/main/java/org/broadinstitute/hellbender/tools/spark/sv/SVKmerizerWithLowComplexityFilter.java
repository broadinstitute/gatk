package org.broadinstitute.hellbender.tools.spark.sv;

import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.tools.spark.sv.SVKmerLong.Base;

/** iterator over kmers with a specified minimum base-wise shannon entropy */
public class SVKmerizerWithLowComplexityFilter extends SVKmerizer {
    private final double minEntropy;
    private final int[] baseCounts;
    private static final Map<Integer, double[]> entropyMap = new ConcurrentHashMap<>();

    public SVKmerizerWithLowComplexityFilter( final byte[] seq, final int kSize, final double minEntropy ) {
        this(new ASCIICharSequence(seq), kSize, minEntropy);
    }

    public SVKmerizerWithLowComplexityFilter( final CharSequence seq, final int kSize, final double minEntropy ) {
        super(kSize, seq);
        this.minEntropy = minEntropy;
        baseCounts = new int[Base.values().length];
        baseCounts[(int)Base.A.value] = kSize;
        this.nextKmer = nextKmer(new SVKmerLong(kSize), 0);
    }

    public static Stream<SVKmer> stream(final CharSequence seq, final int kSize, final double minEntropy ) {
        return StreamSupport.stream(((Iterable<SVKmer>)() ->
                new SVKmerizerWithLowComplexityFilter(seq, kSize, minEntropy)).spliterator(), false);
    }

    public static Stream<SVKmer> stream(final byte[] seq, final int kSize, final double minEntropy ) {
        return StreamSupport.stream(((Iterable<SVKmer>)() ->
                new SVKmerizerWithLowComplexityFilter(seq, kSize, minEntropy)).spliterator(), false);
    }

    @Override
    protected SVKmer nextKmer(SVKmer tmpKmer, int validBaseCount ) {
        final double[] entropies = getEntropies(kSize);
        final int len = seq.length();
        while ( idx < len ) {
            final int firstIdx = (int)tmpKmer.firstBase(kSize).value;
            switch ( seq.charAt(idx) ) {
                case 'a': case 'A':
                    tmpKmer = tmpKmer.successor(Base.A, kSize);
                    baseCounts[firstIdx] -= 1;
                    baseCounts[(int)Base.A.value] += 1;
                    break;
                case 'c': case 'C':
                    tmpKmer = tmpKmer.successor(Base.C, kSize);
                    baseCounts[firstIdx] -= 1;
                    baseCounts[(int)Base.C.value] += 1;
                    break;
                case 'g': case 'G':
                    tmpKmer = tmpKmer.successor(Base.G, kSize);
                    baseCounts[firstIdx] -= 1;
                    baseCounts[(int)Base.G.value] += 1;
                    break;
                case 't': case 'T':
                    tmpKmer = tmpKmer.successor(Base.T, kSize);
                    baseCounts[firstIdx] -= 1;
                    baseCounts[(int)Base.T.value] += 1;
                    break;
                default:
                    validBaseCount = -1;
                    break;
            }
            idx += 1;

            if ( ++validBaseCount >= kSize ) {
                final double entropy = Arrays.stream(baseCounts).mapToDouble(count -> entropies[count]).sum();
                if ( entropy >= minEntropy ) return tmpKmer;
            }
        }
        return null;
    }

    /**
     * This returns an array where the i'th entry is the contribution to the total entropy of a kmer of some
     * base that occurs i times in a kmer of the given size.
     */
    private static double[] getEntropies( final int kSize ) {
        double[] result = entropyMap.get(kSize);
        if ( result == null ) {
            result = new double[kSize+1];
            final double minusLog2 = -Math.log(2);
            // note: result[0] = 0 by definition, and it's already set to that value.
            // the same is true for result[kSize].
            for ( int idx = 1; idx != kSize; ++idx ) {
                final double probability = 1.0 * idx / kSize;
                result[idx] = probability * Math.log(probability) / minusLog2;
            }
            entropyMap.putIfAbsent(kSize, result);
        }
        return result;
    }
}
