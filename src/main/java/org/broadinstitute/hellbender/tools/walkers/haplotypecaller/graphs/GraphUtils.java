package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Utility functions used in the graphs package
 */
public final class GraphUtils {
    private GraphUtils() {}

    /**
     * Compute the maximum shared prefix length of list of bytes.
     *
     * @param listOfBytes a list of bytes with at least one element
     * @return the number of shared bytes common at the start of all bytes
     */
    public static int commonMaximumPrefixLength(final List<byte[]> listOfBytes) {
        final int minLength = GraphUtils.minKmerLength(listOfBytes);
        for ( int i = 0; i < minLength; i++ ) {
            final byte b = listOfBytes.get(0)[i];
            for ( int j = 1; j < listOfBytes.size(); j++ ) {
                if ( b != listOfBytes.get(j)[i] ) {
                    return i;
                }
            }
        }

        return minLength;
    }

    /**
     * Compute the maximum shared suffix length of list of bytes.
     *
     * @param listOfBytes a list of bytes with at least one element
     * @param minLength the min. length among all byte[] in listOfBytes
     * @return the number of shared bytes common at the end of all bytes
     */
    public static int commonMaximumSuffixLength(final List<byte[]> listOfBytes, final int minLength) {
        for ( int suffixLen = 0; suffixLen < minLength; suffixLen++ ) {
            final byte b = listOfBytes.get(0)[listOfBytes.get(0).length - suffixLen - 1];
            for ( int j = 1; j < listOfBytes.size(); j++ ) {
                if ( b != listOfBytes.get(j)[listOfBytes.get(j).length - suffixLen - 1] ) {
                    return suffixLen;
                }
            }
        }
        return minLength;
    }

    /**
     * Get the list of kmers as byte[] from the vertices in the graph
     *
     * @param vertices a collection of vertices
     * @return a list of their kmers in order of the iterator on vertices
     */
    public static List<byte[]> getKmers(final Collection<SeqVertex> vertices) {
        return Utils.nonNull(vertices).stream().map(v -> v.getSequence()).collect(Collectors.toList());
    }

    /**
     * Get the minimum length of a collection of byte[]
     *
     * @param kmers a list of kmers whose .length min we want
     * @return the min of the kmers, if kmers is empty the result is 0
     */
    public static int minKmerLength(final Collection<byte[]> kmers) {
        return Utils.nonNull(kmers).stream().mapToInt(kmer -> kmer.length).min().orElse(0);
    }

}
