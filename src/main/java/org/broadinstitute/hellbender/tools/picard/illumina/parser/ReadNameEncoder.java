package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/**
 * @author mccowan
 */
public interface ReadNameEncoder {
    /**
     * Generates a read name string for the provided cluster. 
     *
     * @param cluster The cluster whose reads are having its name generated
     * @param pairNumber 1 if this is the first of the pair, 2 if it is the second, or null if this not a paired read.
     * @return The read name
     */
    String generateReadName(ClusterData cluster, Integer pairNumber);
}
