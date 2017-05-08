package org.broadinstitute.hellbender.tools.spark.pathseq;

/**
 * Pathogen abundance scores assigned to a taxonomic node and reported by the ClassifyReads tool.
 * See the ClassifyReads tool for scoring documentation.
 */
public final class PSPathogenTaxonScore {
    public double score = 0; //PathSeq abundance score calculated in the ClassifyReads tool
    public double scoreNormalized = 0; //Score normalized to percent of total
    public int reads = 0; //Number of total reads mapped
    public int unambiguous = 0; //Number of reads mapped unamibuously to this node
    public long refLength = 0; //Length of reference in bp

    public final static String outputHeader = "score\tscore_normalized\treads\tunambiguous\treference_length";

    @Override
    public String toString() {
        return score + "\t" + scoreNormalized + "\t" + reads + "\t" + unambiguous + "\t" + refLength;
    }
}
