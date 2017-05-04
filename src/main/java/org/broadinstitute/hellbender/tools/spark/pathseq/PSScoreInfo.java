package org.broadinstitute.hellbender.tools.spark.pathseq;

/**
 * Scores assigned to a taxonomic node, used in ClassifyReads
 */
public final class PSScoreInfo {
    public double score = 0; //PathSeq score
    public double score_normalized = 0; //Score normalized to percent of pathogen reads
    public int reads = 0; //Number of total reads mapped
    public int unambiguous = 0; //Number of reads mapped unamibuously to this node
    public long ref_length = 0; //Length of reference in bp

    public final static String outputHeader = "score\tscore_normalized\treads\tunambiguous\treference_length";

    @Override
    public String toString() {
        return score + "\t" + score_normalized + "\t" + reads + "\t" + unambiguous + "\t" + ref_length;
    }
}
