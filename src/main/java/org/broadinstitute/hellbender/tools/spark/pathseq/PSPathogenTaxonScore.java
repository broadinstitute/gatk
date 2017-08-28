package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Pathogen abundance scores assigned to a taxonomic node and reported by the PathSeqScoreSpark tool.
 * See the PathSeqScoreSpark tool for scoring documentation.
 */
public final class PSPathogenTaxonScore {

    final static String DEFAULT_KINGDOM_ID = "0"; //Non-existent id
    public final static String outputHeader = "kingdom\tscore\tscore_normalized\treads\tunambiguous\treference_length";

    public double selfScore = 0; //Total abundance score assigned directly to this taxon
    public double descendentScore = 0; //Sum of descendents' scores
    public double scoreNormalized = 0; //selfScore + descendentScore, normalized to percent of total selfScores
    public int totalReads = 0; //Number of total reads mapped
    public int unambiguousReads = 0; //Number of reads mapped unamibuously to this node
    public long referenceLength = 0; //Length of reference in bp
    public String kingdomTaxonId = DEFAULT_KINGDOM_ID;

    @Override
    public String toString() {
        return toString(kingdomTaxonId);
    }

    public String toString(final PSTree tree) {
        return toString(tree.getNameOf(kingdomTaxonId));
    }

    private String toString(final String kingdomString) {
        return kingdomString + "\t" + (selfScore + descendentScore) + "\t" + scoreNormalized + "\t" + totalReads + "\t" + unambiguousReads + "\t" + referenceLength;
    }

    public PSPathogenTaxonScore add(final PSPathogenTaxonScore other) {
        final PSPathogenTaxonScore result = new PSPathogenTaxonScore();
        result.selfScore = this.selfScore + other.selfScore;
        result.descendentScore = this.descendentScore + other.descendentScore;
        result.totalReads = this.totalReads + other.totalReads;
        result.unambiguousReads = this.unambiguousReads + other.unambiguousReads;
        if (this.referenceLength != other.referenceLength) {
            throw new GATKException("Cannot add PSPathogenTaxonScores with different reference lengths.");
        }
        if (!this.kingdomTaxonId.equals(other.kingdomTaxonId)) {
            throw new GATKException("Cannot add PSPathogenTaxonScores with different kingdoms.");
        }
        result.referenceLength = this.referenceLength;
        return result;
    }
}
