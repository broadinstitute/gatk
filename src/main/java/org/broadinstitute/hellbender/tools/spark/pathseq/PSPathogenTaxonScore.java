package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Pathogen abundance scores assigned to a taxonomic node and reported by the PathSeqScoreSpark tool.
 * See the PathSeqScoreSpark tool for scoring documentation.
 */
public final class PSPathogenTaxonScore {

    static final int DEFAULT_KINGDOM_ID = 0; //Non-existent id
    static final String outputHeader = String.join("\t","kingdom", "score", "score_normalized", "reads", "unambiguous", "reference_length");

    private double selfScore = 0; //Total abundance score assigned directly to this taxon
    private double descendentScore = 0; //Sum of descendents' scores
    private double scoreNormalized = 0; //selfScore + descendentScore, normalized to percent of total selfScores
    private int totalReads = 0; //Number of total reads mapped
    private int unambiguousReads = 0; //Number of reads mapped unamibuously to this node
    private long referenceLength = 0; //Length of reference in bp
    private int kingdomTaxonId = DEFAULT_KINGDOM_ID;

    @Override
    public String toString() {
        return toString(String.valueOf(kingdomTaxonId));
    }

    public String toString(final PSTree tree) {
        return toString(tree.getNameOf(kingdomTaxonId));
    }

    private String toString(final String kingdomLabel) {
        return String.join("\t", kingdomLabel, String.valueOf(selfScore + descendentScore),
                String.valueOf(scoreNormalized), String.valueOf(totalReads), String.valueOf(unambiguousReads),
                String.valueOf(referenceLength));
    }

    public PSPathogenTaxonScore add(final PSPathogenTaxonScore other) {
        Utils.nonNull(other, "Cannot add taxon score to null");
        if (this.referenceLength != other.referenceLength) {
            throw new GATKException("Cannot add PSPathogenTaxonScores with different reference lengths.");
        }
        if (this.kingdomTaxonId != other.kingdomTaxonId) {
            throw new GATKException("Cannot add PSPathogenTaxonScores with different kingdoms.");
        }
        final PSPathogenTaxonScore result = new PSPathogenTaxonScore();
        result.selfScore = this.selfScore + other.selfScore;
        result.descendentScore = this.descendentScore + other.descendentScore;
        result.totalReads = this.totalReads + other.totalReads;
        result.unambiguousReads = this.unambiguousReads + other.unambiguousReads;
        result.referenceLength = this.referenceLength;
        result.kingdomTaxonId = this.kingdomTaxonId;
        return result;
    }

    public double getSelfScore() {
        return selfScore;
    }

    public double getDescendentScore() {
        return descendentScore;
    }

    public double getScoreNormalized() {
        return scoreNormalized;
    }

    public int getTotalReads() {
        return totalReads;
    }

    public int getUnambiguousReads() {
        return unambiguousReads;
    }

    public long getReferenceLength() {
        return referenceLength;
    }

    public int getKingdomTaxonId() {
        return kingdomTaxonId;
    }

    public void addSelfScore(final double selfScore) {
        Utils.validateArg(selfScore >= 0, "Taxon self score must be non-negative");
        this.selfScore += selfScore;
    }

    public void addDescendentScore(final double descendentScore) {
        Utils.validateArg(descendentScore >= 0, "Taxon descendent score must be non-negative");
        this.descendentScore += descendentScore;
    }

    public void addScoreNormalized(final double scoreNormalized) {
        Utils.validateArg(scoreNormalized >= 0, "Taxon normalized score must be non-negative");
        this.scoreNormalized += scoreNormalized;
    }

    public void addTotalReads(final int totalReads) {
        Utils.validateArg(totalReads >= 0, "Taxon read count must be non-negative");
        this.totalReads += totalReads;
    }

    public void addUnambiguousReads(final int unambiguousReads) {
        Utils.validateArg(unambiguousReads >= 0, "Taxon unambiguous read count must be non-negative");
        this.unambiguousReads += unambiguousReads;
    }

    public void setReferenceLength(final long referenceLength) {
        Utils.validateArg(referenceLength >= 0, "Taxon reference length must be non-negative");
        this.referenceLength = referenceLength;
    }

    public void setKingdomTaxonId(final int kingdomTaxonId) {
        Utils.validateArg(kingdomTaxonId >= 0, "Taxon kingdom must be non-negative");
        this.kingdomTaxonId = kingdomTaxonId;
    }
}
