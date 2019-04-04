package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

public class Datum {
    private final double tumorLogOdds;
    private final double artifactProb;
    private final double nonSequencingErrorProb;
    private final int altCount;
    private final int totalCount;
    private final int indelLength; // positive for insertion, negative for deletion, zero for SNV

    public Datum(final double tumorLogOdds, final double artifactProb, final double nonSomaticProb, final int altCount, final int totalCount, final int indelLength) {
        this.tumorLogOdds = tumorLogOdds;
        this.artifactProb = artifactProb;
        this.altCount = altCount;
        this.totalCount = totalCount;
        this.indelLength = indelLength;
        nonSequencingErrorProb = 1 - (1 - artifactProb) * (1 - nonSomaticProb);
    }

    public double getTumorLogOdds() { return tumorLogOdds; }

    public double getArtifactProb() { return artifactProb; }

    public double getNonSequencingErrorProb() { return nonSequencingErrorProb; }

    public int getAltCount() { return altCount; }

    public int getTotalCount() { return totalCount; }

    public int getIndelLength() { return indelLength; }
}
