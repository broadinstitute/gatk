package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;


import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.Transition;

/** Simple class for holding summary information derived from the orientation bias filter.  Typically, the orientation bias filter will produce a list of these.  Each representing a sample-artifact transition pair.
 *
 */
public class OrientationSampleTransitionSummary {
    private String  sample;
    private Transition artifactMode;
    private Transition artifactModeComplement;
    private double  obq;
    private long  numArtifactMode;
    private long  numArtifactModeFiltered;
    private long  numNotArtifactMode;
    private long  numNonRefPassingVariants;

    public OrientationSampleTransitionSummary(String sample, Transition artifactMode, Transition artifactModeComplement, double obq, long numArtifactMode, long numArtifactModeFiltered, long numNotArtifactMode, long numNonRefPassingVariants) {
        this.sample = sample;
        this.artifactMode = artifactMode;
        this.artifactModeComplement = artifactModeComplement;
        this.obq = obq;
        this.numArtifactMode = numArtifactMode;
        this.numArtifactModeFiltered = numArtifactModeFiltered;
        this.numNotArtifactMode = numNotArtifactMode;
        this.numNonRefPassingVariants = numNonRefPassingVariants;
    }

    public String getSample() {
        return sample;
    }

    public void setSample(String sample) {
        this.sample = sample;
    }

    public Transition getArtifactMode() {
        return artifactMode;
    }

    public void setArtifactMode(Transition artifactMode) {
        this.artifactMode = artifactMode;
    }

    public double getObq() {
        return obq;
    }

    public void setObq(double obq) {
        this.obq = obq;
    }

    public long getNumArtifactMode() {
        return numArtifactMode;
    }

    public void setNumArtifactMode(long numArtifactMode) {
        this.numArtifactMode = numArtifactMode;
    }

    public long getNumArtifactModeFiltered() {
        return numArtifactModeFiltered;
    }

    public void setNumArtifactModeFiltered(long numArtifactModeFiltered) {
        this.numArtifactModeFiltered = numArtifactModeFiltered;
    }

    public long getNumNotArtifactMode() {
        return numNotArtifactMode;
    }

    public void setNumNotArtifactMode(long numNotArtifactMode) {
        this.numNotArtifactMode = numNotArtifactMode;
    }

    public long getNumNonRefPassingVariants() {
        return numNonRefPassingVariants;
    }

    public void setNumNonRefPassingVariants(long numNonRefPassingVariants) {
        this.numNonRefPassingVariants = numNonRefPassingVariants;
    }

    public Transition getArtifactModeComplement() {
        return artifactModeComplement;
    }

    public void setArtifactModeComplement(Transition artifactModeComplement) {
        this.artifactModeComplement = artifactModeComplement;
    }
}
