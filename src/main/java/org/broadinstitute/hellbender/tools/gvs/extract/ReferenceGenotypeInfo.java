package org.broadinstitute.hellbender.tools.gvs.extract;

public class ReferenceGenotypeInfo {
    private String sampleName;
    private int GQ;

    public ReferenceGenotypeInfo(String sampleName, int GQ) {
        this.sampleName = sampleName;
        this.GQ = GQ;
    }

    public String getSampleName() {
        return sampleName;
    }

    public int getGQ() {
        return GQ;
    }
}