package org.broadinstitute.hellbender.tools.gvs.extract;

public class ReferenceGenotypeInfo {
    private final String sampleName;
    private final int sampleId;
    private final int GQ;

    public ReferenceGenotypeInfo(String sampleName, int GQ, int sampleId) {
        this.sampleName = sampleName;
        this.GQ = GQ;
        this.sampleId = sampleId;
    }

    public String getSampleName() {
        return sampleName;
    }

    public int getGQ() {
        return GQ;
    }

    public int getSampleId() {
        return sampleId;
    }
}