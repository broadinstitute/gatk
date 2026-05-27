package org.broadinstitute.hellbender.tools.gvs.extract;

public class ReferenceGenotypeInfo {
    private final String sampleName;
    private final int sampleId;
    private final int GQ;
    private final Integer dp;

    public ReferenceGenotypeInfo(String sampleName, int GQ, int sampleId, Integer dp) {
        this.sampleName = sampleName;
        this.GQ = GQ;
        this.sampleId = sampleId;
        this.dp = dp;
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

    public Integer getDp() {
        return dp;
    }
}