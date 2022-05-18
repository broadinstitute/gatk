package org.broadinstitute.hellbender.tools.gvs.extract;

public class SpanningDeletionRecord extends ReferenceRecord {
    private final String gt;
    private final String gq;

    public SpanningDeletionRecord(long location, long sampleId, int length, String gt, String gq) {
        super(location, sampleId, length, "*");
        this.gt = gt;
        this.gq = gq;
    }

    public String getGT() {
        return gt;
    }

    public String getGQ() {
        return gq;
    }
}
