package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequence;

public final class SAMRecordAndReference {
    private final SAMRecord samRec;
    private final ReferenceSequence refSeq;

    public SAMRecordAndReference(final SAMRecord samRec, final ReferenceSequence refSeq) {
        this.samRec = samRec;
        this.refSeq = refSeq;
    }

    public SAMRecord getSamRecord() {
        return samRec;
    }

    public ReferenceSequence getReferenceSequence() {
        return refSeq;
    }
}
