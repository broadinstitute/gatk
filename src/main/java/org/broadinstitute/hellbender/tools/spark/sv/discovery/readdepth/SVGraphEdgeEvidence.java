package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.BreakpointPair;

public final class SVGraphEdgeEvidence {
    private final EvidenceTargetLink evidenceTargetLink;
    private final VariantContext calledVariant;
    private final BreakpointPair breakpointPair;
    private final EventRecord eventRecord;

    public SVGraphEdgeEvidence() {
        this.evidenceTargetLink = null;
        this.calledVariant = null;
        this.breakpointPair = null;
        this.eventRecord = null;
    }

    public SVGraphEdgeEvidence(final EvidenceTargetLink evidenceTargetLink) {
        this.evidenceTargetLink = evidenceTargetLink;
        this.calledVariant = null;
        this.breakpointPair = null;
        this.eventRecord = null;
    }

    public SVGraphEdgeEvidence(final VariantContext calledVariant) {
        this.evidenceTargetLink = null;
        this.calledVariant = calledVariant;
        this.breakpointPair = null;
        this.eventRecord = null;
    }

    public SVGraphEdgeEvidence(final BreakpointPair breakpointPair) {
        this.evidenceTargetLink = null;
        this.calledVariant = null;
        this.breakpointPair = breakpointPair;
        this.eventRecord = null;
    }

    public SVGraphEdgeEvidence(final EventRecord eventRecord) {
        this.evidenceTargetLink = null;
        this.calledVariant = null;
        this.breakpointPair = null;
        this.eventRecord = eventRecord;
    }

    public EvidenceTargetLink getEvidenceTargetLink() {
        return evidenceTargetLink;
    }

    public VariantContext getCalledVariant() {
        return calledVariant;
    }

    public BreakpointPair getBreakpointPair() {
        return breakpointPair;
    }

    public EventRecord getEventRecord() {
        return eventRecord;
    }

    public String getId() {
        if (evidenceTargetLink != null) return evidenceTargetLink.toString();
        if (calledVariant != null) return calledVariant.getID();
        if (breakpointPair != null) return breakpointPair.toString();
        if (eventRecord != null) return eventRecord.getId();
        return "NA";
    }
}
