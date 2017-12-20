package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.BreakpointPair;

public final class SVGraphEdgeEvidence {
    private final EvidenceTargetLink evidenceTargetLink;
    private final VariantContext calledVariant;
    private final BreakpointPair breakpointPair;

    public SVGraphEdgeEvidence() {
        this.evidenceTargetLink = null;
        this.calledVariant = null;
        this.breakpointPair = null;
    }

    public SVGraphEdgeEvidence(final EvidenceTargetLink evidenceTargetLink) {
        this.evidenceTargetLink = evidenceTargetLink;
        this.calledVariant = null;
        this.breakpointPair = null;
    }

    public SVGraphEdgeEvidence(final VariantContext calledVariant) {
        this.evidenceTargetLink = null;
        this.calledVariant = calledVariant;
        this.breakpointPair = null;
    }

    public SVGraphEdgeEvidence(final BreakpointPair breakpointPair) {
        this.evidenceTargetLink = null;
        this.calledVariant = null;
        this.breakpointPair = breakpointPair;
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
}
