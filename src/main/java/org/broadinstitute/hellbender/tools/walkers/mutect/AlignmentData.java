package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

public class AlignmentData {

    private AlignmentContext alignmentContext;
    private ReferenceContext referenceContext;

    public AlignmentData (AlignmentContext alignmentContext, ReferenceContext referenceContext)
    {
        this.alignmentContext = alignmentContext;
        this.referenceContext = referenceContext;
    }

    public AlignmentContext getAlignmentContext() {
        return alignmentContext;
    }

    public ReferenceContext getReferenceContext() {
        return referenceContext;
    }
}
