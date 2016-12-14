package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

/**
 * Encapsulates an {@link AlignmentContext} with its {@link ReferenceContext} and {@link FeatureContext}.
 */
public class LocusWalkerContext {
    private final AlignmentContext alignmentContext;
    private final ReferenceContext referenceContext;
    private final FeatureContext featureContext;

    public LocusWalkerContext(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        this.alignmentContext = alignmentContext;
        this.referenceContext = referenceContext;
        this.featureContext = featureContext;
    }

    public AlignmentContext getAlignmentContext() {
        return alignmentContext;
    }

    public ReferenceContext getReferenceContext() {
        return referenceContext;
    }

    public FeatureContext getFeatureContext() {
        return featureContext;
    }
}
