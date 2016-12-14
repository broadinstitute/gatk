package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Encapsulates an {@link GATKRead} with its {@link ReferenceContext} and {@link FeatureContext}.
 */
public class ReadWalkerContext {
    private final GATKRead read;
    private final ReferenceContext referenceContext;
    private final FeatureContext featureContext;

    public ReadWalkerContext(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        this.read = read;
        this.referenceContext = referenceContext;
        this.featureContext = featureContext;
    }

    public GATKRead getRead() {
        return read;
    }

    public ReferenceContext getReferenceContext() {
        return referenceContext;
    }

    public FeatureContext getFeatureContext() {
        return featureContext;
    }
}
