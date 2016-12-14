package org.broadinstitute.hellbender.engine.spark;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

/**
 * Encapsulates a {@link VariantContext} with the reads that overlap it (the {@link ReadsContext} and
 * its {@link ReferenceContext} and {@link FeatureContext}.
 */
public class VariantWalkerContext {
    private final VariantContext variant;
    private final ReadsContext readsContext;
    private final ReferenceContext referenceContext;
    private final FeatureContext featureContext;

    public VariantWalkerContext(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        this.variant = variant;
        this.readsContext = readsContext;
        this.referenceContext = referenceContext;
        this.featureContext = featureContext;
    }

    public VariantContext getVariant() {
        return variant;
    }

    public ReadsContext getReadsContext() {
        return readsContext;
    }

    public ReferenceContext getReferenceContext() {
        return referenceContext;
    }

    public FeatureContext getFeatureContext() {
        return featureContext;
    }
}
