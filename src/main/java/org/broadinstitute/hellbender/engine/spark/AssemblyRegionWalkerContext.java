package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

/**
 * Encapsulates an {@link AssemblyRegion} with its {@link ReferenceContext} and {@link FeatureContext}.
 */
public final class AssemblyRegionWalkerContext {
    private final AssemblyRegion assemblyRegion;
    private final ReferenceContext referenceContext;
    private final FeatureContext featureContext;

    public AssemblyRegionWalkerContext(AssemblyRegion assemblyRegion, ReferenceContext referenceContext, FeatureContext featureContext) {
        this.assemblyRegion = assemblyRegion;
        this.referenceContext = referenceContext;
        this.featureContext = featureContext;
    }

    public AssemblyRegion getAssemblyRegion() {
        return assemblyRegion;
    }

    public ReferenceContext getReferenceContext() {
        return referenceContext;
    }

    public FeatureContext getFeatureContext() {
        return featureContext;
    }
}
