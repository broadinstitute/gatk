package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

public abstract class Mutect2VariantFilter extends Mutect2Filter {
    public Mutect2VariantFilter() { }

    public double errorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        final double result = requiredAnnotations().stream().allMatch(vc::hasAttribute) ? calculateErrorProbability(vc, filteringEngine, referenceContext) : 0;
        return Mutect2FilteringEngine.roundFinitePrecisionErrors(result);
    }

    protected abstract double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext);
}
