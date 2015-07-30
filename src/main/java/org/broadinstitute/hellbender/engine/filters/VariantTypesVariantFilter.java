package org.broadinstitute.hellbender.engine.filters;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Set;

/**
 * Keep only variants with any of these variant types.
 */
public final class VariantTypesVariantFilter implements VariantFilter {
    private static final long serialVersionUID = 1L;

    private final Set<VariantContext.Type> sampleTypes;

    public VariantTypesVariantFilter(Set<VariantContext.Type> includeTypes) {
        Utils.nonNull(includeTypes);
        sampleTypes = includeTypes;
    }

    @Override
    public boolean test(final VariantContext vc) {
        final VariantContext.Type vcSampleType = vc.getType();
        return sampleTypes.contains(vcSampleType);
    }
}
