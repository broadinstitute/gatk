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
    private final boolean ignoreNonRef;

    public VariantTypesVariantFilter(Set<VariantContext.Type> includeTypes) {
        this(includeTypes, false);
    }

    public VariantTypesVariantFilter(Set<VariantContext.Type> includeTypes, final boolean ignoreNonRefAlleles) {
        Utils.nonNull(includeTypes);
        sampleTypes = includeTypes;
        ignoreNonRef = ignoreNonRefAlleles;
    }

    @Override
    public boolean test(final VariantContext vc) {
        final VariantContext.Type vcSampleType = vc.getType(ignoreNonRef);
        return sampleTypes.contains(vcSampleType);
    }
}
