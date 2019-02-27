package org.broadinstitute.hellbender.engine.filters;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * Collects common variant filters.
 */
public final class VariantFilterLibrary {
    public static VariantFilter ALLOW_ALL_VARIANTS = new AllowAllVariantsVariantFilter();
    public static VariantFilter NOT_SV_OR_SYMBOLIC = new NotSymbolicOrSVVariantFilter();
    public static VariantFilter PASSES_FILTERS = new PassesFiltersVariantFilter();

    /** Do not filter out any variants. */
    public static class AllowAllVariantsVariantFilter implements VariantFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final VariantContext variant) { return true; }
    }

    /** Filter out any variants that are symbolic or SV. */
    public static class NotSymbolicOrSVVariantFilter implements VariantFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final VariantContext variant) {
            return !variant.isSymbolicOrSV();
        }
    }

    /** Filter out any variants that fail (variant-level) filters. */
    public static class PassesFiltersVariantFilter implements VariantFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final VariantContext variant) {
            return !variant.isFiltered();
        }
    }

}
