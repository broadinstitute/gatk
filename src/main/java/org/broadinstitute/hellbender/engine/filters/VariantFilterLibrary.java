package org.broadinstitute.hellbender.engine.filters;

/**
 * Collects common variant filters.
 */
public final class VariantFilterLibrary {
    public static VariantFilter ALLOW_ALL_VARIANTS = variant -> true;
}
