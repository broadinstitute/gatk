package org.broadinstitute.hellbender.engine.filters;

import java.io.Serializable;

/**
 * Collects common variant filters.
 */
public final class VariantFilterLibrary {
    public static VariantFilter ALLOW_ALL_VARIANTS = (VariantFilter & Serializable) variant -> true;
}
