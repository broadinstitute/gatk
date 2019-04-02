package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.Map;
import java.util.Set;

/**
 * A rule to match against the Funcotations from a variant within a {@link FuncotationFilter}.
 */
interface FuncotationFiltrationRule {

    /**
     * Check if a set of Funcotations matches this rule.
     */
    boolean checkRule(final Set<Map.Entry<String, String>> funcotations, VariantContext variant);
}
