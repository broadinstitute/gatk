package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.Map;

/**
 * A rule to match against a variant within a {@link FuncotationFilter}.
 */
interface FuncotationFiltrationRule {

    /**
     * Check if a variant and the Funcotations for one of its transcripts match this rule.
     */
    boolean checkRule(final VariantContext variant, final Map<String, String> prunedTranscriptFuncotations);
}
