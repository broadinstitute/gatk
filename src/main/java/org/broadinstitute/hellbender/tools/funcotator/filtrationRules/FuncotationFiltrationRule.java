package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import java.util.Map;

/**
 * A rule to match against the Funcotations from a variant within a {@link FuncotationFilter}.
 */
interface FuncotationFiltrationRule {

    /**
     * Check if a set of Funcotations matches this rule.
     */
    boolean checkRule(final Map<String, String> funcotations);
}
