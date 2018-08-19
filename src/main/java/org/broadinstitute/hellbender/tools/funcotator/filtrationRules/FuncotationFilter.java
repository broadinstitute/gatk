package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.Map;

/**
 * A filter to apply to Funcotations in {@link org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations}.
 *
 * Filters can define an arbitrary number of rules which must match on the Funcotations of a variant in order
 * for that variant to "pass". Passing variants will be annotated with the filter's name in the output VCF.
 */
public abstract class FuncotationFilter {

    /**
     * The INFO annotation value which should be added to all variants which pass this filter.
     */
    private final String filterName;

    FuncotationFilter(final String filterName) {
        this.filterName = filterName;
    }

    public String getFilterName() {
        return filterName;
    }

    /**
     * Check all of this filter's rules against a set of Funcotations.
     *
     * @param prunedTranscriptFuncotations Funcotation values of a single transcript. Assumed to have
     *                                     been "pruned" to remove null / empty values. Never {@code null}
     * @return true if the Funcotations match all of this filter's rules, and false otherwise
     */
    public Boolean checkFilter(final Map<String, String> prunedTranscriptFuncotations) {
        Utils.nonNull(prunedTranscriptFuncotations);

        return getRules().stream()
                .map(rule -> rule.checkRule(prunedTranscriptFuncotations))
                .reduce(Boolean::logicalAnd)
                .orElse(false);
    }

    /**
     * Build the collection of rules which must match to pass this filter.
     */
    abstract List<FuncotationFiltrationRule> getRules();
}
