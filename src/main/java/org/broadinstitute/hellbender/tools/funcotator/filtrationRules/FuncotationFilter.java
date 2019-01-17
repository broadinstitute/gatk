package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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
    public Boolean checkFilter(final Set<Map.Entry<String, String>> prunedTranscriptFuncotations) {
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

    protected Stream<String> matchOnKeyOrDefault(Set<Map.Entry<String, String>> funcotations, String key, String defaultValue) {
        return getMatchesOrDefault(funcotations, entry -> entry.getKey().equals(key), defaultValue);
    }

    protected Stream<String> getMatchesOrDefault(Set<Map.Entry<String, String>> funcotations, Predicate<Map.Entry<String, String>> matcher, String defaultValue) {
        Set<String> matched = funcotations.stream().filter(matcher).map(Map.Entry::getValue).collect(Collectors.toSet());
        if (matched.isEmpty()) {
            return Stream.of(defaultValue);
        } else {
            return matched.stream();
        }
    }

    protected Boolean containsKey(Set<Map.Entry<String, String>> funcotations, String key) {
        return funcotations.stream().anyMatch(entry -> entry.getKey().equals(key));
    }

}
