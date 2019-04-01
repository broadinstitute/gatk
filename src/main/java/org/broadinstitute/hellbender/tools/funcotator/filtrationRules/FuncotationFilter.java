package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import htsjdk.variant.variantcontext.VariantContext;
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
     * @param variant VariantContext of this transcript.
     *
     * @return true if the Funcotations match all of this filter's rules, and false otherwise
     */
    public Boolean checkFilter(final Set<Map.Entry<String, String>> prunedTranscriptFuncotations, final VariantContext variant) {
        Utils.nonNull(prunedTranscriptFuncotations);

        return getRules().stream()
                .map(rule -> rule.checkRule(prunedTranscriptFuncotations, variant))
                .reduce(Boolean::logicalAnd)
                .orElse(false);
    }

    /**
     * Build the collection of rules which must match to pass this filter.
     */
    abstract List<FuncotationFiltrationRule> getRules();


    /**
     * Given a set of extracted funcotations of a variant, return a Stream of values that belong to the key.
     * @param funcotations A Set of Map.Entry key-value pairs of funcotations
     * @param key The funcotations to get
     * @param defaultValue A default value to get if there are no matches
     * @return A Stream of matched values
     */
    protected Stream<String> matchOnKeyOrDefault(Set<Map.Entry<String, String>> funcotations, String key, String defaultValue) {
        return getMatchesOrDefault(funcotations, entry -> entry.getKey().equals(key), defaultValue);
    }

    /**
     * Given a set of extracted funcotations, return a Stream of values whose entry matches the supplied predicate
     * @param funcotations A Set of Map.Entry key-value pairs of funcotations
     * @param matcher A predicate taking a Map.Entry to filter on
     * @param defaultValue A default value to get if there are no matches
     * @return A Stream of matched values
     */
    protected Stream<String> getMatchesOrDefault(Set<Map.Entry<String, String>> funcotations, Predicate<Map.Entry<String, String>> matcher, String defaultValue) {
        Set<String> matched = funcotations.stream().filter(matcher).map(Map.Entry::getValue).collect(Collectors.toSet());
        if (matched.isEmpty()) {
            return Stream.of(defaultValue);
        } else {
            return matched.stream();
        }
    }

    /**
     * Does the set of funcotations contain the key supplied
     * @param funcotations A Set of Map.Entry key-value pairs of funcotations
     * @param key The key to look for
     * @return Whether or not the key exists in the funcotations
     */
    protected Boolean containsKey(Set<Map.Entry<String, String>> funcotations, String key) {
        return funcotations.stream().anyMatch(entry -> entry.getKey().equals(key));
    }

}
