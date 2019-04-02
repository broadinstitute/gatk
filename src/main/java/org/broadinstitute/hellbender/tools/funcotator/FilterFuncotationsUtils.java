package org.broadinstitute.hellbender.tools.funcotator;

import java.util.AbstractMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class FilterFuncotationsUtils {
    /**
     * Munge a FuncotationMap into a stream of sets of pairs of funcotation name / funcotation value.
     */
    public static Stream<Set<Map.Entry<String, String>>> getTranscriptFuncotations(final FuncotationMap funcotationMap) {
        return funcotationMap.getTranscriptList().stream()
                .map(funcotationMap::get)
                .map(funcotations -> funcotations.stream()
                        .flatMap(FilterFuncotationsUtils::extractFuncotationFields)
                        .filter(entry -> entry.getValue() != null && !entry.getValue().isEmpty())
                        .collect(Collectors.toSet()));
    }

    /**
     * Parse the entries in a Funcotation into a stream of map entries.
     */
    private static Stream<Map.Entry<String, String>> extractFuncotationFields(final Funcotation funcotation) {
        return funcotation.getFieldNames().stream()
                .map(name -> new AbstractMap.SimpleEntry<>(name, funcotation.getField(name)));
    }
}
