package org.broadinstitute.hellbender.tools.funcotator;


import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Class for managing aliases and querying {@link Funcotation} to determine fields.
 *
 * Typically, you would use this to map output columns to input funcotation fields for tsv generation.
 */
public class AliasProvider {
    /**
     * Maps a field name to a list of alternate field names.
     *
     * For example, "CONTIG" -> ["CONTIG", "chr", "Chromsome", "contig"]
     */
    final private LinkedHashMap<String, List<String>> aliasMap;

    /**
     *
     * @param aliasMap key name to list of potential alternate key names.  No need to include the key name itself in
     *                 the value list.
     */
    public AliasProvider(final LinkedHashMap<String, List<String>> aliasMap) {
        Utils.nonNull(aliasMap);

        // Take all of the aliases and prepend the key to list of values.
        this.aliasMap = aliasMap.entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,e-> createAliasList(aliasMap, e.getKey()),
                        (x1, x2) -> {
                            throw new IllegalArgumentException("Should not be able to have duplicate field names."); },
                        LinkedHashMap::new ));
    }

    // Make sure that the column name itself is at the front of the alias list.
    private static List<String> createAliasList(final LinkedHashMap<String, List<String>> columnNameToAliasMap, final String columnName) {
        final List<String> result = new ArrayList<>();
        result.add(columnName);
        result.addAll(columnNameToAliasMap.get(columnName));
        return result;
    }

    /**
     * Search the given funcotations (derived from the funoctation map and transcript key) for the highest priority alias that
     *  appears for the given field name.
     *
     * Notes:
     * - This method assumes that the funcotation map has the same fields regardless of allele.
     *
     * @param fieldName field name for which to find an alias in the given funcotations.  Never {@code null}
     * @param txToFuncotationMap {@link FuncotationMap} to search for field names.  Never {@code null}
     * @param txId transcript ID to use in the search.   Never {@code null}
     * @return the highest priority alias field name found in the funcotation map.
     *  Returns empty string ("") if nothing is found.
     */
    private String findFieldNameInFuncotations(final String fieldName, final FuncotationMap txToFuncotationMap, final String txId) {
        Utils.nonNull(fieldName);
        Utils.nonNull(txToFuncotationMap);
        Utils.nonNull(txId);

        if (!aliasMap.keySet().contains(fieldName)) {
            return "";
        }
        final List<String> candidateFieldNames = aliasMap.get(fieldName);
        final Set<String> funcotationFieldNames = txToFuncotationMap.getFieldNames(txId);
        return candidateFieldNames.stream()
                .filter(funcotationFieldNames::contains)
                .findFirst().orElse("");
    }

    /**
     * Create a mapping that links the fields of this alias provider to the funcotation fields in the given
     *  {@link FuncotationMap}.
     *
     * @param txToFuncotationMap Never {@code null}
     * @param txId Transcript ID.  Never {@code null}
     * @return A mapping of the fields in this alias provider to the funcotation fields contained in the FuncotationMap,
     * transcript ID pairing.  The keys will be the same keys from the alias map that this instance was initialized with.
     * Value for a given key will be empty string ("") if no matching funcotation field exists.
     * Never {@code null}
     */
    public LinkedHashMap<String, String> createColumnNameToFieldNameMap(final FuncotationMap txToFuncotationMap, final String txId) {
        return aliasMap.keySet().stream()
                .map(k -> Pair.of(k, findFieldNameInFuncotations(k, txToFuncotationMap, txId)))
                .collect(Collectors.toMap(Pair::getLeft, Pair::getRight,
                        (x1, x2) -> {
                            throw new IllegalArgumentException("Should not be able to have duplicate field names."); },
                        LinkedHashMap::new ));
    }

    /**
     * Get the fields that have aliases.
     *
     * @return Never {@code null}
     */
    public List<String> getFields() {
        return new ArrayList<>(aliasMap.keySet());
    }
}
