package org.broadinstitute.hellbender.tools.funcotator;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * TODO: Docs
 *
 * Supports ordering of the transcripts.
 */
public class FuncotationMap {

    final private Map<String, Map<String, String>> txToFuncotationNameToFuncotationValue = new LinkedHashMap<>();

    private FuncotationMap() {}

    /** TODO: Update docs
     * Given a VCF header and an annotation line, return a map of name to funcotation.
     *
     * @param funcotationHeaderKeys The funcotation keys from the description of the funcotation info field.
     *                              @see FuncotatorUtils#extractFuncotatorKeysFromHeaderDescription(String).  Never {@code null}
     * @param funcotationAttributeAsString The value in a funcotation info field.  Never {@code null}
     * @return map of the fields specified in the description to the values in the attribute.  Never {@code null}
     */
    public void add(final String transcriptIdFuncotationName, final String[] funcotationHeaderKeys, final String funcotationAttributeAsString) {
        Utils.nonNull(funcotationHeaderKeys);
        Utils.nonNull(funcotationAttributeAsString);

        final String[] valuesSplit = StringUtils.splitByWholeSeparatorPreserveAllTokens(funcotationAttributeAsString, "|");
        final String[] values = funcotationAttributeAsString.startsWith("|") ? ArrayUtils.subarray(valuesSplit, 1, valuesSplit.length): valuesSplit;

        if (funcotationHeaderKeys.length != values.length) {
            throw new GATKException.ShouldNeverReachHereException("Could not parse FUNCOTATION field properly.");
        }
        final Map<String, String> funcotationNameToFuncotationValue = IntStream.range(0, funcotationHeaderKeys.length).boxed()
                .collect(Collectors.toMap(i -> funcotationHeaderKeys[i], i -> values[i]));
        final String transcriptId = funcotationNameToFuncotationValue.get(transcriptIdFuncotationName);
        txToFuncotationNameToFuncotationValue.put(transcriptId, funcotationNameToFuncotationValue);
    }

    /**
     * TODO: Docs
     * @param transcriptId
     * @param funcotationName
     * @return {@code null} if no funcotation found.
     */
    public String get(final String transcriptId, final String funcotationName) {
        return txToFuncotationNameToFuncotationValue.getOrDefault(transcriptId, new LinkedHashMap<>()).get(funcotationName);
    }

    /** TODO: Docs
     *
     * @param transcriptId
     * @return
     */
    public Map<String, String> getByTranscript(final String transcriptId){
        return txToFuncotationNameToFuncotationValue.get(transcriptId);
    }

    /**
     * TODO: Docs
     * @param funcotationName
     * @return Map of transcript ID to the value for the funcotation.
     */
    public Map<String, String> getByFuncotationName(final String funcotationName) {
        return txToFuncotationNameToFuncotationValue.entrySet().stream()
                .collect(Collectors.toMap(e -> e.getKey(), e -> e.getValue().get(funcotationName)));
    }

    /** TODO: Docs
     *
     * @param transcriptIdFuncotationName
     * @param funcotationHeaderKeys
     * @param funcotationAttributeAsString
     * @return
     */
    public static FuncotationMap create(final String transcriptIdFuncotationName, final String[] funcotationHeaderKeys, final String funcotationAttributeAsString) {
        final FuncotationMap result = createEmpty();
        result.add(transcriptIdFuncotationName, funcotationHeaderKeys, funcotationAttributeAsString);
        return result;
    }

    /** TODO: Docs
     *
     * @return
     */
    public static FuncotationMap createEmpty() {
        return new FuncotationMap();
    }

    public List<String> keyList() {
        return new ArrayList<>(txToFuncotationNameToFuncotationValue.keySet());
    }
}
