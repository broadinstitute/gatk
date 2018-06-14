package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


/**
 * A linked map of transcript IDs to funcotations.  Also, supports some querying.
 *
 * Supports ordering of the transcripts.
 *
 * Multiple GencodeFuncotations for the same transcript are prohibited.
 *
 * Not thread-safe.
 */
public class FuncotationMap {

    public final static String NO_TRANSCRIPT_AVAILABLE_KEY = "no_transcript";

    /** Standard Logger.  */
    protected static final Logger logger = LogManager.getLogger(FuncotationMap.class);


    final private Map<String, LinkedHashSet<Funcotation>> txToFuncotations = new LinkedHashMap<>();

    private FuncotationMap() {}


    /**
     * @param transcriptId the specified transcript ID.  Use {@see NO_TRANSCRIPT_AVAILABLE_KEY} if there are no transcripts.  Never {@code null}
     * @return A list of the Gencode Funcotations only. Empty list, if nothing found.  Never {@code null}
     */
    public List<GencodeFuncotation> getGencodeFuncotations(final String transcriptId) {
        Utils.nonNull(transcriptId);
        return txToFuncotations.getOrDefault(transcriptId, new LinkedHashSet<>()).stream()
                .filter(FuncotatorUtils::isGencodeFuncotation).map(f-> (GencodeFuncotation) f).collect(Collectors.toList());
    }

    /**
     * Get all funcotations for the given transcript.  Note that the alleles will be mixed.
     * @param transcriptId the specified transcript ID.  Use {@see NO_TRANSCRIPT_AVAILABLE_KEY} if there are no transcripts.
     * @return Empty list, if nothing found.  Never {@code null}
     */
    public List<Funcotation> get(final String transcriptId) {
        Utils.nonNull(transcriptId);
        return new ArrayList<>(txToFuncotations.getOrDefault(transcriptId, new LinkedHashSet<>()));
    }

    /**
     * @param transcriptId the specified transcript ID.  Use {@see NO_TRANSCRIPT_AVAILABLE_KEY} if there are no transcripts.  Never {@code null}
     * @param fieldName The field name to search.  Never {@code null}
     * @param allele Only return fields from funcotations with the specified allele.  Never {@code null}
     * @return Value of the given field for the transcript ID and allele.  Return {@code null} if field not found.
     */
    public String getFieldValue(final String transcriptId, final String fieldName, final Allele allele) {
        Utils.nonNull(transcriptId);
        Utils.nonNull(fieldName);
        Utils.nonNull(allele);
        final List<String> values = txToFuncotations.getOrDefault(transcriptId, new LinkedHashSet<>()).stream()
                .filter(f -> f.hasField(fieldName))
                .filter(f -> f.getAltAllele().equals(allele))
                .map(f -> f.getField(fieldName))
                .collect(Collectors.toList());
        if (values.size() > 1) {
            throw new GATKException.ShouldNeverReachHereException("Found more than one value for " + transcriptId + ", "
                    + allele + ", " + fieldName);
        }
        if (values.size() == 0) {
            return null;
        } else {
            return values.get(0);
        }
    }

    /**
     * @return An empty FuncotationMap.  Never {@code null}
     */
    private static FuncotationMap createEmpty() {
        return new FuncotationMap();
    }

    /**
     * @return A FuncotationMap with only one transcript ID, which is {@see NO_TRANSCRIPT_AVAILABLE_KEY}.  Never {@code null}
     */
    public static FuncotationMap createNoTranscriptInfo(final List<Funcotation> funcotations) {
        final FuncotationMap result = createEmpty();
        result.add(NO_TRANSCRIPT_AVAILABLE_KEY, funcotations);
        return result;
    }

    /**  This method checks for transcript duplications in the list of GencodeFuncotations.
     * Recommended to gather all gencode funcotations before calling this method.
     * @param gencodeFuncotations  Gencode funcotations, each for a unique transcript.  If duplicate transcript IDs are
     *                             found, error is thrown.  Never {@code null}
     * @return A new FuncotationMap created only from the given list of GencodeFuncotations.  Never {@code null}
     */
    public static FuncotationMap createFromGencodeFuncotations(final List<GencodeFuncotation> gencodeFuncotations) {

        Utils.nonNull(gencodeFuncotations);
        Utils.validateArg(!areDuplicateTranscriptIDsFound(gencodeFuncotations), "Duplicate transcript ID entries were found in input: " +
                gencodeFuncotations.stream().map(gf -> gf.getAnnotationTranscript()).collect(Collectors.joining(",")));
        final FuncotationMap result = createEmpty();
        gencodeFuncotations.forEach(f -> result.addWithoutGencodeCheck(f.getAnnotationTranscript(), f));
        return result;
    }

    /**
     * Add the given funcotations to the given transcript ID.
     * TODO: See https://github.com/broadinstitute/gatk/issues/4850 since the enforcement/need of the No-Gencode Funcotations rule is a hack.
     * @param txId  Never {@code null}
     * @param funcotations  Cannot have any Gencode funcotations or error is thrown.  Never {@code null}
     */
    public void add(final String txId, final List<Funcotation> funcotations) {
        Utils.nonNull(txId);
        Utils.nonNull(funcotations);
        if (FuncotatorUtils.areAnyGencodeFuncotation(funcotations)) {
            throw new GATKException.ShouldNeverReachHereException( "At this time, a Gencode Funcotation cannot be added to a FuncotationMap.  If you see this error message, please contact the GATK dev team with a forum post.");
        }
        addWithoutGencodeCheck(txId, funcotations);
    }

    private void addWithoutGencodeCheck(final String txId, final List<Funcotation> funcotations) {
        final LinkedHashSet<Funcotation> existingFuncotationsToUpdate = txToFuncotations.getOrDefault(txId, new LinkedHashSet<>());
        existingFuncotationsToUpdate.addAll(funcotations);
        txToFuncotations.put(txId, existingFuncotationsToUpdate);
    }

    private void addWithoutGencodeCheck(final String txId, final Funcotation funcotation) {
        final LinkedHashSet<Funcotation> existingFuncotationsToUpdate = txToFuncotations.getOrDefault(txId, new LinkedHashSet<>());
        existingFuncotationsToUpdate.add(funcotation);
        txToFuncotations.put(txId, existingFuncotationsToUpdate);
    }

    /** Add the given funcotation to the given transcript ID.
     * TODO: https://github.com/broadinstitute/gatk/issues/4850 since the enforcement/need of the No-Gencode Funcotations rule is a hack.
     * @param txId Never {@code null}
     * @param funcotation Cannot be a Gencode funcotation or error is thrown.  Never {@code null}
     */
    public void add(final String txId, final Funcotation funcotation) {
        Utils.nonNull(txId);
        Utils.nonNull(funcotation);
        if (FuncotatorUtils.isGencodeFuncotation(funcotation)) {
            throw new GATKException.ShouldNeverReachHereException( "At this time, a Gencode Funcotation cannot be added to a FuncotationMap.  If you see this error message, please contact the GATK dev team with a forum post.");
        }
        addWithoutGencodeCheck(txId, funcotation);
    }


    /**
     * Get the list of transcripts in order.
     *
     * @return Never {@code null}
     */
    public List<String> getTranscriptList() {
        return new ArrayList<>(txToFuncotations.keySet());
    }

    /**
     * Get the transcripts as a of transcripts in order.
     *
     * @return Never {@code null}
     */
    public LinkedHashSet<String> getTranscriptSet() {
        return new LinkedHashSet<>(txToFuncotations.keySet());
    }

    /** Create a FuncotationMap where all Funcotations will be TableFuncotations.  This is useful for parsing existing
     *   VCFs.
     * Only renders for a single allele.
     * See {@link FuncotatorUtils#extractFuncotatorKeysFromHeaderDescription(String)} for getting the funcotation keys.
     *
     * @param transcriptFieldName The field name to use for transcript IDs.  Use {@see NO_TRANSCRIPT_AVAILABLE_KEY} if unknown.
     *                            If not in the funcotation keys, then the Funcotation map will be created with one transcript ID, {@see NO_TRANSCRIPT_AVAILABLE_KEY}
     *                            Never {@code null}
     * @param funcotationKeys The ordered keys of the funcotation field.  Never {@code null}
     * @param funcotationAttributeForSingleAllele  The funcotation attribute from a VCF, split for a single Allele.  Never {@code null}
     * @param altAllele The alternate allele for the created funcotations.  Never {@code null}
     * @param datasourceName The datasource name to use for all of the created funcotatinos.  Never {@code null}
     * @return a funcotation map.  Note that no funcotations will be GencodeFuncotations.  Never {@code null}
     */
    public static FuncotationMap createAsAllTableFuncotationsFromVcf(final String transcriptFieldName, final String[] funcotationKeys,
                                                                     final String funcotationAttributeForSingleAllele, final Allele altAllele,
                                                                     final String datasourceName) {
        Utils.nonNull(transcriptFieldName);
        Utils.nonNull(funcotationKeys);
        Utils.nonNull(funcotationAttributeForSingleAllele);
        Utils.nonNull(altAllele);
        Utils.nonNull(datasourceName);

        final FuncotationMap result = createEmpty();
        final String[] funcotationAttributeForSingleAlleleByTranscript = StringUtils.splitByWholeSeparator(funcotationAttributeForSingleAllele, VcfOutputRenderer.END_TRANSCRIPT_DELIMITER +
                VcfOutputRenderer.ALL_TRANSCRIPT_DELIMITER + VcfOutputRenderer.START_TRANSCRIPT_DELIMITER);

        for (final String funcotationAttribute : funcotationAttributeForSingleAlleleByTranscript) {
            final String[] values = StringUtils.splitByWholeSeparatorPreserveAllTokens(funcotationAttribute, VcfOutputRenderer.FIELD_DELIMITER);
            if (values[0].startsWith(VcfOutputRenderer.START_TRANSCRIPT_DELIMITER)) {
                values[0] = values[0].replace(VcfOutputRenderer.START_TRANSCRIPT_DELIMITER, "");
            }
            if (values[values.length - 1].endsWith(VcfOutputRenderer.END_TRANSCRIPT_DELIMITER)) {
                values[values.length - 1] = values[values.length - 1].replace(VcfOutputRenderer.END_TRANSCRIPT_DELIMITER, "");
            }
            if (values.length != funcotationKeys.length) {
                logger.error("Keys:  " + StringUtils.join(funcotationKeys, ", "));
                logger.error("Values:  " + StringUtils.join(values, ", "));
                throw new GATKException.ShouldNeverReachHereException("Cannot parse the funcotation attribute.  Num values: " + values.length + "   Num keys: " + funcotationKeys.length);
            }
            final Map<String, String> simpleNameValuePairs = IntStream.range(0, values.length).boxed().collect(Collectors
                    .toMap(i -> funcotationKeys[i], i-> values[i]));

            final List<String> valuesAsList = Arrays.asList(funcotationKeys).stream().map(k -> simpleNameValuePairs.get(k)).collect(Collectors.toList());
            result.add(simpleNameValuePairs.getOrDefault(transcriptFieldName, NO_TRANSCRIPT_AVAILABLE_KEY), TableFuncotation.create(Arrays.asList(funcotationKeys), valuesAsList, altAllele, datasourceName, null));
        }
        return result;
    }

    private static boolean areDuplicateTranscriptIDsFound(final List<GencodeFuncotation> gencodeFuncotations) {
        return gencodeFuncotations.size() != new HashSet<>(gencodeFuncotations).size();
    }
}
