package org.broadinstitute.hellbender.tools.funcotator;

import com.esotericsoftware.kryo.Kryo;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.spark.GATKRegistrator;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
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
    List<GencodeFuncotation> getGencodeFuncotations(final String transcriptId) {
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
     * @return Value of the given field for the transcript ID and allele.  Return {@code null} if field not found in any
     *  funcotation.  Note that if the funcotations support the given field name, but the variant did not overlap any
     *  records, an empty string will be returned.
     */
    public String getFieldValue(final String transcriptId, final String fieldName, final Allele allele) {
        Utils.nonNull(transcriptId);
        Utils.nonNull(fieldName);
        Utils.nonNull(allele);
        final Set<String> values = txToFuncotations.getOrDefault(transcriptId, new LinkedHashSet<>()).stream()
                .filter(f -> f.hasField(fieldName))
                .filter(f -> f.getAltAllele().equals(allele))
                .map(f -> f.getField(fieldName))
                .collect(Collectors.toSet());
        if (values.size() > 1) {
            throw new UserException.BadInput("Found more than one unique value for the tuple {" + transcriptId + ", "
                    + allele + ", " + fieldName + "}: " + values.stream().collect(Collectors.joining(", ")));
        }
        if (values.size() == 0) {
            return null;
        } else {
            return values.iterator().next();
        }
    }

    /**
     * @return An empty FuncotationMap.  Never {@code null}
     */
    @VisibleForTesting
    public static FuncotationMap createEmpty() {
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
                gencodeFuncotations.stream().map(GencodeFuncotation::getAnnotationTranscript).collect(Collectors.joining(",")));
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

        if (FuncotatorUtils.areAnyGencodeFuncotation(funcotations) && (txToFuncotations.size() > 0)) {
            throw new GATKException.ShouldNeverReachHereException( "At this time, a Gencode Funcotation cannot be added to a non-empty FuncotationMap.  If you see this error message, please contact the GATK dev team with a forum post.");
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
     * Get the list of transcripts, in order.
     *
     * @return Never {@code null}
     */
    public List<String> getTranscriptList() {
        return new ArrayList<>(txToFuncotations.keySet());
    }

    /**
     * Get the transcripts as a set, but preserve the order.
     *
     * @return Never {@code null}
     */
    @VisibleForTesting
    LinkedHashSet<String> getTranscriptSet() {
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

            final List<String> valuesAsList = Arrays.stream(funcotationKeys).map(simpleNameValuePairs::get).collect(Collectors.toList());
            result.add(simpleNameValuePairs.getOrDefault(transcriptFieldName, NO_TRANSCRIPT_AVAILABLE_KEY), TableFuncotation.create(Arrays.asList(funcotationKeys), valuesAsList, altAllele, datasourceName, null));
        }
        return result;
    }

    private static boolean areDuplicateTranscriptIDsFound(final List<GencodeFuncotation> gencodeFuncotations) {
        return gencodeFuncotations.size() != new HashSet<>(gencodeFuncotations).size();
    }

    /**
     * Get all field names found in the funcotations for the given transcript ID.
     *
     * @param transcriptId transcript ID.  Never {@code null}
     * @return Field names in the funcotations.  Empty set if no funcotations (or funcotations have no fields).
     * Never {@code null}
     */
    public Set<String> getFieldNames(final String transcriptId) {
        Utils.nonNull(transcriptId);

        final LinkedHashSet<Funcotation> funcotations =  txToFuncotations.getOrDefault(transcriptId, new LinkedHashSet<>());
        return funcotations.stream().map(Funcotation::getFieldNames).flatMap(LinkedHashSet::stream).collect(Collectors.toSet());
    }

    /**
     * See {@link FuncotationMap#getFieldNames(String)}, but this returns field names for all transcripts and all alleles.
     * @return  See {@link FuncotationMap#getFieldNames(String)}
     */
    @VisibleForTesting
    Set<String> getFieldNames() {
        final List<String> txIds = getTranscriptList();
        final LinkedHashSet<String> result = new LinkedHashSet<>();
        txIds.forEach(txId -> result.addAll(getFieldNames(txId)));
        return result;
    }

    /**
     * See {@link FuncotationMap#getFieldNames(String)}, but this returns field names for a single transcript and a
     *  single allele.
     * @param transcriptId See {@link FuncotationMap#getFieldNames(String)}
     * @param allele Never {@code null}
     * @return See {@link FuncotationMap#getFieldNames(String)}
     */
    private Set<String> getFieldNames(final String transcriptId, final Allele allele) {
        Utils.nonNull(transcriptId);
        Utils.nonNull(allele);

        final LinkedHashSet<Funcotation> funcotations =  txToFuncotations.getOrDefault(transcriptId, new LinkedHashSet<>());
        return funcotations.stream()
                .filter(f -> f.getAltAllele().equals(allele))
                .map(Funcotation::getFieldNames)
                .flatMap(LinkedHashSet::stream).collect(Collectors.toCollection(LinkedHashSet::new));
    }

    /**
     * Get all the alleles in all of the funcotations for a given transcript ID.
     *
     * @param transcriptId Never {@code null}
     * @return a set of alleles that are contained in all funcotations associated with the given transcriptId.  Never {@code null}
     * Will return empty list if there are no funcotations associated with the given transcriptId.
     */
    public Set<Allele> getAlleles(final String transcriptId) {
        final LinkedHashSet<Funcotation> funcotations =  txToFuncotations.getOrDefault(transcriptId, new LinkedHashSet<>());
        return funcotations.stream().map(Funcotation::getAltAllele).collect(Collectors.toSet());
    }

    /**
     * @return whether all transcript-allele combinations have the same fields in the corresponding funcotations.
     */
    public boolean doAllTxAlleleCombinationsHaveTheSameFields() {

        // First get every field seen in this funcotation map.
        final Set<String> allFields = getFieldNames();

        // Then get all txIds
        final List<String> txIds = getTranscriptList();

        final List<Pair<String,Allele>> txAlleleCombos = new ArrayList<>();
        for (final String txId : txIds) {
            getAlleles(txId).forEach(a -> txAlleleCombos.add(Pair.of(txId, a)));
        }

        // For each transcript-allele combo, get the fields
        return txAlleleCombos.stream().allMatch(p -> getFieldNames(p.getLeft(), p.getRight()).equals(allFields));
    }

    /**
     * Copy creation.
     * @param funcotationMap Never {@code null}
     * @return a copy of the input.  Never {@code null}
     */
    public static FuncotationMap create(final FuncotationMap funcotationMap) {
        Utils.nonNull(funcotationMap);
        final Kryo kryo = new Kryo();

        // The krypo instance is modified in-place with this call.
        GATKRegistrator.registerFuncotationMapDependencies(kryo);

        // Register this class to be serialized.
        kryo.register(FuncotationMap.class);

        return kryo.copy(funcotationMap);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final FuncotationMap that = (FuncotationMap) o;

        return txToFuncotations.equals(that.txToFuncotations);
    }

    @Override
    public int hashCode() {
        return txToFuncotations.hashCode();
    }
}
