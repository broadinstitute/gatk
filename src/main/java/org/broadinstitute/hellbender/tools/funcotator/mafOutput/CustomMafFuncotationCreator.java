package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.TumorNormalPair;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;


/**
 * Produces custom MAF fields (e.g. t_alt_count) into a Funcotation that can be used by the {@link org.broadinstitute.hellbender.tools.funcotator.mafOutput.MafOutputRenderer}
 */
public class CustomMafFuncotationCreator {

    public final static List<String> COUNT_FIELD_NAMES = Arrays.asList(
            MafOutputRendererConstants.FieldName_t_alt_count,
            MafOutputRendererConstants.FieldName_t_ref_count,
            MafOutputRendererConstants.FieldName_n_alt_count,
            MafOutputRendererConstants.FieldName_n_ref_count,
            MafOutputRendererConstants.FieldName_tumor_f);

    public final static String MAF_DBSNP_VAL_STATUS_FIELD = "custom_dbsnp_val_status";
    public final static String BY_FREQ = "byFrequency";
    public final static String BY_1KG = "by1000genomes";
    public final static Map<String, String> DBSNP_VALIDATION_VALUE_MAP =
            ImmutableMap.<String,String>builder()
                    .put(MafOutputRendererConstants.DBSNP_DS_NAME + "_VLD", BY_FREQ)
                    .put(MafOutputRendererConstants.DBSNP_DS_NAME + "_KGVAL", BY_1KG)
                    .put(MafOutputRendererConstants.DBSNP_DS_NAME + "_KGPROD", BY_1KG)
                    .put(MafOutputRendererConstants.DBSNP_DS_NAME + "_KGPilot1", BY_1KG)
                    .put(MafOutputRendererConstants.DBSNP_DS_NAME + "_KGPilot123", BY_1KG)
                    .put(MafOutputRendererConstants.DBSNP_DS_NAME + "_KGPilot3", BY_1KG)
                    .put(MafOutputRendererConstants.DBSNP_DS_NAME + "_KGPhase1", BY_1KG)
                    .put(MafOutputRendererConstants.DBSNP_DS_NAME + "_KGPhase3", BY_1KG)
                    .build();
    public static final String MAF_DBSNP_VAL_STATUS_DELIMITER = ";";

    /**
     * Create the MAF count fields (See {@link CustomMafFuncotationCreator#COUNT_FIELD_NAMES}) for a single pair in the given variant.
     *
     * If no samples can be linked to a tumor and normal combination, then blank values are generated for the count fields.
     * @param variant Never {@code null}
     * @param tnPair Never {@code null}
     * @return List of funcotations for the pair.  This list will usually be of length one, unless multiallelics are involved.
     *  Never {@code null}  Field values can be blank:
     *  - if there are no AD or AF format fields
     *  - if the specified pair is not in the variant sample list.
     */
    private static List<Funcotation> createCustomMafCountFields(final VariantContext variant, final TumorNormalPair tnPair) {
        Utils.nonNull(variant);
        Utils.nonNull(tnPair);
        final List<Funcotation> result = new ArrayList<>();

        for (int i = 0; i < variant.getAlternateAlleles().size(); i++) {
            final Allele allele = variant.getAlternateAllele(i);
            final String tumorSampleName = tnPair.getTumor();
            final String normalSampleName = tnPair.getNormal();

            // This code assumes that the AD field is of Type "R"
            final Genotype tumorGenotype = variant.getGenotype(tumorSampleName);
            final boolean hasTumorAD = (tumorGenotype != null) && (tumorGenotype.hasAD());
            final boolean hasTumorAF = (tumorGenotype != null) && (tumorGenotype.hasAnyAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY));
            final String tAltCount = hasTumorAD ? Integer.toString(tumorGenotype.getAD()[i + 1]) : "";
            final String tRefCount = hasTumorAD ? Integer.toString(tumorGenotype.getAD()[0]) : "";
            final String tumorFAllAlleles = hasTumorAF ? tumorGenotype.getAnyAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY).toString() : "";
            final String tumorF = hasTumorAF ? StringUtils.split(tumorFAllAlleles, ",")[i] : "";

            final Genotype normalGenotype = variant.getGenotype(normalSampleName);
            final boolean hasNormalAD = (normalGenotype != null) && (normalGenotype.hasAD());
            final String nAltCount = hasNormalAD ? Integer.toString(normalGenotype.getAD()[i + 1]) : "";
            final String nRefCount = hasNormalAD ? Integer.toString(normalGenotype.getAD()[0]) : "";

            final List<String> fieldValues = Arrays.asList(
                    tAltCount,
                    tRefCount,
                    nAltCount,
                    nRefCount,
                    tumorF);

            result.add(TableFuncotation.create(COUNT_FIELD_NAMES, fieldValues, allele, MafOutputRendererConstants.MAF_COUNT_RENDERING_DATASOURCE_DUMMY_NAME, createCustomMafCountFieldsMetadata()));
        }

        return result;
    }

    private static FuncotationMetadata createCustomMafCountFieldsMetadata() {
        return VcfFuncotationMetadata.create(Arrays.asList(
                new VCFInfoHeaderLine(COUNT_FIELD_NAMES.get(0), VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Number of alternate reads in the tumor."),
                new VCFInfoHeaderLine(COUNT_FIELD_NAMES.get(1), VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Number of reference reads in the tumor."),
                new VCFInfoHeaderLine(COUNT_FIELD_NAMES.get(2), VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Number of alternate reads in the normal."),
                new VCFInfoHeaderLine(COUNT_FIELD_NAMES.get(3), VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Number of reference reads in the normal."),
                new VCFInfoHeaderLine(COUNT_FIELD_NAMES.get(4), VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele fractions of alternate alleles in the tumor.")
        ));
    }

    private static FuncotationMetadata createCustomMafDbSnpFieldsMetadata() {
        return VcfFuncotationMetadata.create(Collections.singletonList(
                new VCFInfoHeaderLine(MAF_DBSNP_VAL_STATUS_FIELD, VCFHeaderLineCount.A, VCFHeaderLineType.String, "dbSNP Validation status.")
        ));
    }

    /**
     * Create the fields in {@link CustomMafFuncotationCreator#COUNT_FIELD_NAMES} for the pair given in the variant.
     *
     * @param variant  Never {@code null}
     * @param tnPairs Should be length 0 or 1.  Others beyond the first are ignored.  See https://github.com/broadinstitute/gatk/issues/4912  Never {@code null}.
     * @return Funcotations for the fields in {@link CustomMafFuncotationCreator#COUNT_FIELD_NAMES}.  If there are no pairs as the input, return an empty list.
     * Never {@code null}
     */
    public static List<Funcotation> createCustomMafCountFields(final VariantContext variant, final List<TumorNormalPair> tnPairs) {
        Utils.nonNull(variant);
        Utils.nonNull(tnPairs);
        if (tnPairs.size() == 0) {
            return createEmptyCountFields(variant.getAlternateAlleles());
        }

        // Note that we assume that there is only one pair here.  When more pairs are possible, return a flattened list.
        // TODO: We want to support more than one pair eventually.  See https://github.com/broadinstitute/gatk/issues/4912
        return CustomMafFuncotationCreator.createCustomMafCountFields(variant, tnPairs.get(0));
    }

    /**
     * Produce one empty Funcotation per alternate allele.  Each will have the fields that should be produced by this class,
     *  but the values will be blank.
     * @param alleles alternate alleles to create the funcotations for.  Never {@code null}.
     * @return Never {@code null}.  Empty list if input list is empty.
     */
    private static List<Funcotation> createEmptyCountFields(final List<Allele> alleles) {
        Utils.nonNull(alleles);
        final List<Funcotation> result = new ArrayList<>();

        for (final Allele allele : alleles) {
            result.add(TableFuncotation.create(COUNT_FIELD_NAMES, Collections.nCopies(COUNT_FIELD_NAMES.size(), ""), allele, MafOutputRendererConstants.MAF_COUNT_RENDERING_DATASOURCE_DUMMY_NAME, createCustomMafCountFieldsMetadata()));
        }
        return result;
    }

    /**
     * Creates a list of funcotations for each input dbSNP funcotation.  Input funcotations that are not dbSNP are ignored.
     *
     * Note that the determination is by using a default dbSNP name {@link MafOutputRendererConstants#DBSNP_DS_NAME}
     *
     * @param funcotations Never {@code null}
     * @return list of the custom dbSNP values.  If input is empty,  {@code null}, or input contains no discernible
     *  dbSNP funcotations, return an empty list.  The alt allele will be the same as the corresponding input funcotation.
     */
    public static List<Funcotation> createCustomMafDbSnpFields(final List<Funcotation> funcotations) {
        Utils.nonNull(funcotations);
        if (funcotations.size() == 0) {
            return Collections.emptyList();
        }

        // Only process funcotations that came from a db snp funcotation factory.
        final List<Funcotation> sourceDbSnpFuncotationsToGenerate = funcotations.stream()
                .filter(f -> f.getDataSourceName().startsWith(MafOutputRendererConstants.DBSNP_DS_NAME))
                .collect(Collectors.toList());

        return sourceDbSnpFuncotationsToGenerate.stream()
                .map(CustomMafFuncotationCreator::createDbSnpCustomFuncotation)
                .collect(Collectors.toList());
    }

    private static Funcotation createDbSnpCustomFuncotation(final Funcotation dbSnpFuncotation) {
        return TableFuncotation.create(
                Collections.singletonList(MAF_DBSNP_VAL_STATUS_FIELD),
                Collections.singletonList(DBSNP_VALIDATION_VALUE_MAP.entrySet().stream()
                    .filter(e -> dbSnpFuncotation.getFieldOrDefault(e.getKey(), "").toLowerCase().contains("true"))
                    .map(e -> DBSNP_VALIDATION_VALUE_MAP.getOrDefault(e.getKey(), ""))
                    .distinct()
                    .collect(Collectors.joining(MAF_DBSNP_VAL_STATUS_DELIMITER))),
                dbSnpFuncotation.getAltAllele(),
                MafOutputRendererConstants.MAF_DBSNP_RENDERING_DATASOURCE_DUMMY_NAME,
                createCustomMafDbSnpFieldsMetadata());
    }
}
