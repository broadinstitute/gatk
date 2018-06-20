package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.samples.Sex;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Stream;

@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "Functional Annotation Filtration",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class FilterFuncotations extends VariantWalker {
    private static Logger logger = LogManager.getLogger(FilterFuncotations.class);
    //==================================================================================================================
    // Arguments:

    //-----------------------------------------------------
    // Required args:

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which annotated variants should be written.")
    protected File outputFile;

    @Argument(
            shortName = "G",
            fullName = "gender",
            doc = "Sample gender for this vcf.")
    protected Sex gender;

    private VariantContextWriter outputVcfWriter;
    private String[] funcotationKeys;
    private List<FuncotationFilter> funcotationFilters = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        registerFilters();
        VCFHeader vcfHeader = getHeaderForVariants();
        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        if (funcotationHeaderLine != null) {
            funcotationKeys = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());
            outputVcfWriter = createVCFWriter(outputFile);
            funcotationFilters.forEach(filter -> vcfHeader.addMetaDataLine(new VCFFilterHeaderLine(filter.getFilterName())));
            outputVcfWriter.writeHeader(vcfHeader);
        } else {
            logger.error("Input VCF does not have Funcotator annotations.");
        }
    }

    private void registerFilters() {
        funcotationFilters.add(new ClinVarFilter(gender));
        funcotationFilters.add(new LofFilter());
        funcotationFilters.add(new LmmFilter());
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        outputVcfWriter.add(applyFilters(variant));
    }

    private VariantContext applyFilters(VariantContext variant) {
        Map<Allele, FuncotationMap> funcs = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                funcotationKeys, variant, "FilterFuncs", "FilterFuncs"
        );
        VariantContextBuilder variantContextBuilder = new VariantContextBuilder(variant);
        funcs.values().forEach(funcotationMap ->
                    funcotationFilters.forEach(filter -> {
                        if (filter.checkFilter(variant, funcotationMap)) {
                            variantContextBuilder.filter(filter.getFilterName());
                        }
                    }));
        return variantContextBuilder.make();
    }


    @Override
    public void closeTool() {
        if (outputVcfWriter != null) {
            outputVcfWriter.close();
        }
    }
}

abstract class FuncotationFiltrationRule {
    private static Logger logger = LogManager.getLogger(ClinVarFilter.class);
    private final String[] fieldNames;
    private final String ruleName;

    // Boolean flag denotes that this rule is optionally applied
    // for example:
    //  `If participant is female flag a het variant in the GLA gene`
    // In this case we don't want to apply a logical and as the test will flag all males as false. We just want to
    // optionally apply the filter to any variant matching this rule
    private final boolean optional;

    FuncotationFiltrationRule(String ruleName, boolean optional, String... fieldNames) {
        this.optional = optional;
        this.fieldNames = fieldNames;
        this.ruleName = ruleName;
    }

    abstract boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap);

    boolean applyRule(VariantContext variant, FuncotationMap funcotationMap) {
        Map<String, Stream<String>> fieldValuesMap = new HashMap<>();

        // get all funcotations for all transcripts
        Stream<Funcotation> allTranscriptValues
                = funcotationMap.getTranscriptList().stream().map(funcotationMap::get).flatMap(Collection::stream);

        // for each field name the filter cares about collect a list of values from the funcotationMap
        Arrays.stream(fieldNames).forEach(fieldName ->
                fieldValuesMap.put(fieldName, allTranscriptValues.map(funcotation ->
                        funcotation.getField(fieldName)).filter(fieldValue -> Objects.nonNull(fieldValue) && !fieldValue.isEmpty())));

        if (fieldValuesMap.isEmpty()) {
            return false;
        } else {
            return optionallyLog(ruleFunction(variant, fieldValuesMap), variant);
        }
    }

    boolean optionallyLog(Boolean result, VariantContext variant) {
        if (result) logger.warn(String.format("Matched Rule: %s For Variant %s", ruleName, variant));
        return result;
    }

    public boolean isOptional() {
        return optional;
    }

    public boolean notOptional() {
        return !optional;
    }
}

//A
class ClinVarFilter extends FuncotationFilter {
    private static final String CLIN_VAR_VCF_GENEINFO = "ClinVar_VCF_GENEINFO";
    private static final String CLIN_VAR_VCF_CLNSIG = "ClinVar_VCF_CLNSIG";
    private final Sex gender;

    ClinVarFilter(Sex gender) {
        super("CLINVAR");
        this.gender = gender;
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        final List<FuncotationFiltrationRule> clinVarFiltrationRules = new ArrayList<>();

        // 1) The gene name must be on the ACMG59 list (American College of Medical Genomics).
        //TODO

        // 2) ClinVar annotations specifies Pathogenicity or Likely pathogenic.
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-pathogenic", false, CLIN_VAR_VCF_CLNSIG) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap) {
                return fieldValueMap.get(CLIN_VAR_VCF_CLNSIG)
                        .map(clinicalSignificance ->
                                (clinicalSignificance.contains("Pathogenic") || clinicalSignificance.contains("Likely_pathogenic")))
                        .reduce(Boolean::logicalOr)
                        .orElse(false);
            }
        });

        // 3) Frequency: Max Minor Allele Freq is ≤5% in GnoMAD (ExAC for Proof of Concept)
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-MAF", false, CLIN_VAR_VCF_AF_EXAC) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap) {
                return fieldValueMap.get(CLIN_VAR_VCF_AF_EXAC)
                        .map(exacMaf -> Double.valueOf(exacMaf) <= 0.05)
                        .reduce(Boolean::logicalAnd)
                        .orElse(false);
            }
        });

        // 4) If participant is female flag a het variant in the GLA gene (x-linked) [edge case that needs more detail]
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-option-female-het-GAL", true, CLIN_VAR_VCF_GENEINFO) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap) {
                return fieldValueMap.get(CLIN_VAR_VCF_GENEINFO)
                        .map(geneInfo -> (gender == Sex.FEMALE
                                && geneInfo.substring(0, geneInfo.indexOf(":")).equals("GLA")
                                && variant.getHetCount() > 0))
                        .reduce(Boolean::logicalOr)
                        .orElse(false);
            }
        });
        return clinVarFiltrationRules;
    }
}

class LofFilter extends FuncotationFilter {

    LofFilter() {
        super("LOF");
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        // 1) 1) Variant classification is FRAME_SHIFT_*, NONSENSE, START_CODON_DEL, and SPLICE_SITE
        // (within 2 bases on either side of exon or intron) on any transcript.
        // TODO

        // 2) LoF is disease mechanism (that is do not flag genes where LoF is not part of disease mechanism e.g. RyR1)
        // - create static list
        // TODO

        // 3) Frequency: Max Minor Allele Freq is ≤1% in GnoMAD (ExAC for Proof of Concept)

        final List<FuncotationFiltrationRule> lofFiltrationRules = new ArrayList<>();

        lofFiltrationRules.add(new FuncotationFiltrationRule("LOF-MAF", false, CLIN_VAR_VCF_AF_EXAC) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap) {
                return fieldValueMap.get(CLIN_VAR_VCF_AF_EXAC)
                        .map(exacMaf -> Double.valueOf(exacMaf) <= 0.01)
                        .reduce(Boolean::logicalAnd)
                        .orElse(false);
            }
        });
        return lofFiltrationRules;
    }
}

class LmmFilter extends FuncotationFilter {

    LmmFilter() {
        super("LMM");
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        List<FuncotationFiltrationRule> lmmFiltrationRules = new ArrayList<>();

        // 1) LMM gives us a list of all path/LP variants they have seen. We flag any variant that appears on this
        // list regardless of GnoMAD freq. (optional for Proof of Concept)

        return lmmFiltrationRules;
    }
}

abstract class FuncotationFilter {
    static final String CLIN_VAR_VCF_AF_EXAC = "ClinVar_VCF_AF_EXAC";

    private final String filterName;

    FuncotationFilter(String filterName) {
        this.filterName = filterName;
    }

    Boolean checkFilter(VariantContext variant, FuncotationMap funcotationMap) {
        // apply non-optional rules first with a logical and
        Boolean nonOptionalResult = getRules().stream().filter(FuncotationFiltrationRule::notOptional).map(rule ->
                rule.applyRule(variant, funcotationMap)).reduce(Boolean::logicalAnd).orElse(false);

        // apply optional rules with a logical or
        return Boolean.logicalOr(nonOptionalResult, getRules().stream().filter(FuncotationFiltrationRule::isOptional).map(rule ->
                rule.applyRule(variant, funcotationMap)).reduce(Boolean::logicalOr).orElse(false));
    }

    abstract List<FuncotationFiltrationRule> getRules();

    public String getFilterName() {
        return filterName;
    }
}