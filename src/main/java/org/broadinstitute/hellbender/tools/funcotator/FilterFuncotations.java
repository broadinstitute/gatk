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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
                funcotationMap.getTranscriptList().forEach(transcriptId -> {
                    funcotationFilters.forEach(filter -> {
                        if (filter.checkFilter(variant, funcotationMap.get(transcriptId))) {
                            variantContextBuilder.filter(filter.getFilterName());
                        }
                    });
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
    private final String[] fieldNames;

    FuncotationFiltrationRule(String... fieldNames) {
        this.fieldNames = fieldNames;
    }

    abstract boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap);

    boolean applyRule(VariantContext variant, Funcotation funcotation) {
        Map<String, String> fieldValueMap = new HashMap<>();

        Arrays.stream(fieldNames).forEach(fieldName -> fieldValueMap.put(fieldName, funcotation.getField(fieldName)));
        fieldValueMap.values().removeIf(value -> (value == null || value.isEmpty()));

        if (fieldValueMap.isEmpty()) {
            return false;
        } else {
            return ruleFunction(variant, fieldValueMap);
        }
    }
}

//A
class ClinVarFilter extends FuncotationFilter {
    private static final String CLIN_VAR_VCF_GENEINFO = "ClinVar_VCF_GENEINFO";
    private static final String CLIN_VAR_VCF_CLNSIG = "ClinVar_VCF_CLNSIG";
    private static final String CLIN_VAR_VCF_AF_EXAC = "ClinVar_VCF_AF_EXAC";
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
        clinVarFiltrationRules.add(new FuncotationFiltrationRule(CLIN_VAR_VCF_CLNSIG) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap) {
                String clinicalSignificance = fieldValueMap.get(CLIN_VAR_VCF_CLNSIG);
                return clinicalSignificance.contains("Pathogenic") || clinicalSignificance.contains("Likely_pathogenic");
            }
        });

        // 3) Frequency: Max Minor Allele Freq is ≤5% in GnoMAD (ExAC for Proof of Concept)
        clinVarFiltrationRules.add(new FuncotationFiltrationRule(CLIN_VAR_VCF_AF_EXAC) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap) {
                Double alleleFreqExac = Double.valueOf(fieldValueMap.get(CLIN_VAR_VCF_AF_EXAC));
                return alleleFreqExac <= 0.05;
            }
        });

        // 4) If participant is female flag a het variant in the GLA gene (x-linked) [edge case that needs more detail]
        clinVarFiltrationRules.add(new FuncotationFiltrationRule(CLIN_VAR_VCF_GENEINFO, CLIN_VAR_VCF_CLNSIG) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap) {
                String geneInfo = fieldValueMap.get(CLIN_VAR_VCF_GENEINFO);
                return Sex.FEMALE == gender && geneInfo.substring(0, geneInfo.indexOf(":")).equals("GLA")
                        && variant.getHetCount() > 0;
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
        return new ArrayList<>();
    }
}

class LmmFilter extends FuncotationFilter {

    LmmFilter() {
        super("LMM");
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        return new ArrayList<>();
    }
}

abstract class FuncotationFilter {

    private final String filterName;

    FuncotationFilter(String filterName) {
        this.filterName = filterName;
    }

    Boolean checkFilter(VariantContext variant, List<Funcotation> funcotations) {
        //check each funcotation against each rule and reduce
        return funcotations.stream().map(funcotation ->
                getRules().stream().map(rule ->
                        rule.applyRule(variant, funcotation)).reduce(false, (a, b) -> a || b))
                .reduce(false, (a, b) -> a || b);
    }

    abstract List<FuncotationFiltrationRule> getRules();

    public String getFilterName() {
        return filterName;
    }
}