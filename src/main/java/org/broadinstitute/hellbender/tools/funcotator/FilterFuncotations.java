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
import java.io.IOException;
import java.nio.file.Files;
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

    enum ReferenceVersion {
        hg19(19), hg38(27);

        final int gencodeVersion;

        ReferenceVersion(int gencodeVersion) {
            this.gencodeVersion = gencodeVersion;
        }
    }

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
            fullName =  FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME,
            doc = "The version of the Human Genome reference which was used to Funcotate the input VCF."
    )
    private ReferenceVersion referenceVersion;

    @Argument(
            shortName = "G",
            fullName = "gender",
            doc = "Sample gender for this vcf.")
    protected Sex gender;

    @Argument(
            shortName = "A",
            fullName = "ACMG59",
            doc = "American College of Medical Genomics list.  https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/"
    )
    protected File acmg59ListFile;

    @Argument(
            shortName = "F",
            fullName = "loflist",
            doc = "List of genes with LoF is disease mechanism."
    )
    protected File lofListFile;

    private VariantContextWriter outputVcfWriter;
    private String[] funcotationKeys;
    private List<FuncotationFilter> funcotationFilters = new ArrayList<>();
    private List<String> acmg59Genes;
    private List<String> lofGenes;

    @Override
    public void onTraversalStart() {
        try {
            acmg59Genes = Files.readAllLines(acmg59ListFile.toPath());
            lofGenes = Files.readAllLines(lofListFile.toPath());
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
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void registerFilters() {
        funcotationFilters.add(new ClinVarFilter(gender, acmg59Genes));
        funcotationFilters.add(new LofFilter(referenceVersion, lofGenes));
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

    FuncotationFiltrationRule(String ruleName, String... fieldNames) {
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

        return !fieldValuesMap.isEmpty() && optionallyLog(ruleFunction(variant, fieldValuesMap), variant);
    }

    private boolean optionallyLog(Boolean result, VariantContext variant) {
        if (result) logger.warn(String.format("Matched Rule: %s For Variant %s", ruleName, variant));
        return result;
    }
}

abstract class FuncotationFilter {
    static final String CLIN_VAR_VCF_AF_EXAC = "ClinVar_VCF_AF_EXAC";
    static final String CLIN_VAR_VCF_GENEINFO = "ClinVar_VCF_GENEINFO";
    private final String filterName;

    FuncotationFilter(String filterName) {
        this.filterName = filterName;
    }

    Boolean checkFilter(VariantContext variant, FuncotationMap funcotationMap) {
        return getRules().stream()
                .map(rule -> rule.applyRule(variant, funcotationMap))
                .reduce(Boolean::logicalAnd)
                .orElse(false);
    }

    abstract List<FuncotationFiltrationRule> getRules();

    public String getFilterName() {
        return filterName;
    }
}

//A
class ClinVarFilter extends FuncotationFilter {
    private static final String CLIN_VAR_VCF_CLNSIG = "ClinVar_VCF_CLNSIG";
    private final Sex gender;
    private final List<String> acmg59Genes;

    ClinVarFilter(Sex gender, List<String> acmg59Genes) {
        super("CLINVAR");
        this.gender = gender;
        this.acmg59Genes = acmg59Genes;
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        final List<FuncotationFiltrationRule> clinVarFiltrationRules = new ArrayList<>();

        // 1) The gene name must be on the ACMG59 list (American College of Medical Genomics).
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-ACMG59", CLIN_VAR_VCF_GENEINFO) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap) {
                return fieldValueMap.get(CLIN_VAR_VCF_GENEINFO)
                        .map(geneInfo -> acmg59Genes.contains(geneInfo.substring(0, geneInfo.indexOf(":"))))
                        .reduce(Boolean::logicalOr)
                        .orElse(false);
            }
        });

        // 2) ClinVar annotations specifies Pathogenicity or Likely pathogenic.
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-pathogenic", CLIN_VAR_VCF_CLNSIG) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap) {
                return fieldValueMap.get(CLIN_VAR_VCF_CLNSIG)
                        .map(clinicalSignificance ->
                                clinicalSignificance.contains("Pathogenic") || clinicalSignificance.contains("Likely_pathogenic"))
                        .reduce(Boolean::logicalOr)
                        .orElse(false);
            }
        });

        // 3) Frequency: Max Minor Allele Freq is ≤5% in GnoMAD (ExAC for Proof of Concept)
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-MAF", CLIN_VAR_VCF_AF_EXAC) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap) {
                return fieldValueMap.get(CLIN_VAR_VCF_AF_EXAC)
                        .map(exacMaf -> Double.valueOf(exacMaf) <= 0.05)
                        .reduce(Boolean::logicalOr)
                        .orElse(false);
            }
        });

        // 4) If participant is female flag a het variant in the GLA gene (x-linked) [edge case that needs more detail]
        if (gender == Sex.FEMALE) {
            clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-option-female-het-GAL", CLIN_VAR_VCF_GENEINFO) {
                @Override
                boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap) {
                    return variant.getHetCount() > 0 && fieldValueMap.get(CLIN_VAR_VCF_GENEINFO)
                            .map(geneInfo -> geneInfo.startsWith("GLA:"))
                            .reduce(Boolean::logicalOr)
                            .orElse(false);
                }
            });
        }
        return clinVarFiltrationRules;
    }
}

class LofFilter extends FuncotationFilter {

    private static final String FRAME_SHIFT_PREFIX = "FRAME_SHIFT_";
    private static final List<String> CONSTANT_LOF_CLASSIFICATIONS = Arrays.asList("NONSENSE", "START_CODON_DEL", "SPLICE_SITE");

    private final String classificationFuncotation;
    private final List<String> lofGenes;

    LofFilter(FilterFuncotations.ReferenceVersion ref, List<String> lofGenes) {
        super("LOF");
        this.classificationFuncotation = "Gencode_" + ref.gencodeVersion + "_variantClassification";
        this.lofGenes = lofGenes;
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        final List<FuncotationFiltrationRule> lofFiltrationRules = new ArrayList<>();
        // 1) 1) Variant classification is FRAME_SHIFT_*, NONSENSE, START_CODON_DEL, and SPLICE_SITE
        // (within 2 bases on either side of exon or intron) on any transcript.
        // TODO
        lofFiltrationRules.add(new FuncotationFiltrationRule("LOF-CLASS", classificationFuncotation) {

            @Override
            boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap) {
                return fieldValueMap.get(classificationFuncotation)
                        .map(classification ->
                                classification.startsWith(FRAME_SHIFT_PREFIX) || CONSTANT_LOF_CLASSIFICATIONS.contains(classification))
                        .reduce(Boolean::logicalOr)
                        .orElse(false);
            }
        });

        // 2) LoF is disease mechanism (that is do not flag genes where LoF is not part of disease mechanism e.g. RyR1)
        // - create static list
        lofFiltrationRules.add(new FuncotationFiltrationRule("Lof-mechanism", CLIN_VAR_VCF_GENEINFO) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap) {
                return fieldValueMap.get(CLIN_VAR_VCF_GENEINFO)
                        .map(geneInfo -> lofGenes.contains(geneInfo.substring(0, geneInfo.indexOf(":"))))
                        .reduce(Boolean::logicalOr)
                        .orElse(false);
            }
        });

        // 3) Frequency: Max Minor Allele Freq is ≤1% in GnoMAD (ExAC for Proof of Concept)
        lofFiltrationRules.add(new FuncotationFiltrationRule("LOF-MAF", CLIN_VAR_VCF_AF_EXAC) {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, Stream<String>> fieldValueMap) {
                return fieldValueMap.get(CLIN_VAR_VCF_AF_EXAC)
                        .map(exacMaf -> Double.valueOf(exacMaf) <= 0.01)
                        .reduce(Boolean::logicalOr)
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
