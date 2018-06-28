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
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
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
    protected ReferenceVersion referenceVersion;

    @Argument(
            fullName = "acmg59-list",
            doc = "American College of Medical Genomics list. https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/"
    )
    protected File acmg59ListFile;

    @Argument(
            fullName = "lof-list",
            doc = "List of genes with LoF is disease mechanism."
    )
    protected File lofListFile;

    private VariantContextWriter outputVcfWriter;
    private String[] funcotationKeys;
    private List<FuncotationFilter> funcotationFilters = new ArrayList<>();
    private List<String> acmg59Genes;
    private List<String> lofGenes;
    // FIXME: Read from VCF header.
    private Sex gender = Sex.UNKNOWN;

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
        funcotationFilters.add(new ClinVarFilter(acmg59Genes));
        if (gender == Sex.FEMALE) {
            funcotationFilters.add(new ClinVarFemaleFilter());
        }
        funcotationFilters.add(new LofFilter(referenceVersion, lofGenes));
        funcotationFilters.add(new LmmFilter());
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        outputVcfWriter.add(applyFilters(variant));
    }

    private VariantContext applyFilters(VariantContext variant) {
        Map<Allele, FuncotationMap> funcs = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                funcotationKeys, variant, "Gencode_" + referenceVersion.gencodeVersion + "_annotationTranscript", "FILTER"
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

    enum ExacSubPopulation {
        AFR, AMR, EAS, FIN, NFE, OTH, SAS
    }

    private static Logger logger = LogManager.getLogger(ClinVarFilter.class);
    private final String ruleName;

    FuncotationFiltrationRule(String ruleName) {
        this.ruleName = ruleName;
    }

    abstract boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap);

    boolean applyRule(VariantContext variant, FuncotationMap funcotationMap) {

        final Stream<Map<String, String>> funcotationsByTranscript = funcotationMap.getTranscriptList().stream()
                .map(funcotationMap::get).map(funcotations ->
                        funcotations.stream()
                                .flatMap(this::extractFuncotationFields)
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)));

        return funcotationsByTranscript.anyMatch(funcotationValues ->
                !funcotationValues.isEmpty() && optionallyLog(ruleFunction(variant, funcotationValues), variant));
    }

    private Stream<Map.Entry<String, String>> extractFuncotationFields(final Funcotation funcotation) {
        return funcotation.getFieldNames().stream()
                .map(name -> new AbstractMap.SimpleEntry<>(name, funcotation.getField(name)));
    }

    private boolean optionallyLog(Boolean result, VariantContext variant) {
        if (result) logger.warn(String.format("Matched Rule: %s For Variant %s", ruleName, variant));
        return result;
    }

    Stream<Double> getMaxMinorAlleleFreqs(int alleleCount, final Map<String, String> funcotations) {
        final double[] maxMafsByAllele = new double[alleleCount];
        for (int i = 0; i < alleleCount; i++) {
            maxMafsByAllele[i] = Double.NaN;
        }

        Arrays.stream(ExacSubPopulation.values()).forEach(subpop -> {
            final Optional<String> alleleCountsString = Optional.ofNullable(funcotations.get("ExAC_AC_" + subpop.name()));
            alleleCountsString.ifPresent(countsString -> {
                final String[] alleleCounts = countsString.split("_[^_]+_");
                final int chromCount = Integer.valueOf(funcotations.getOrDefault("ExAC_AN_" + subpop.name(), "0"));

                for (int i = 0; i < alleleCount; i++) {
                    final double maf = Double.valueOf(alleleCounts[i]) / chromCount;
                    if (maxMafsByAllele[i] < maf) {
                        maxMafsByAllele[i] = maf;
                    }
                }
            });
        });

        return Arrays.stream(maxMafsByAllele).boxed();
    }
}

abstract class FuncotationFilter {
    static final String CLIN_VAR_VCF_CLNSIG = "ClinVar_VCF_CLNSIG";
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


class ClinVarFilter extends FuncotationFilter {
    private final List<String> acmg59Genes;

    ClinVarFilter(List<String> acmg59Genes) {
        super("CLINVAR");
        this.acmg59Genes = acmg59Genes;
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        final List<FuncotationFiltrationRule> clinVarFiltrationRules = new ArrayList<>();

        // 1) The gene name must be on the ACMG59 list (American College of Medical Genomics).
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-ACMG59") {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap) {
                final String geneInfo = fieldValueMap.getOrDefault(CLIN_VAR_VCF_GENEINFO, "");
                return acmg59Genes.contains(geneInfo.substring(0, geneInfo.indexOf(":")));
            }
        });

        // 2) ClinVar annotations specifies Pathogenicity or Likely pathogenic.
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-pathogenic") {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap) {
                final String significance = fieldValueMap.getOrDefault(CLIN_VAR_VCF_CLNSIG, "");
                return significance.contains("Pathogenic") || significance.contains("Likely_pathogenic");
            }
        });

        // 3) Frequency: Max Minor Allele Freq is ≤5% in GnoMAD (ExAC for Proof of Concept)
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-MAF") {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap) {
                return getMaxMinorAlleleFreqs(variant.getAlternateAlleles().size(), fieldValueMap)
                        .anyMatch(d -> d <= 0.05);
            }
        });
        return clinVarFiltrationRules;
    }
}

class ClinVarFemaleFilter extends FuncotationFilter {

    ClinVarFemaleFilter() {
        super("CLINVAR");
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        // 4) If participant is female flag a het variant in the GLA gene (x-linked) [edge case that needs more detail]
        final FuncotationFiltrationRule rule = new FuncotationFiltrationRule("ClinVar-option-female-het-GAL") {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap) {
                return variant.getHetCount() > 0 &&
                        fieldValueMap.getOrDefault(CLIN_VAR_VCF_GENEINFO, "").startsWith("GLA:");
            }
        };
        return Collections.singletonList(rule);
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
        lofFiltrationRules.add(new FuncotationFiltrationRule("LOF-class") {

            @Override
            boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap) {
                final String classification = fieldValueMap.getOrDefault(classificationFuncotation, "");
                return classification.startsWith(FRAME_SHIFT_PREFIX) || CONSTANT_LOF_CLASSIFICATIONS.contains(classification);
            }
        });

        // 2) LoF is disease mechanism (that is do not flag genes where LoF is not part of disease mechanism e.g. RyR1)
        // - create static list
        lofFiltrationRules.add(new FuncotationFiltrationRule("LOF-mechanism") {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap) {
                final String geneInfo = fieldValueMap.getOrDefault(CLIN_VAR_VCF_GENEINFO, "");
                return lofGenes.contains(geneInfo.substring(0, geneInfo.indexOf(":")));
            }
        });

        // 3) Frequency: Max Minor Allele Freq is ≤1% in GnoMAD (ExAC for Proof of Concept)
        lofFiltrationRules.add(new FuncotationFiltrationRule("LOF-MAF") {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap) {
                return getMaxMinorAlleleFreqs(variant.getAlternateAlleles().size(), fieldValueMap)
                        .anyMatch(d -> d <= 0.01);
            }
        });
        return lofFiltrationRules;
    }
}

class LmmFilter extends FuncotationFilter {

    private static final String LMM_FLAGGED = "LMMKnown_LMM_Flagged";

    LmmFilter() {
        super("LMM");
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        // 1) LMM gives us a list of all path/LP variants they have seen. We flag any variant that appears on this
        // list regardless of GnoMAD freq. (optional for Proof of Concept)
        final FuncotationFiltrationRule rule = new FuncotationFiltrationRule("LMM-path-LP") {
            @Override
            boolean ruleFunction(VariantContext variant, Map<String, String> fieldValueMap) {
                return Boolean.valueOf(fieldValueMap.getOrDefault(LMM_FLAGGED, "false"));
            }
        };
        return Collections.singletonList(rule);
    }
}
