package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@CommandLineProgramProperties(
        summary = FilterFuncotations.SUMMARY,
        oneLineSummary = FilterFuncotations.ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class FilterFuncotations extends VariantWalker {

    static final String ONE_LINE_SUMMARY = "Filter variants based on clinically-significant Funcotations.";
    static final String SUMMARY = ONE_LINE_SUMMARY +
            " Proof-of-concept hard-coded to look for specific Funcotations from ClinVar, ExAC, and LMM.";

    static final String CLINSIG_RULE_KEY = "CLINSIG";
    static final String NOT_CLINSIG_FILTER = "NOT_" + CLINSIG_RULE_KEY;

    enum ReferenceVersion {
        hg19(19), hg38(27);

        final int gencodeVersion;

        ReferenceVersion(int gencodeVersion) {
            this.gencodeVersion = gencodeVersion;
        }
    }

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

    private VariantContextWriter outputVcfWriter;
    private String[] funcotationKeys;
    private List<FuncotationFilter> funcotationFilters = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        registerFilters();
        final VCFHeader vcfHeader = getHeaderForVariants();

        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        if (funcotationHeaderLine != null) {
            funcotationKeys = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());
            outputVcfWriter = createVCFWriter(outputFile);
            vcfHeader.addMetaDataLine(new VCFFilterHeaderLine(NOT_CLINSIG_FILTER, "Filter for clinically insignificant variants."));
            vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(CLINSIG_RULE_KEY, 1, VCFHeaderLineType.String,
                    "Rule(s) which caused this annotation to be flagged as clinically significant."));
            outputVcfWriter.writeHeader(vcfHeader);
        } else {
            throw new UserException.BadInput("Input VCF does not have Funcotator annotations.");
        }
    }

    private void registerFilters() {
        funcotationFilters.add(new ClinVarFilter());
        funcotationFilters.add(new LofFilter(referenceVersion));
        funcotationFilters.add(new LmmFilter());
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        outputVcfWriter.add(applyFilters(variant));
    }

    private VariantContext applyFilters(VariantContext variant) {

        final Set<String> matchingFilters = new HashSet<>();
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(variant);

        final Map<Allele, FuncotationMap> funcs = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                funcotationKeys, variant, "Gencode_" + referenceVersion.gencodeVersion + "_annotationTranscript", "FILTER");

        funcs.values().forEach(funcotationMap -> {
            final Stream<Map<String, String>> transcriptFuncotations = funcotationMap.getTranscriptList().stream()
                    .map(funcotationMap::get)
                    .map(funcotations -> funcotations.stream()
                            .flatMap(this::extractFuncotationFields)
                            .filter(entry -> entry.getValue() != null && !entry.getValue().isEmpty())
                            .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)));

            transcriptFuncotations.forEach(funcotations -> {
                final Set<String> matches = funcotationFilters.stream()
                        .filter(f -> f.checkFilter(variant, funcotations))
                        .map(FuncotationFilter::getFilterName)
                        .collect(Collectors.toSet());
                matchingFilters.addAll(matches);
            });
        });

        String clinicalSignificance = matchingFilters.isEmpty() ? "NONE" : String.join(",", matchingFilters);
        variantContextBuilder.attribute(CLINSIG_RULE_KEY, clinicalSignificance);

        if (matchingFilters.isEmpty()) {
            variantContextBuilder.filter(NOT_CLINSIG_FILTER);
        } else {
            variantContextBuilder.passFilters();
        }
        return variantContextBuilder.make();
    }

    private Stream<Map.Entry<String, String>> extractFuncotationFields(final Funcotation funcotation) {
        return funcotation.getFieldNames().stream()
                .map(name -> new AbstractMap.SimpleEntry<>(name, funcotation.getField(name)));
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

    FuncotationFiltrationRule(final String ruleName) {
        this.ruleName = ruleName;
    }

    boolean checkRule(final VariantContext variant,final Map<String, String> prunedTranscriptFuncotations) {
        return !prunedTranscriptFuncotations.isEmpty() &&
                optionallyLog(applyRule(variant, prunedTranscriptFuncotations), variant);
    }

    abstract boolean applyRule(final VariantContext variant, final Map<String, String> prunedTranscriptFuncotations);

    private boolean optionallyLog(final boolean result, final VariantContext variant) {
        if (result) logger.debug(String.format("Matched Rule: %s For Variant %s", ruleName, variant));
        return result;
    }

    double getMaxMinorAlleleFreq(final Map<String, String> funcotations) {
        return Arrays.stream(ExacSubPopulation.values())
                .filter(subpop -> funcotations.containsKey("ExAC_AC_" + subpop.name()))
                .map(subpop -> {
                    final Double ac = Double.valueOf(funcotations.get("ExAC_AC_" + subpop.name()));
                    final Integer an = Integer.valueOf(funcotations.get("ExAC_AN_" + subpop.name()));

                    if (an == 0) {
                        // If a variant has never been seen in ExAC, report it as 0% MAF.
                        return 0d;
                    } else {
                        return ac / an;
                    }
                })
                .max(Double::compareTo)
                .orElse(0d);
    }
}

abstract class FuncotationFilter {

    private final String filterName;

    FuncotationFilter(final String filterName) {
        this.filterName = filterName;
    }

    Boolean checkFilter(final VariantContext variant, final Map<String, String> prunedTranscriptFuncotations) {
        return getRules().stream()
                .map(rule -> rule.checkRule(variant, prunedTranscriptFuncotations))
                .reduce(Boolean::logicalAnd)
                .orElse(false);
    }

    abstract List<FuncotationFiltrationRule> getRules();

    public String getFilterName() {
        return filterName;
    }
}


class ClinVarFilter extends FuncotationFilter {

    private static final String ACMG_DISEASE_FUNCOTATION = "ACMG_recommendation_Disease_Name";
    private static final String CLIN_VAR_VCF_CLNSIG = "ClinVar_VCF_CLNSIG";

    ClinVarFilter() {
        super("CLINVAR");
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        final List<FuncotationFiltrationRule> clinVarFiltrationRules = new ArrayList<>();

        // 1) The gene name must be on the ACMG59 list (American College of Medical Genomics).
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-ACMG59") {
            @Override
            boolean applyRule(final VariantContext variant, final Map<String, String> prunedTranscriptFuncotations) {
                return prunedTranscriptFuncotations.containsKey(ACMG_DISEASE_FUNCOTATION);
            }
        });

        // 2) ClinVar annotations specifies Pathogenicity or Likely pathogenic.
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-pathogenic") {
            @Override
            boolean applyRule(final VariantContext variant, final Map<String, String> prunedTranscriptFuncotations) {
                final String significance = prunedTranscriptFuncotations.getOrDefault(CLIN_VAR_VCF_CLNSIG, "");
                return significance.contains("Pathogenic") || significance.contains("Likely_pathogenic");
            }
        });

        // 3) Frequency: Max Minor Allele Freq is ≤5% in GnoMAD (ExAC for Proof of Concept)
        clinVarFiltrationRules.add(new FuncotationFiltrationRule("ClinVar-MAF") {
            @Override
            boolean applyRule(final VariantContext variant, final Map<String, String> prunedTranscriptFuncotations) {
                return getMaxMinorAlleleFreq(prunedTranscriptFuncotations) <= 0.05;
            }
        });
        return clinVarFiltrationRules;
    }
}

class LofFilter extends FuncotationFilter {

    private static final String LOF_GENE_FUNCOTATION = "ACMGLMMLof_LOF_Mechanism";
    private static final String FRAME_SHIFT_PREFIX = "FRAME_SHIFT_";
    private static final List<String> CONSTANT_LOF_CLASSIFICATIONS = Arrays.asList("NONSENSE", "START_CODON_DEL", "SPLICE_SITE");

    private final String classificationFuncotation;

    LofFilter(final FilterFuncotations.ReferenceVersion ref) {
        super("LOF");
        this.classificationFuncotation = "Gencode_" + ref.gencodeVersion + "_variantClassification";
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        final List<FuncotationFiltrationRule> lofFiltrationRules = new ArrayList<>();
        // 1) 1) Variant classification is FRAME_SHIFT_*, NONSENSE, START_CODON_DEL, and SPLICE_SITE
        // (within 2 bases on either side of exon or intron) on any transcript.
        // TODO
        lofFiltrationRules.add(new FuncotationFiltrationRule("LOF-class") {

            @Override
            boolean applyRule(final VariantContext variant, final Map<String, String> prunedTranscriptFuncotations) {
                final String classification = prunedTranscriptFuncotations.getOrDefault(classificationFuncotation, "");
                return classification.startsWith(FRAME_SHIFT_PREFIX) || CONSTANT_LOF_CLASSIFICATIONS.contains(classification);
            }
        });

        // 2) LoF is disease mechanism (that is do not flag genes where LoF is not part of disease mechanism e.g. RyR1)
        // - create static list
        lofFiltrationRules.add(new FuncotationFiltrationRule("LOF-mechanism") {
            @Override
            boolean applyRule(final VariantContext variant, final Map<String, String> prunedTranscriptFuncotations) {
                return prunedTranscriptFuncotations.getOrDefault(LOF_GENE_FUNCOTATION, "NO").equals("YES");
            }
        });

        // 3) Frequency: Max Minor Allele Freq is ≤1% in GnoMAD (ExAC for Proof of Concept)
        lofFiltrationRules.add(new FuncotationFiltrationRule("LOF-MAF") {
            @Override
            boolean applyRule(final VariantContext variant, final Map<String, String> prunedTranscriptFuncotations) {
                return getMaxMinorAlleleFreq(prunedTranscriptFuncotations) <= 0.01;
            }
        });
        return lofFiltrationRules;
    }
}

class LmmFilter extends FuncotationFilter {

    private static final String LMM_FLAGGED = "LMMKnown_LMM_FLAGGED";

    LmmFilter() {
        super("LMM");
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        // 1) LMM gives us a list of all path/LP variants they have seen. We flag any variant that appears on this
        // list regardless of GnoMAD freq. (optional for Proof of Concept)
        final FuncotationFiltrationRule rule = new FuncotationFiltrationRule("LMM-path-LP") {
            @Override
            boolean applyRule(final VariantContext variant, final Map<String, String> prunedTranscriptFuncotations) {
                return Boolean.valueOf(prunedTranscriptFuncotations.getOrDefault(LMM_FLAGGED, "false"));
            }
        };
        return Collections.singletonList(rule);
    }
}
