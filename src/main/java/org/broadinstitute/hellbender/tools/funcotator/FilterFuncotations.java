package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.TwoPassVariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.ArHetvarFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.ArHomvarFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.ClinVarFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.FuncotationFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.LmmFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.LofFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.TwoPassFuncotationFilter;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Filter variants based on clinically-significant Funcotations.
 *
 * This proof-of-concept tool is an example for how to parse and use the VCF output of Funcotator.
 * It's currently hard-coded to look for specific {@link Funcotation}s from:
 * <ul>
 *     <li><a href="http://www.clinvar.com/">ClinVar</a></li>
 *     <li><a href="http://personalizedmedicine.partners.org/laboratory-for-molecular-medicine/">Laboratory for Molecular Medicine (LMM)</a></li>
 * </ul>
 * It also looks for Funcotations from whichever of the following data sets is specified by the user:
 * <ul>
 *     <li><a href="http://exac.broadinstitute.org/">Exome Aggregation Consortium (ExAC)</a></li>
 *     <li><a href="http://gnomad.broadinstitute.org/">Genome Aggregation Database (gnomAD)</a></li>
 * </ul>
 */
@CommandLineProgramProperties(
        summary = FilterFuncotations.SUMMARY,
        oneLineSummary = FilterFuncotations.ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class FilterFuncotations extends TwoPassVariantWalker {

    static final String ONE_LINE_SUMMARY = "Filter variants based on clinically-significant Funcotations.";
    static final String SUMMARY = ONE_LINE_SUMMARY +
            "\nThis proof-of-concept tool is an example for how to parse and use the VCF output of Funcotator." +
            "\nCurrently hard-coded to look for specific Funcotations from:" +
            "\n  * ClinVar (http://www.clinvar.com/)" +
            "\n  * Laboratory for Molecular Medicine (LMM) (http://personalizedmedicine.partners.org/laboratory-for-molecular-medicine/)" +
            "\nAlso looks for Funcotations from whichever of the following is specified by the user:" +
            "\n  * Exome Aggregation Consortium (ExAC) (http://exac.broadinstitute.org/)" +
            "\n  * Genome Aggregation Database (gnomAD) (http://gnomad.broadinstitute.org/)";

    /**
     * The version of the Human Genome reference which was used when Funcotating the input VCF.
     *
     * Used to derive names of Gencode Funcotations.
     */
    public enum Reference {
        b37(19), hg19(19), hg38(27);

        private final int gencodeVersion;

        Reference(int gencodeVersion) {
            this.gencodeVersion = gencodeVersion;
        }

        public int getGencodeVersion() {
            return gencodeVersion;
        }
    }

    /**
     * The allele frequency data source that was used when Funcotating the input VCF.
     */
    public enum AlleleFrequencyDataSource {
        exac, gnomad
    }

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which filtered variants should be written.")
    protected GATKPath outputFile;

    @Argument(
            fullName =  FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME,
            doc = "The version of the Human Genome reference which was used to Funcotate the input VCF."
    )
    protected Reference reference;

    @Argument(
            fullName = FuncotatorArgumentDefinitions.ALLELE_FREQUENCY_DATA_SOURCE_NAME,
            doc = "The allele frequency data source (ExAC or gnomAD) that was used to Funcotate the input VCF."
    )
    protected AlleleFrequencyDataSource afDataSource;

    private VariantContextWriter outputVcfWriter;
    private String[] funcotationKeys;
    private final List<TwoPassFuncotationFilter> firstPassFilters = new ArrayList<>();
    private final List<FuncotationFilter> secondPassFilters = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        final VCFHeader vcfHeader = getHeaderForVariants();

        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        if (funcotationHeaderLine != null) {
            funcotationKeys = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());
            outputVcfWriter = createVCFWriter(outputFile);
            vcfHeader.addMetaDataLine(new VCFFilterHeaderLine(FilterFuncotationsConstants.NOT_CLINSIG_FILTER,
                    FilterFuncotationsConstants.NOT_CLINSIG_FILTER_DESCRIPTION));
            vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(FilterFuncotationsConstants.CLINSIG_INFO_KEY, 1,
                    VCFHeaderLineType.String, FilterFuncotationsConstants.CLINSIG_INFO_KEY_DESCRIPTION));
            outputVcfWriter.writeHeader(vcfHeader);
        } else {
            throw new UserException.BadInput("Could not extract Funcotation keys from " +
                    VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME + " field in input VCF header.");
        }

        registerFilters();
    }

    private void registerFilters() {
        FuncotationFilter clinVarFilter = new ClinVarFilter(afDataSource);
        FuncotationFilter lofFilter = new LofFilter(reference, afDataSource);
        FuncotationFilter lmmFilter = new LmmFilter();
        FuncotationFilter homvarFilter = new ArHomvarFilter(reference);
        TwoPassFuncotationFilter hetvarFilter = new ArHetvarFilter(reference, funcotationKeys);
        firstPassFilters.add(hetvarFilter);
        secondPassFilters.add(clinVarFilter);
        secondPassFilters.add(lofFilter);
        secondPassFilters.add(lmmFilter);
        secondPassFilters.add(homvarFilter);
        secondPassFilters.add(hetvarFilter);

    }

    @Override
    public void firstPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        for (TwoPassFuncotationFilter filter: firstPassFilters) {
            filter.firstPassApply(variant);
        }
    }

    @Override
    protected void afterFirstPass() {
        for (TwoPassFuncotationFilter filter: firstPassFilters) {
            filter.afterFirstPass();
        }
    }

    @Override
    public void secondPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final Set<String> matchingFilters = getMatchingFilters(variant, secondPassFilters);
        outputVcfWriter.add(applyFilters(variant, matchingFilters));
    }

    /**
     * Collect the names of the {@link FuncotationFilter}s matching the Funcotations of the given variant.
     *
     * The filter will be treated as a match if it matches Funcotations for any of the transcripts in the
     * variant's Funcotation map.
     */
    private Set<String> getMatchingFilters(final VariantContext variant, final List<FuncotationFilter> funcotationFilters) {
        final Set<String> matchingFilters = new HashSet<>();


        final Map<Allele, FuncotationMap> funcs = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                funcotationKeys, variant, "Gencode_" + reference.gencodeVersion + "_annotationTranscript", "FAKE_SOURCE");

        funcs.values().forEach(funcotationMap -> {
            FilterFuncotationsUtils.getTranscriptFuncotations(funcotationMap).forEach(funcotations -> {
                final Set<String> matches = funcotationFilters.stream()
                        .filter(f -> f.checkFilter(funcotations, variant))
                        .map(FuncotationFilter::getFilterName)
                        .collect(Collectors.toSet());
                matchingFilters.addAll(matches);
            });
        });

        return matchingFilters;
    }

    /**
     * Mark a variant as matching a set of Funcotation filters, or as matching no filters.
     */
    private VariantContext applyFilters(final VariantContext variant, final Set<String> matchingFilters) {
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(variant);
        final boolean isSignificant = !matchingFilters.isEmpty();

        // Mark the individual filters that make the variant significant, if any.
        final String clinicalSignificance = isSignificant ?
                String.join(FilterFuncotationsConstants.FILTER_DELIMITER, matchingFilters) :
                FilterFuncotationsConstants.CLINSIG_INFO_NOT_SIGNIFICANT;
        variantContextBuilder.attribute(FilterFuncotationsConstants.CLINSIG_INFO_KEY, clinicalSignificance);

        // Also set the filter field for insignificant variants, to make it easier for
        // downstream tools to extract out the interesting data.
        if (isSignificant) {
            variantContextBuilder.passFilters();
        } else {
            variantContextBuilder.filter(FilterFuncotationsConstants.NOT_CLINSIG_FILTER);
        }

        return variantContextBuilder.make();
    }

    @Override
    public void closeTool() {
        if (outputVcfWriter != null) {
            outputVcfWriter.close();
        }
    }
}
