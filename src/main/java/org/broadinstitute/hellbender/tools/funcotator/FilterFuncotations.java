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
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.ClinVarFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.FuncotationFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.LmmFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.LofFilter;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Filter variants based on clinically-significant Funcotations.
 *
 * This proof-of-concept tool is an example for how to parse and use the VCF output of Funcotator.
 * It's currently hard-coded to look for specific {@link Funcotation}s from:
 * <ul>
 *     <li><a href="http://www.clinvar.com/">ClinVar</a></li>
 *     <li><a href="http://exac.broadinstitute.org/">Exome Aggregation Consortium (ExAC)</a></li>
 *     <li><a href="http://personalizedmedicine.partners.org/laboratory-for-molecular-medicine/">Laboratory for Molecular Medicine (LMM)</a></li>
 * </ul>
 */
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
            "\nThis proof-of-concept tool is an example for how to parse and use the VCF output of Funcotator." +
            "\nCurrently hard-coded to look for specific Funcotations from:" +
            "\n  * ClinVar (http://www.clinvar.com/)" +
            "\n  * Exome Aggregation Consortium (ExAC) (http://exac.broadinstitute.org/)" +
            "\n  * Laboratory for Molecular Medicine (LMM) (http://personalizedmedicine.partners.org/laboratory-for-molecular-medicine/)";

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

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which filtered variants should be written.")
    protected File outputFile;

    @Argument(
            fullName =  FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME,
            doc = "The version of the Human Genome reference which was used to Funcotate the input VCF."
    )
    protected Reference reference;

    private VariantContextWriter outputVcfWriter;
    private String[] funcotationKeys;
    private final List<FuncotationFilter> funcotationFilters = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        registerFilters();
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
    }

    private void registerFilters() {
        funcotationFilters.add(new ClinVarFilter());
        funcotationFilters.add(new LofFilter(reference));
        funcotationFilters.add(new LmmFilter());
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        outputVcfWriter.add(applyFilters(variant, getMatchingFilters(variant)));
    }

    /**
     * Collect the names of the {@link FuncotationFilter}s matching the Funcotations of the given variant.
     *
     * The filter will be treated as a match if it matches Funcotations for any of the transcripts in the
     * variant's Funcotation map.
     */
    private Set<String> getMatchingFilters(final VariantContext variant) {
        final Set<String> matchingFilters = new HashSet<>();


        final Map<Allele, FuncotationMap> funcs = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                funcotationKeys, variant, "Gencode_" + reference.gencodeVersion + "_annotationTranscript", "FAKE_SOURCE");

        funcs.values().forEach(funcotationMap -> {
            final Stream<Map<String, String>> transcriptFuncotations = funcotationMap.getTranscriptList().stream()
                    .map(funcotationMap::get)
                    .map(funcotations -> funcotations.stream()
                            .flatMap(this::extractFuncotationFields)
                            .filter(entry -> entry.getValue() != null && !entry.getValue().isEmpty())
                            .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)));

            transcriptFuncotations.forEach(funcotations -> {
                final Set<String> matches = funcotationFilters.stream()
                        .filter(f -> f.checkFilter(funcotations))
                        .map(FuncotationFilter::getFilterName)
                        .collect(Collectors.toSet());
                matchingFilters.addAll(matches);
            });
        });

        return matchingFilters;
    }

    /**
     * Parse the entries in a Funcotation into a stream of map entries.
     */
    private Stream<Map.Entry<String, String>> extractFuncotationFields(final Funcotation funcotation) {
        return funcotation.getFieldNames().stream()
                .map(name -> new AbstractMap.SimpleEntry<>(name, funcotation.getField(name)));
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
