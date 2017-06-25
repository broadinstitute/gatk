package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.tools.exome.FilterByOrientationBias;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Filter SNVs and indels from a Mutect2 callset.
 *
 * <p>
 *     FilterMutectCalls encapsulates GATK3 MuTect2's filtering functionality.
 *     GATK4 Mutect2 retains variant calling and some prefiltering.
 *     Separating calling and filtering functionalities into two tools better enables an iterative filtering process
 *     that allows for context-specific optimizations. To filter further based on sequence context artifacts,
 *     additionally use {@link FilterByOrientationBias}.
 * </p>
 *
 * <p>
 *     Filtering thresholds for both normal_artifact_lod (default threshold 0.0) and tumor_lod (default threshold 5.3) can be set in this tool.
 *     If the normal artifact log odds is larger than the threshold, then FilterMutectCalls applies the artifact-in-normal filter.
 *     For matched normal analyses with tumor contamination in the normal, consider increasing the normal_artifact_lod threshold.
 *     If the tumor log odds is smaller than the threshold, then FilterMutectCalls filters the variant.
 * </p>
 *
 * <p>
 *     If given a --contaminationTable file, e.g. results from
 *     {@link org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination}, the tool will additionally
 *     filter on contamination fractions. Alternatively, provide a numerical fraction to filter with --contamination_fraction_to_filter.
 * </p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" FilterMutectCalls \
 *   -V tumor_matched_m2_snvs_indels.vcf.gz \
 *   -contaminationTable contamination.table \
 *   -O tumor_matched_m2_filtered.vcf.gz
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Filter somatic SNVs and indels called by Mutect2",
        oneLineSummary = "Filter somatic SNVs and indels called by Mutect2",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class FilterMutectCalls extends VariantWalker {

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file", optional=false)
    private final String outputVcf = null;

    @ArgumentCollection
    protected M2FiltersArgumentCollection MTFAC = new M2FiltersArgumentCollection();

    private VariantContextWriter vcfWriter;

    private final List<VariantContext> unfilteredCalls = new ArrayList<>();


    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = inputHeader.getMetaDataInSortedOrder().stream()
                .filter(line -> !line.getKey().equals(Mutect2FilteringEngine.FILTERING_STATUS_VCF_KEY)) //remove header line from Mutect2 stating that calls are unfiltered.
                .collect(Collectors.toSet());
        headerLines.add(new VCFHeaderLine(Mutect2FilteringEngine.FILTERING_STATUS_VCF_KEY, "These calls have been filtered by " + FilterMutectCalls.class.getSimpleName() + " to label false positives with a list of failed filters and true positives with PASS."));

        GATKVCFConstants.MUTECT_FILTER_NAMES.stream().map(GATKVCFHeaderLines::getFilterLine).forEach(headerLines::add);

        headerLines.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public Object onTraversalSuccess() {
        final String tumorSample = getHeaderForVariants().getMetaDataLine(Mutect2Engine.TUMOR_SAMPLE_KEY_IN_VCF_HEADER).getValue();
        final Mutect2FilteringEngine filteringEngine = new Mutect2FilteringEngine(MTFAC, tumorSample);
        // TODO: implement sophisticated filtering
        for (final VariantContext vc : unfilteredCalls) {
            final VariantContextBuilder vcb = new VariantContextBuilder(vc);
            vcb.filters(filteringEngine.calculateFilters(MTFAC, vc));
            vcfWriter.add(vcb.make());
        }
        return "SUCCESS";
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        unfilteredCalls.add(vc);
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
