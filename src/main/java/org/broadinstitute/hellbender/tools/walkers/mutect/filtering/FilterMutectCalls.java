package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.exome.FilterByOrientationBias;
import org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>Filter variants in a Mutect2 VCF callset.</p>
 *
 * <p>
 *     FilterMutectCalls encapsulates GATK3 MuTect2's filtering functionality and adds additional filters.
 *     Thresholds for filters are contained in {@link M2FiltersArgumentCollection} and described in
 *     <a href='https://github.com/broadinstitute/gatk/tree/master/docs/mutect/mutect.pdf' target='_blank'>https://github.com/broadinstitute/gatk/tree/master/docs/mutect/mutect.pdf</a>.
 *     To filter based on sequence context artifacts, see {@link FilterByOrientationBias}.
 * </p>
 * <p>
 *     Filtering thresholds for both normal-artifact-lod (default threshold 0.0) and tumor-lod (default threshold 5.3) can be set in this tool.
 *     If the normal artifact log odds is larger than the threshold, then FilterMutectCalls applies the artifact-in-normal filter.
 *     For matched normal analyses with tumor contamination in the normal, consider increasing the normal-artifact-lod threshold.
 *     If the tumor log odds is smaller than the threshold, then FilterMutectCalls filters the variant.
 * </p>
 * <p>
 *     If given a --contamination-table file, e.g. results from
 *     {@link CalculateContamination}, the tool will additionally
 *     filter on contamination fractions. Alternatively, provide a numerical fraction to filter with the --contamination argument.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk FilterMutectCalls \
 *   -V somatic.vcf.gz \
 *   --contamination-table contamination.table \
 *   -O filtered.vcf.gz
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Filter somatic SNVs and indels called by Mutect2",
        oneLineSummary = "Filter somatic SNVs and indels called by Mutect2",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public final class FilterMutectCalls extends MultiplePassVariantWalker {

    public static final String FILTERING_STATS_LONG_NAME = "filtering-stats";

    public static final String FILTERING_STATUS_VCF_KEY = "filtering_status";

    public static final String FILTERING_STATS_EXTENSION = ".filteringStats.tsv";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName =StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file", optional=false)
    private final String outputVcf = null;

    @Argument(fullName = Mutect2.MUTECT_STATS_SHORT_NAME, doc="The Mutect stats file output by Mutect2", optional=true)
    private final String statsTable = null;

    @Argument(fullName = FILTERING_STATS_LONG_NAME, doc="The output filtering stats file", optional=true)
    private final String filteringStatsOutput = null;

    @ArgumentCollection
    protected M2FiltersArgumentCollection MTFAC = new M2FiltersArgumentCollection();

    private VariantContextWriter vcfWriter;

    private Mutect2FilteringEngine filteringEngine;

    private static final int NUMBER_OF_LEARNING_PASSES = 2;

    @Override
    protected int numberOfPasses() { return NUMBER_OF_LEARNING_PASSES + 1; }

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = inputHeader.getMetaDataInSortedOrder().stream()
                .filter(line -> !line.getKey().equals(FILTERING_STATUS_VCF_KEY)) //remove header line from Mutect2 stating that calls are unfiltered.
                .collect(Collectors.toSet());
        headerLines.add(new VCFHeaderLine(FILTERING_STATUS_VCF_KEY, "These calls have been filtered by " + FilterMutectCalls.class.getSimpleName() + " to label false positives with a list of failed filters and true positives with PASS."));

        GATKVCFConstants.MUTECT_FILTER_NAMES.stream().map(GATKVCFHeaderLines::getFilterLine).forEach(headerLines::add);

        headerLines.addAll(getDefaultToolVCFHeaderLines());

        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(vcfHeader);


        final File mutect2StatsTable = new File(statsTable == null ? drivingVariantFile + Mutect2.DEFAULT_STATS_EXTENSION : statsTable);
        filteringEngine = new Mutect2FilteringEngine(MTFAC, vcfHeader, mutect2StatsTable);
        if (!mutect2StatsTable.exists()) {
            logger.warn("Mutect stats table " + mutect2StatsTable + " not found.  Filtering will proceed without this information.");
        }
    }

    @Override
    protected void nthPassApply(final VariantContext variant,
                                final ReadsContext readsContext,
                                final ReferenceContext referenceContext,
                                final FeatureContext featureContext,
                                final int n) {
        if (n < NUMBER_OF_LEARNING_PASSES) {
            filteringEngine.accumulateData(variant);
        } else if (n == NUMBER_OF_LEARNING_PASSES) {
            vcfWriter.add(filteringEngine.applyFiltersAndAccumulateOutputStats(variant));
        } else {
            throw new GATKException.ShouldNeverReachHereException("This walker should never reach (zero-indexed) pass " + n);
        }
    }

    @Override
    protected void afterNthPass(final int n) {
        if (n < NUMBER_OF_LEARNING_PASSES) {
            filteringEngine.learnParameters();
        } else if (n == NUMBER_OF_LEARNING_PASSES) {
            final File filteringStatsFile = new File(filteringStatsOutput != null ? filteringStatsOutput : outputVcf + FILTERING_STATS_EXTENSION);
            filteringEngine.writeFilteringStats(filteringStatsFile);
        } else {
            throw new GATKException.ShouldNeverReachHereException("This walker should never reach (zero-indexed) pass " + n);
        }
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}
