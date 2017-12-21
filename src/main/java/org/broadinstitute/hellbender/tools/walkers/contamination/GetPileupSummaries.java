package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Given a sites VCF of common population SNPs and a BAM file, summarizes alt and ref counts along with population allele frequency.
 *
 * <p>
 *     The resulting table is the input for {@link CalculateContamination}.
 *     The sites VCF, e.g. gnomAD resource file, contains population allele frequencies (AF) in the INFO field.
 *     Note the default maximum population allele frequency (--maximumPopulationAlleleFrequency or -maxAF) is set to 0.2,
 *     which limits sites the tool considers to those in the variants resource file that have allele frequencies (AF) of 0.2 or less.
 * </p>
 *
 * <p>The full gnomAD is not necessary.  Rather, a much smaller eight-column sites-only vcf restricted to commonly-variant germline SNPs
 * and containing no INFO field other than allele frequency (AF) gives identical results and runs faster.
 * See the GATK Resource Bundle for an example human file.</p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk GetPileupSummaries \
 *   -I tumor.bam \
 *   -L intervals.list \
 *   -V variants_for_contamination.vcf.gz \
 *   -O pileups.table
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Calculate pileup statistics for inferring contamination",
        oneLineSummary = "Calculate pileup statistics for inferring contamination",
        programGroup = CoverageAnalysisProgramGroup.class)
@BetaFeature
@DocumentedFeature
public class GetPileupSummaries extends MultiVariantWalker {

    public static final String MAX_SITE_AF_LONG_NAME = "maximum-population-allele-frequency";
    public static final String MIN_SITE_AF_LONG_NAME = "minimum-population-allele-frequency";
    public static final String MAX_SITE_AF_SHORT_NAME = "max-af";
    public static final String MIN_SITE_AF_SHORT_NAME = "min-af";

    private static final double DEFAULT_MIN_POPULATION_AF = 0.01;
    private static final double DEFAULT_MAX_POPULATION_AF = 0.2;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table", optional=false)
    private File outputTable;

    @Argument(fullName = MIN_SITE_AF_LONG_NAME,
            shortName = MIN_SITE_AF_SHORT_NAME,
            doc = "Minimum population allele frequency of sites to consider.  A low value increases accuracy at the expense of speed.", optional = true)
    private double minPopulationAlleleFrequency = DEFAULT_MIN_POPULATION_AF;

    @Argument(fullName = MAX_SITE_AF_LONG_NAME,
            shortName = MAX_SITE_AF_SHORT_NAME,
            doc = "Maximum population allele frequency of sites to consider.", optional = true)
    private double maxPopulationAlleleFrequency = DEFAULT_MAX_POPULATION_AF;

    private final List<PileupSummary> pileupSummaries = new ArrayList<>();

    private VariantContext lastVariant = null;

    private boolean sawVariantsWithoutAlleleFrequency = false;
    private boolean sawVariantsWithAlleleFrequency = false;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return false;
    }

    private static final int MAPPING_QUALITY_THRESHOLD = 50;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>();
        filters.add(new MappingQualityReadFilter(MAPPING_QUALITY_THRESHOLD));
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.PRIMARY_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE);
        filters.add(ReadFilterLibrary.GOOD_CIGAR);
        filters.add(new WellformedReadFilter());
        return filters;
    }

    @Override
    public void onTraversalStart() {
        final boolean alleleFrequencyInHeader = getHeaderForVariants().getInfoHeaderLines().stream()
                .anyMatch(line -> line.getID().equals(VCFConstants.ALLELE_FREQUENCY_KEY));
        if (!alleleFrequencyInHeader) {
            throw new UserException.BadInput("Population vcf does not have an allele frequency (AF) info field in its header.");
        }
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        // if we input multiple sources of variants, ignore repeats
        if (lastVariant != null && vc.getStart() == lastVariant.getStart()) {
            return;
        } else if ( vc.isBiallelic() && vc.isSNP() && alleleFrequencyInRange(vc) ) {
            final ReadPileup pileup = GATKProtectedVariantContextUtils.getPileup(vc, readsContext);
            pileupSummaries.add(new PileupSummary(vc, pileup));
        }
        lastVariant = vc;
    }

    @Override
    public Object onTraversalSuccess() {
        if (sawVariantsWithoutAlleleFrequency && !sawVariantsWithAlleleFrequency) {
            throw new UserException.BadInput("No variants in population vcf had an allele frequency (AF) field.");
        }
        PileupSummary.writePileupSummaries(pileupSummaries, outputTable);
        return "SUCCESS";
    }

    private boolean alleleFrequencyInRange(final VariantContext vc) {
        if (!vc.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
            if (!sawVariantsWithoutAlleleFrequency) {
                logger.warn(String.format("Variant context at %s:%d lacks allele frequency (AF) field.", vc.getContig(), vc.getStart()));
                sawVariantsWithoutAlleleFrequency = true;
            }
            return false;
        } else {
            sawVariantsWithAlleleFrequency = true;
            final double alleleFrequency = vc.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, -1.0);
            return minPopulationAlleleFrequency < alleleFrequency && alleleFrequency < maxPopulationAlleleFrequency;
        }
    }
}
