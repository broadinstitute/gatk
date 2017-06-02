package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Given a sites vcf of common SNPs and a bam file, get summary of alt and ref counts along with population allele frequency.
 *
 * The output of this tool is the input for {@link CalculateContamination}.  Use them as follows:
 *
 * java -jar gatk.jar GetPileupSummaries -V ExAC.vcf -V thousand_genomes.vcf -I input.bam -O pileups.table
 * java -jar gatk.jar CalculateContamination -I pileups.table -O contamination.table
 *
 * Created by David Benjamin on 2/15/17.
 */
@CommandLineProgramProperties(
        summary = "Calculate pileup statistics for inferring contamination",
        oneLineSummary = "Calculate pileup statistics for inferring contamination",
        programGroup = VariantProgramGroup.class)
public class GetPileupSummaries extends MultiVariantWalker {

    public final String MAX_SITE_AF_LONG_NAME = "maximumPopulationAlleleFrequency";
    public final String MIN_SITE_AF_LONG_NAME = "minimumPopulationAlleleFrequency";
    public final String MAX_SITE_AF_SHORT_NAME = "maxAF";
    public final String MIN_SITE_AF_SHORT_NAME = "minAF";

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

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return false;
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
        PileupSummary.writePileupSummaries(pileupSummaries, outputTable);
        return "SUCCESS";
    }

    private boolean alleleFrequencyInRange(final VariantContext vc) {
        if (!vc.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
            return false;
        } else {
            final double alleleFrequency = vc.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, -1.0);
            return minPopulationAlleleFrequency < alleleFrequency && alleleFrequency < maxPopulationAlleleFrequency;
        }
    }
}
