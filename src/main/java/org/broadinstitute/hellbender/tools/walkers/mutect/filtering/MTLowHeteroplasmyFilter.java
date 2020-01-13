package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.TwoPassVariantWalker;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;

import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.LOW_HET_FILTER_NAME;

@CommandLineProgramProperties(
        summary = "If too many low heteroplasmy sites pass other filters, then filter all low heteroplasmy sites",
        oneLineSummary = "If too many low het sites, filter all low het sites",
        programGroup = VariantFilteringProgramGroup.class
)
public class MTLowHeteroplasmyFilter extends TwoPassVariantWalker {

    public static final String MIN_LOW_HET_SITES_LONG_NAME = "min-low-het-sites";
    public static final String LOW_HET_THRESHOLD_LONG_NAME = "low-het-threshold";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF file")
    private String outputVcf = null;

    @Argument(fullName = MIN_LOW_HET_SITES_LONG_NAME,
            doc = "Number of low het sites allowed to pass other filters before filtering out all low het sites. Default is 5",
            optional=true)
    private int minLowHetSites = 3;

    @Argument(fullName = LOW_HET_THRESHOLD_LONG_NAME,
            doc = "Threshold for determining a low heteroplasmy site. Default is 0.1",
            optional=true)
    private final double lowHetThreshold = 0.1;

    private boolean failedLowHet = false;
    private int unfilteredLowHetSites = 0;

    private VariantContextWriter vcfWriter;

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        // TODO why isn't it being added in the GATKVCFHeaderLines
        inputHeader.addMetaDataLine(new VCFFilterHeaderLine(LOW_HET_FILTER_NAME, "All low heteroplasmy sites are filtered when at least x low het sites pass all other filters"));
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(inputHeader);
    }

    @Override
    protected void firstPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        // if the site is not filtered, but it is low het increment counter
        if (variant.isNotFiltered() && isLowHeteroplasmy(variant)) {
            unfilteredLowHetSites++;
        }
    }

    @Override
    protected void afterFirstPass() {
        failedLowHet = unfilteredLowHetSites > minLowHetSites;
    }

    @Override
    protected void secondPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        if (failedLowHet && isLowHeteroplasmy(variant)) {
            vcb.filter(GATKVCFConstants.LOW_HET_FILTER_NAME);
        }
        vcfWriter.add(vcb.make());
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

    protected boolean isLowHeteroplasmy(VariantContext v) {
        // does 0.0 make sense for orElse?
        return v.getGenotypes().stream().map(g -> lowestAF(g)).min(Double::compareTo).orElse(0.0) < lowHetThreshold;
    }

    protected double lowestAF(Genotype g) {
        int[] depths = g.getAD();
        return MathUtils.arrayMin(depths)/ MathUtils.sum(depths);
    }
}
