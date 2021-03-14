package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.TwoPassVariantWalker;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "If too many low heteroplasmy sites pass other filters, then filter all low heteroplasmy sites",
        oneLineSummary = "If too many low het sites, filter all low het sites",
        programGroup = VariantFilteringProgramGroup.class
)
public class MTLowHeteroplasmyFilterTool extends TwoPassVariantWalker {

    public static final String MAX_ALLOWED_LOW_HETS_LONG_NAME = "max-allowed-low-hets";
    public static final String LOW_HET_THRESHOLD_LONG_NAME = "low-het-threshold";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF file")
    private GATKPath outputVcf = null;

    @Argument(fullName = MAX_ALLOWED_LOW_HETS_LONG_NAME,
            doc = "Number of low het sites allowed to pass other filters before filtering out all low het sites. Default is 3",
            optional=true)
    private int maxAllowedLowHets = 3;

    @Argument(fullName = LOW_HET_THRESHOLD_LONG_NAME,
            doc = "Threshold for determining a low heteroplasmy site. Default is 0.1",
            optional=true)
    private final double lowHetThreshold = 0.1;

    private boolean failedLowHet = false;
    private int unfilteredLowHetSites = 0;

    private VariantContextWriter vcfWriter;

    @Override
    public void onTraversalStart() {
        final VCFHeader header = getHeaderForVariants();
        header.addMetaDataLine(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_HET_FILTER_NAME));
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(header);
    }

    @Override
    protected void firstPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        // if the site is not filtered but it is low het, increment counter
        if (variant.isNotFiltered() && isSiteLowHeteroplasmy(variant)) {
            unfilteredLowHetSites++;
        }
    }

    @Override
    protected void afterFirstPass() {
        failedLowHet = unfilteredLowHetSites > maxAllowedLowHets;
    }

    @Override
    protected void secondPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        if (failedLowHet) {
            List<Boolean> appliedFilter = areAllelesArtifacts(variant);
            if (!appliedFilter.contains(Boolean.FALSE)) {
                vcb.filter(filterName());
            }
            if (appliedFilter.contains(Boolean.TRUE)) {
                vcb.attribute(GATKVCFConstants.AS_FILTER_STATUS_KEY, AlleleFilterUtils.getMergedASFilterString(variant, appliedFilter, filterName()));
            }
        }
        vcfWriter.add(vcb.make());
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

    protected String filterName() {
        return GATKVCFConstants.LOW_HET_FILTER_NAME;
    }

    public List<Double> getData(Genotype g) {
        return Arrays.stream(VariantContextGetters.getAttributeAsDoubleArray(g, GATKVCFConstants.ALLELE_FRACTION_KEY, () -> null, Double.MAX_VALUE)).boxed().collect(Collectors.toList());
    }

    protected boolean isSiteLowHeteroplasmy(VariantContext v) {
        return !(v.getGenotypes().stream().flatMap(g -> getData(g).stream().filter(x -> x < lowHetThreshold)).collect(Collectors.toList()).isEmpty());
    }

    protected List<Boolean> areAllelesArtifacts(final VariantContext vc) {
        LinkedHashMap<Allele, List<Double>> dataByAllele = Mutect2AlleleFilter.getAltDataByAllele(vc, g -> g.hasExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY), this::getData, null);
        return dataByAllele.values().stream()
                .map(afList -> afList.stream().max(Double::compareTo).orElse(0.0) < lowHetThreshold).collect(Collectors.toList());
    }
}
