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
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.TwoPassVariantWalker;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
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
    private String outputVcf = null;

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
        vcfWriter = createVCFWriter(new File(outputVcf));
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

    public List<Integer> getData(Genotype g) {
        return Arrays.stream(g.getAD()).boxed().collect(Collectors.toList());
    }

    protected boolean isSiteLowHeteroplasmy(VariantContext v) {
        return v.getGenotypes().stream().map(g -> lowestAF(g)).min(Double::compareTo).orElse(0.0) < lowHetThreshold;
    }

    protected List<Boolean> areAllelesArtifacts(final VariantContext vc) {
        VariantContextBuilder vcb = new VariantContextBuilder(vc);
        LinkedHashMap<Allele, List<Integer>> dataByAllele = Mutect2AlleleFilter.getDataByAllele(vc, Genotype::hasAD, this::getData, null);
        Integer total = dataByAllele.values().stream().map(alleleCounts -> alleleCounts.stream().max(Integer::compareTo).orElse(0)).mapToInt(Integer::intValue).sum();
        return dataByAllele.entrySet().stream()
                .filter(entry -> !vc.getReference().equals(entry.getKey()))
                .map(entry -> (entry.getValue().stream().max(Integer::compareTo).orElse(0) / (double) total) < lowHetThreshold).collect(Collectors.toList());
    }

    protected double lowestAF(Genotype g) {
        List<Integer> depths = getData(g);
        return Collections.min(depths.subList(1, depths.size())) / (double) depths.stream().mapToInt(Integer::intValue).sum();
    }
}
