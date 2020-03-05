package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = VariantFilteringProgramGroup.class
)

public class NuMTFilterTool extends VariantWalker {
    public static final String MEDIAN_AUTOSOMAL_COVERAGE_LONG_NAME = "autosomal-coverage";
    public static final String MAX_NUMT_COPIES_IN_AUTOSOME_LONG_NAME = "max-numt-autosomal-copies";
    private static final double DEFAULT_MEDIAN_AUTOSOMAL_COVERAGE = 0;
    private static final double DEFAULT_MAX_NUMT_AUTOSOMAL_COPIES = 4;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF file")
    private String outputVcf = null;

    @Argument(fullName = MEDIAN_AUTOSOMAL_COVERAGE_LONG_NAME, optional = true, doc = "Median autosomal coverage for filtering potential NuMTs when calling on mitochondria.")
    public double medianAutosomalCoverage = DEFAULT_MEDIAN_AUTOSOMAL_COVERAGE;

    @Argument(fullName = MAX_NUMT_COPIES_IN_AUTOSOME_LONG_NAME, optional = true, doc = "Max expected NUMT copies in autosome used for filtering potential NuMTs when calling on mitochondria.")
    public double maxNuMTAutosomalCopies = DEFAULT_MAX_NUMT_AUTOSOMAL_COPIES;


    private VariantContextWriter vcfWriter;
    private static final double LOWER_BOUND_PROB = .01;
    private int maxAltDepthCutoff = 0;


    @Override
    public void onTraversalStart() {
        final VCFHeader header = getHeaderForVariants();
        header.addMetaDataLine(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.POSSIBLE_NUMT_FILTER_NAME));
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(header);
        if (maxNuMTAutosomalCopies > 0 && medianAutosomalCoverage > 0) {
            final PoissonDistribution autosomalCoverage = new PoissonDistribution(medianAutosomalCoverage * maxNuMTAutosomalCopies / 2.0);
            maxAltDepthCutoff = autosomalCoverage.inverseCumulativeProbability(1 - LOWER_BOUND_PROB);
        }
    }

    public List<Integer> getData(Genotype g) {
        return Arrays.stream(g.getAD()).boxed().collect(Collectors.toList());
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        LinkedHashMap<Allele, List<Integer>> dataByAllele = Mutect2AlleleFilter.getDataByAllele(variant, Genotype::hasAD, this::getData, null);
        List<Boolean> appliedFilter = dataByAllele.entrySet().stream()
                .filter(entry -> !variant.getReference().equals(entry.getKey()))
                .map(entry -> entry.getValue().stream().max(Integer::compare).orElse(0) < maxAltDepthCutoff).collect(Collectors.toList());
        if (!appliedFilter.contains(Boolean.FALSE)) {
            vcb.filter(filterName());
        }
        if (appliedFilter.contains(Boolean.TRUE)) {
            vcb.attribute(GATKVCFConstants.AS_FILTER_STATUS_KEY, AlleleFilterUtils.getMergedASFilterString(variant, appliedFilter, filterName()));
        }
        vcfWriter.add(vcb.make());
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

    public String filterName() {
        return GATKVCFConstants.POSSIBLE_NUMT_FILTER_NAME;
    }

}
