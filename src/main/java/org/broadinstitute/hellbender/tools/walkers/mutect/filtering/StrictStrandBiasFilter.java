package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasBySample;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class StrictStrandBiasFilter extends HardFilter { //HardAlleleFilter<Integer> {
    private final int minReadsOnEachStrand;

    public StrictStrandBiasFilter(final int minReadsOnEachStrand) {
        this.minReadsOnEachStrand = minReadsOnEachStrand;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    public Predicate<Genotype> checkPreconditions() {
        return g -> g.hasExtendedAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
    }

    public List<Integer> getData(Genotype g) {
        int[] data = GATKProtectedVariantContextUtils.getAttributeAsIntArray(g, GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, () -> null, 0);
        return Arrays.stream(data).boxed().collect(Collectors.toList());
    }

//    @Override
//    public List<Boolean> areAllelesArtifacts(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
//        final MutableInt altForwardCount = new MutableInt(0);
//        final MutableInt altReverseCount = new MutableInt(0);

//        LinkedHashMap<Allele, List<Integer>> dataByAllele = getDataByAllele(vc, checkPreconditions(), this::getData, filteringEngine);
//        return dataByAllele.entrySet().stream()
//                .filter(entry -> !entry.getKey().isSymbolic() && !vc.getReference().equals(entry.getKey()))
//                .map(entry -> minReadsOnEachStrand > 0 && entry.getValue().stream().min(Integer::compare).orElse(0) < minReadsOnEachStrand).collect(Collectors.toList());


//        vc.getGenotypes().stream().filter(filteringEngine::isTumor)
//                .filter(g -> g.hasExtendedAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY))
//                .forEach(g -> {
//                    final int[] strandBiasCounts = GATKProtectedVariantContextUtils.getAttributeAsIntArray(g, GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, () -> null, 0);
//                    altForwardCount.add(StrandBiasBySample.getAltForwardCountFromFlattenedContingencyTable(strandBiasCounts));
//                    altReverseCount.add(StrandBiasBySample.getAltReverseCountFromFlattenedContingencyTable(strandBiasCounts));
//                });
//
//    // filter if there is no alt evidence in the forward or reverse strand
//        return Math.min(altForwardCount.getValue(), altReverseCount.getValue()) < minReadsOnEachStrand;
//}

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        if (minReadsOnEachStrand == 0) {
            return false;
        }

        final MutableInt altForwardCount = new MutableInt(0);
        final MutableInt altReverseCount = new MutableInt(0);

        vc.getGenotypes().stream().filter(filteringEngine::isTumor)
                .filter(g -> g.hasExtendedAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY))
                .forEach(g -> {
                    final int[] strandBiasCounts = VariantContextGetters.getAttributeAsIntArray(g, GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, () -> null, 0);
                    altForwardCount.add(StrandBiasBySample.getAltForwardCountFromFlattenedContingencyTable(strandBiasCounts));
                    altReverseCount.add(StrandBiasBySample.getAltReverseCountFromFlattenedContingencyTable(strandBiasCounts));
                });

        // filter if there is no alt evidence in the forward or reverse strand
        return Math.min(altForwardCount.getValue(), altReverseCount.getValue()) < minReadsOnEachStrand;
    }


    @Override
    public String filterName() {
        return GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.emptyList(); }
}
