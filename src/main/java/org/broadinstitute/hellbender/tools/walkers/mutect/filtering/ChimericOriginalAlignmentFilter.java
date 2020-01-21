package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class ChimericOriginalAlignmentFilter extends HardAlleleFilter {
    private final double maxNuMTFraction;

    public ChimericOriginalAlignmentFilter(final double maxNuMTFraction) {
        this.maxNuMTFraction = maxNuMTFraction;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    public Predicate<Genotype> checkPreconditions() {
        return Genotype::hasAD;
    }

    public List<Integer> getData(Genotype g) {
        return Arrays.stream(g.getAD()).boxed().collect(Collectors.toList());
    }

    @Override
    public List<Boolean> areAllelesArtifacts(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        if(!vc.isBiallelic()) {
            return Collections.emptyList();
        }
        final int nonMitochondrialOriginalAlignmentCount = vc.getAttributeAsInt(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY, 0);
        LinkedHashMap<Allele, List<Integer>> dataByAllele = getDataByAllele(vc, checkPreconditions(), this::getData, filteringEngine);
        return dataByAllele.entrySet().stream()
                .filter(entry -> !vc.getReference().equals(entry.getKey()))
                .map(entry -> (double) nonMitochondrialOriginalAlignmentCount / entry.getValue().stream().mapToInt(Integer::intValue).sum() > maxNuMTFraction).collect(Collectors.toList());
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.CHIMERIC_ORIGINAL_ALIGNMENT_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY); }
}
