package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.stream.Collectors;

public class MinAlleleFractionFilter extends HardAlleleFilter {
    private final double minAf;

    public MinAlleleFractionFilter(final double minAf) {
        this.minAf = minAf;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    public List<Double> getAltData(final Genotype g) {
        double[] data = VariantContextGetters.getAttributeAsDoubleArray(g, GATKVCFConstants.ALLELE_FRACTION_KEY, () -> null, 1.0);
        return Doubles.asList(data);
    }

    @Override
    public List<Boolean> areAllelesArtifacts(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        LinkedHashMap<Allele, List<Double>> dataByAllele = getAltDataByAllele(vc, g -> g.hasExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY) && filteringEngine.isTumor(g), this::getAltData, filteringEngine);
        return dataByAllele.entrySet().stream()
                .filter(entry -> !vc.getReference().equals(entry.getKey()))
                .map(entry -> entry.getValue().stream().max(Double::compare).orElse(1.0) < minAf).collect(Collectors.toList());
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME;
    }

    @Override
    protected List<String> requiredInfoAnnotations() { return Collections.emptyList(); }
}
