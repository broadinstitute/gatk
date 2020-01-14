package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang.mutable.MutableBoolean;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public class MinAlleleFractionFilter extends HardAlleleFilter<Double> {
    private final double minAf;

    public MinAlleleFractionFilter(final double minAf) {
        this.minAf = minAf;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    public Predicate<Genotype> checkPreconditions() {
        return g -> g.hasExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY);
    }

    public List<Double> getAltData(Genotype g) {
        double[] data = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(g, GATKVCFConstants.ALLELE_FRACTION_KEY, () -> null, 1.0);
        return Arrays.stream(data).boxed().collect(Collectors.toList());
    }

    @Override
    public List<Boolean> areAllelesArtifacts(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        LinkedHashMap<Allele, List<Double>> dataByAllele = getAltDataByAllele(vc, checkPreconditions(), this::getAltData, filteringEngine);
        return dataByAllele.entrySet().stream()
                .filter(entry -> /*!entry.getKey().isSymbolic() &&*/ !vc.getReference().equals(entry.getKey()))
                .map(entry -> entry.getValue().stream().max(Double::compare).orElse(1.0) < minAf).collect(Collectors.toList());
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.emptyList(); }
}
