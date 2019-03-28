package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang.mutable.MutableBoolean;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.IntStream;

public class MinAlleleFractionFilter extends HardFilter {
    private final double minAf;

    public MinAlleleFractionFilter(final double minAf) {
        this.minAf = minAf;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        return vc.getGenotypes().stream().filter(filteringEngine::isTumor)
                .filter(g -> g.hasExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY))
                .anyMatch(g -> {
                    final double[] alleleFractions = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(g, GATKVCFConstants.ALLELE_FRACTION_KEY, () -> null, 1.0);
                    final int numRealAlleles = vc.hasSymbolicAlleles() ? alleleFractions.length - 1 : alleleFractions.length;
                    final OptionalDouble max = IntStream.range(0, numRealAlleles).mapToDouble(a -> alleleFractions[a]).max();
                    return max.getAsDouble() < minAf;
                });
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.emptyList(); }
}
