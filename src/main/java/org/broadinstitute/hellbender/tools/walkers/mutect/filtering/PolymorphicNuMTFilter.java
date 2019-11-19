package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang.mutable.MutableBoolean;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.OptionalInt;
import java.util.stream.IntStream;

public class PolymorphicNuMTFilter extends HardFilter {
    private static final double LOWER_BOUND_PROB = .01;
//    private static final double MULTIPLE_COPIES_MULTIPLIER = 2.0;
    private final int maxAltDepthCutoff;

    public PolymorphicNuMTFilter(final double maxNuMTCopies){
        if (maxNuMTCopies != 0) {
            final PoissonDistribution autosomalCoverage = new PoissonDistribution(maxNuMTCopies);
            maxAltDepthCutoff = autosomalCoverage.inverseCumulativeProbability(1 - LOWER_BOUND_PROB);
        } else {
            maxAltDepthCutoff = 0;
        }
    }

    @Override
    public ErrorType errorType() { return ErrorType.NON_SOMATIC; }

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        return vc.getGenotypes().stream().filter(filteringEngine::isTumor)
                .filter(Genotype::hasAD)
                .anyMatch(g -> {
                    final int[] alleleDepths = g.getAD();
                    final int numRealAlleles = vc.hasSymbolicAlleles() ? alleleDepths.length - 1 : alleleDepths.length;
                    //Start at first alternate allele depth (the ref allele is first)
                    final OptionalInt max = IntStream.range(1, numRealAlleles).map(a -> alleleDepths[a]).max();
                    return max.getAsInt() < maxAltDepthCutoff;
                });
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.POTENTIAL_POLYMORPHIC_NUMT_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.emptyList(); }
}
