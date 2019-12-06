package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.*;
import java.util.stream.Collectors;

public abstract class HardAlleleFilter<T> extends Mutect2AlleleFilter<T> {
    @Override
    public List<Double> calculateErrorProbabilityForAlleles(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        List<Boolean>  alleleArtifacts = areAllelesArtifacts(vc, filteringEngine, referenceContext);
        // only set values for alleles returned
        return alleleArtifacts.stream().map(value -> value ? 1.0 : 0.0).collect(Collectors.toList());
    }

    public abstract List<Boolean> areAllelesArtifacts(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext);

    // the posterior of a hard filter is 0 or 1, hence there's no reason to annotate it
    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.empty();
    }
}
