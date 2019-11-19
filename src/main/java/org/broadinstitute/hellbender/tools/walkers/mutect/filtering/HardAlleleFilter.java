package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.*;

public abstract class HardAlleleFilter<T> extends Mutect2AlleleFilter<T> {
    public void calculateErrorProbabilityForAlleles(LinkedHashMap<Allele, Double> results, final LinkedHashMap<Allele, List<T>> dataByAllele, final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        Map<Allele, Boolean>  alleleArtifacts = areAllelesArtifacts(dataByAllele, vc, filteringEngine, referenceContext);
        // only set values for alleles returned
        alleleArtifacts.forEach((key, value) -> results.put(key, value ? 1.0 : 0.0));
    }

    public abstract Map<Allele, Boolean> areAllelesArtifacts(final LinkedHashMap<Allele, List<T>> dataByAllele, final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext);

    // the posterior of a hard filter is 0 or 1, hence there's no reason to annotate it
    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.empty();
    }
}
