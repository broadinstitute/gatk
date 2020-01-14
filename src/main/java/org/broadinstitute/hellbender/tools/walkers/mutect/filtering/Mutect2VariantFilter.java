package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public abstract class Mutect2VariantFilter extends Mutect2Filter {
    public Mutect2VariantFilter() { }

    @Override
    public List<Double> errorProbabilities(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        int numAltAlleles = vc.getNAlleles() - 1;
        final double result = requiredAnnotations().stream().allMatch(vc::hasAttribute) ? calculateErrorProbability(vc, filteringEngine, referenceContext) : 0.0;
        ArrayList<Double> resultList = new ArrayList<>(numAltAlleles);
        Collections.fill(resultList, Mutect2FilteringEngine.roundFinitePrecisionErrors(result));
        return resultList;

    }

    protected abstract double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext);
}
