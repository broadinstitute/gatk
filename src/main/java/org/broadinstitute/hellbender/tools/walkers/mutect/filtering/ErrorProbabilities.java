package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.Arrays;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public final class ErrorProbabilities {
    private final Map<Mutect2VariantFilter, Double> probabilitiesByFilter;
    private final EnumMap<ErrorType, Double> probabilitiesByType;
    private final double errorProbability;


    public ErrorProbabilities(final List<Mutect2VariantFilter> filters, final VariantContext vc, final Mutect2FilteringEngine filteringEngine, final ReferenceContext referenceContext) {
        probabilitiesByFilter = filters.stream().collect(Collectors.toMap(f -> f, f -> f.errorProbability(vc, filteringEngine, referenceContext)));
        probabilitiesByType = Arrays.stream(ErrorType.values()).collect(Collectors.toMap(v -> v, v -> 0.0, (a,b) -> a, () -> new EnumMap<>(ErrorType.class)));
        filters.forEach(f -> probabilitiesByType.compute(f.errorType(), (type,prob) -> Math.max(prob, probabilitiesByFilter.get(f))));

        // treat errors of different types as independent
        double trueProbability = 1;
        for (final double errorProb : probabilitiesByType.values()) {
            trueProbability *= (1 - errorProb);
        }

        errorProbability = Mutect2FilteringEngine.roundFinitePrecisionErrors(1 - trueProbability);
    }

    public double getErrorProbability() { return errorProbability; }
    public double getTechnicalArtifactProbability() { return probabilitiesByType.get(ErrorType.ARTIFACT); }
    public double getNonSomaticProbability() { return probabilitiesByType.get(ErrorType.NON_SOMATIC); }
    public Map<Mutect2VariantFilter, Double> getProbabilitiesByFilter() { return probabilitiesByFilter; }


}
