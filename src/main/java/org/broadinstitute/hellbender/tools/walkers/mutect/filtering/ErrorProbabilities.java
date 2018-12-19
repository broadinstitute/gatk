package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public final class ErrorProbabilities {
    private final Map<Mutect2VariantFilter, Double> probabilitiesByFilter;
    private final EnumMap<ErrorType, Double> probabilitiesByType;
    private final double errorProbability;


    public ErrorProbabilities(final List<Mutect2VariantFilter> filters, final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        probabilitiesByFilter = filters.stream().collect(Collectors.toMap(f -> f, f -> f.errorProbability(vc, filteringEngine)));

        probabilitiesByType = new EnumMap<>(ErrorType.class);
        for (final ErrorType type : ErrorType.values()) {
            probabilitiesByType.put(type, 0.0);
        }

        for (final Mutect2VariantFilter filter : filters) {
            final double probability = filter.errorProbability(vc, filteringEngine);
            probabilitiesByFilter.put(filter, probability);
            probabilitiesByType.compute(filter.errorType(), (type,prob) -> Math.max(prob, probability));
        }

        // treat errors of different types as independent
        double trueProbability = 1;
        for (final double errorProb : probabilitiesByType.values()) {
            trueProbability *= (1 - errorProb);
        }

        errorProbability = 1 - trueProbability;
    }

    public double getErrorProbability() { return errorProbability; }
    public double getTechnicalArtifactProbability() { return probabilitiesByType.get(ErrorType.ARTIFACT); }
    public double getNonSomaticProbability() { return probabilitiesByType.get(ErrorType.NON_SOMATIC); }
    public Map<Mutect2VariantFilter, Double> getProbabilitiesByFilter() { return probabilitiesByFilter; }


}
