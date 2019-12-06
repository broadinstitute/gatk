package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.*;
import java.util.function.DoubleSupplier;
import java.util.function.Function;
import java.util.stream.Collectors;

public final class ErrorProbabilities {
    private final Map<Mutect2VariantFilter, Double> probabilitiesByFilter;
    private final LinkedHashMap<Mutect2AlleleFilter<?>, List<Double>> probabilitiesByFilterAndAllele;
    private final EnumMap<ErrorType, Double> probabilitiesByType;
    private final double errorProbability;


    public ErrorProbabilities(final List<Mutect2VariantFilter> variantFilters, final List<Mutect2AlleleFilter<?>> alleleFilters, final VariantContext vc, final Mutect2FilteringEngine filteringEngine, final ReferenceContext referenceContext) {
        probabilitiesByFilter = variantFilters.stream().collect(Collectors.toMap(Function.identity(), f -> f.errorProbability(vc, filteringEngine, referenceContext)));
        probabilitiesByFilterAndAllele = alleleFilters.stream().collect(Collectors.toMap(Function.identity(), f -> f.errorProbability(vc, filteringEngine, referenceContext), (a,b) -> a, () -> new LinkedHashMap<>()));
        probabilitiesByType = Arrays.stream(ErrorType.values()).collect(Collectors.toMap(v -> v, v -> 0.0, (a,b) -> a, () -> new EnumMap<>(ErrorType.class)));
        variantFilters.forEach(f -> probabilitiesByType.compute(f.errorType(), (type,prob) -> Math.max(prob, probabilitiesByFilter.get(f))));
//        alleleFilters.forEach(f -> probabilitiesByType.compute(f.errorType(), (type,prob) ->
//                Math.max(prob,
//                        probabilitiesByFilterAndAllele.get(f).values().stream().filter(d -> !d.isNaN()).max(Double::compare).orElseGet(() -> 0.0))));

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
    public Map<Mutect2AlleleFilter<?>, List<Double>> getProbabilitiesByFilterAndAllele() { return probabilitiesByFilterAndAllele; }
}
