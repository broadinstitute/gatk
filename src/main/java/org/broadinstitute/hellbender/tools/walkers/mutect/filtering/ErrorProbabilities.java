package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static java.util.stream.Collectors.*;

public final class ErrorProbabilities {
    private final LinkedHashMap<Mutect2Filter, List<Double>> alleleProbabilitiesByFilter;
    private final Map<ErrorType, List<Double>> probabilitiesByTypeAndAllele;
    private final List<Double> combinedErrorProbabilitiesByAllele;
    private final int numAltAlleles;


    public ErrorProbabilities(final List<Mutect2Filter> filters, final VariantContext vc, final Mutect2FilteringEngine filteringEngine, final ReferenceContext referenceContext) {
        numAltAlleles = vc.getAlternateAlleles().size();
        alleleProbabilitiesByFilter = filters.stream()
                .collect(toMap(
                        Function.identity(),
                        f -> f.errorProbabilities(vc, filteringEngine, referenceContext),
                        (a, b) -> a, LinkedHashMap::new))
                // remove filters that were not applied. i.e. returned empty list
                .entrySet().stream().filter(entry -> !entry.getValue().isEmpty())
                .collect(toMap(Map.Entry::getKey, Map.Entry::getValue, (a, b) -> a, LinkedHashMap::new));
        LinkedHashMap<ErrorType, List<List<Double>>> probabilitiesByAllelesForEachFilter = alleleProbabilitiesByFilter.entrySet().stream().collect(
                groupingBy(entry -> entry.getKey().errorType(), LinkedHashMap::new, mapping(entry -> entry.getValue(), toList())));
        probabilitiesByAllelesForEachFilter.replaceAll((k, v) -> ErrorProbabilities.transpose(v));

        probabilitiesByTypeAndAllele = probabilitiesByAllelesForEachFilter.entrySet().stream().collect(toMap(
                Map.Entry::getKey,
                entry -> entry.getValue().stream().map(alleleList -> alleleList.stream().max(Double::compare).orElse(0.0)).collect(Collectors.toList()),
                (a,b) -> a, HashMap::new));


        // treat errors of different types as independent
        combinedErrorProbabilitiesByAllele = transpose(probabilitiesByTypeAndAllele.values().stream().collect(toList()))
                .stream().map(
                        alleleProbabilities -> alleleProbabilities.stream().map(p -> 1.0 - p).reduce(1.0, (a, b) -> a * b)).collect(Collectors.toList());
        combinedErrorProbabilitiesByAllele.replaceAll(trueProb -> Mutect2FilteringEngine.roundFinitePrecisionErrors(1.0 - trueProb));
    }

    public List<Double> getCombinedErrorProbabilities() { return combinedErrorProbabilitiesByAllele; }
    public List<Double> getTechnicalArtifactProbabilities() { return probabilitiesByTypeAndAllele.get(ErrorType.ARTIFACT); }
    public List<Double> getNonSomaticProbabilities() { return probabilitiesByTypeAndAllele.get(ErrorType.NON_SOMATIC); }
    public Map<Mutect2Filter, List<Double>> getProbabilitiesByFilter() { return alleleProbabilitiesByFilter; }

    public Map<Mutect2Filter, List<Double>> getProbabilitiesForAlleleFilters() {
        return getPartitionedProbabilitiesByFilter(false);
    }

    public Map<Mutect2Filter, Double> getProbabilitiesForVariantFilters() {
        return getPartitionedProbabilitiesByFilter(false).entrySet().stream()
                .filter(entry -> entry.getValue() != null && !entry.getValue().isEmpty())
                .collect(toMap(entry -> entry.getKey(), entry -> entry.getValue().get(0)));
    }

    private Map<Mutect2Filter, List<Double>> getPartitionedProbabilitiesByFilter(boolean variantOnly) {
        Map<Boolean, LinkedHashMap<Mutect2Filter, List<Double>>> groups =
                alleleProbabilitiesByFilter.entrySet().stream().collect(Collectors.partitioningBy(
                        entry -> entry.getKey().getClass().isInstance(Mutect2VariantFilter.class),
                        toMap(Map.Entry::getKey, Map.Entry::getValue, (a,b) -> a, LinkedHashMap::new)));
        return groups.get(variantOnly);
    }

    // TODO would this be useful in a util class somewhere?
    private static <T> List<List<T>> transpose(List<List<T>> list) {
        final int N = list.stream().mapToInt(l -> l.size()).max().orElse(-1);
        List<Iterator<T>> iterList = list.stream().map(it->it.iterator()).collect(toList());
        return IntStream.range(0, N)
                .mapToObj(n -> iterList.stream()
                        .filter(it -> it.hasNext())
                        .map(m -> m.next())
                        .collect(toList()))
                .collect(toList());
    }}
