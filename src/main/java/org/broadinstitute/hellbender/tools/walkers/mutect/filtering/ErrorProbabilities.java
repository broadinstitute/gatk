package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static java.util.stream.Collectors.*;

public final class ErrorProbabilities {
    private  LinkedHashMap<Mutect2Filter, List<Double>> alleleProbabilitiesByFilter;
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

        // if vc has symbolic alleles, remove them from each filter list
        alleleProbabilitiesByFilter.replaceAll((k, v) -> GATKVariantContextUtils.removeDataForSymbolicAltAlleles(vc, v));
        final LinkedHashMap<ErrorType, List<List<Double>>> probabilitiesByAllelesForEachFilter = alleleProbabilitiesByFilter.entrySet().stream().collect(
                groupingBy(entry -> entry.getKey().errorType(), LinkedHashMap::new, mapping(entry -> entry.getValue(), toList())));
        // convert the data so we have a list of probabilities by allele instead of filter
        probabilitiesByAllelesForEachFilter.replaceAll((k, v) -> ErrorProbabilities.transpose(v));

        // foreach error type, get the max probability for each allele
        probabilitiesByTypeAndAllele = probabilitiesByAllelesForEachFilter.entrySet().stream().collect(toMap(
                Map.Entry::getKey,
                entry -> entry.getValue().stream().map(alleleList -> alleleList.stream().max(Double::compare).orElse(0.0)).collect(Collectors.toList()),
                (a,b) -> a, HashMap::new));


        // treat errors of different types as independent
        // transpose the lists of allele probabilities, so it is now a list per allele that contains the prob for each type
        // combine allele-wise
        combinedErrorProbabilitiesByAllele = transpose(probabilitiesByTypeAndAllele.values().stream().collect(toList()))
                .stream().map(
                        alleleProbabilities -> alleleProbabilities.stream().map(p -> 1.0 - p).reduce(1.0, (a, b) -> a * b)).collect(Collectors.toList());
        combinedErrorProbabilitiesByAllele.replaceAll(trueProb -> Mutect2FilteringEngine.roundFinitePrecisionErrors(1.0 - trueProb));
    }

    public List<Double> getCombinedErrorProbabilities() { return combinedErrorProbabilitiesByAllele; }
    public List<Double> getTechnicalArtifactProbabilities() { return probabilitiesByTypeAndAllele.get(ErrorType.ARTIFACT); }
    public List<Double> getNonSomaticProbabilities() { return probabilitiesByTypeAndAllele.get(ErrorType.NON_SOMATIC); }
    public Map<Mutect2Filter, List<Double>> getProbabilitiesByFilter() { return alleleProbabilitiesByFilter; }

    // helper functions for the few operations that still differ depending on whether the filter
    // is per variant or allele
    public LinkedHashMap<Mutect2Filter, List<Double>> getProbabilitiesForAlleleFilters() {
        return getPartitionedProbabilitiesByFilter(false);
    }

    public LinkedHashMap<Mutect2Filter, Double> getProbabilitiesForVariantFilters() {
        return getPartitionedProbabilitiesByFilter(true).entrySet().stream()
                .filter(entry -> entry.getValue() != null && !entry.getValue().isEmpty())
                .collect(toMap(entry -> entry.getKey(), entry -> entry.getValue().get(0), (a,b) -> b, LinkedHashMap::new));
    }

    private LinkedHashMap<Mutect2Filter, List<Double>> getPartitionedProbabilitiesByFilter(boolean variantOnly) {
        final Map<Boolean, LinkedHashMap<Mutect2Filter, List<Double>>> groups =
                alleleProbabilitiesByFilter.entrySet().stream().collect(Collectors.partitioningBy(
                        entry -> Mutect2VariantFilter.class.isAssignableFrom(entry.getKey().getClass()),
                        toMap(Map.Entry::getKey, Map.Entry::getValue, (a,b) -> a, LinkedHashMap::new)));
        return groups.get(variantOnly);
    }

    public static <T> List<List<T>> transpose(List<List<T>> list) {
        // all lists need to be the same size
        if (list.isEmpty()) {
            return list;
        }
        Utils.validateArg(list.stream().map(List::size).distinct().count() == 1, "lists are not the same size");
        final List<Iterator<T>> iterList = list.stream().map(it -> it.iterator()).collect(toList());
        return IntStream.range(0, list.get(0).size())
                .mapToObj(n -> iterList.stream()
                        .filter(it -> it.hasNext())
                        .map(m -> m.next())
                        .collect(toList()))
                .collect(toList());
    }
}
