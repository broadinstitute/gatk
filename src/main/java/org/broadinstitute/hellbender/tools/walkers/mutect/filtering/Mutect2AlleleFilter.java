package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Base class for filters that apply at the allele level. This includes helper functions that
 * convert many lists of data for with data for each allele, to lists of data grouped together by the allele
 */
public abstract class Mutect2AlleleFilter extends Mutect2Filter {

    /**
     * The call to use when the data includes data for the reference
     * @param vc the variant context to get the genotypes from
     * @param preconditions a predicate that must pass for data to be returned
     * @param getData a function that returns the list of data for a genotype. the size must match the number of allele in the dataByAllele map
     * @param filteringEngine
     * @param <T> the type of the data returned
     * @return a map the list of values for each allele
     */
    public static <T> LinkedHashMap<Allele, List<T>> getDataByAllele(final VariantContext vc, Predicate<Genotype> preconditions, Function<Genotype, List<T>> getData, final Mutect2FilteringEngine filteringEngine) {
        // create and initialize a map with all the alleles in the vc as keys and new, empty lists as values
        LinkedHashMap<Allele, List<T>> dataByAllele = vc.getAlleles().stream().collect(Collectors.toMap(Function.identity(), allele -> new ArrayList<>(), (a, b) -> a, () -> new LinkedHashMap<>()));
        return combineDataByAllele(dataByAllele, vc, preconditions, getData, filteringEngine);
    }

    /**
     * The call to use when the data is only for the alternate alleles
     * @param vc the variant context to get the genotypes from
     * @param preconditions a predicate that must pass for data to be returned
     * @param getAltData a function that returns the list of data for a genotype. the size must match the number of allele in the dataByAllele map
     * @param filteringEngine
     * @param <T> the type of the data returned
     * @return a map with the list of values for each allele
     */
    public static <T> LinkedHashMap<Allele, List<T>> getAltDataByAllele(final VariantContext vc, Predicate<Genotype> preconditions, Function<Genotype, List<T>> getAltData, final Mutect2FilteringEngine filteringEngine) {
        // create and initialize a map with all the alleles in the vc as keys and new, empty lists as values
        LinkedHashMap<Allele, List<T>> dataByAllele = vc.getAlternateAlleles().stream().collect(Collectors.toMap(Function.identity(), allele -> new ArrayList<>(), (a, b) -> a, () -> new LinkedHashMap<>()));
        return combineDataByAllele(dataByAllele, vc, preconditions, getAltData, filteringEngine);
    }

    /**
     * Helper function that combines data from multiple genotypes to a list of data for each allele
     * @param dataByAllele the map to return with the list of values separated by allele
     * @param vc the variant context to get the genotypes from
     * @param preconditions a predicate that must pass for data to be returned
     * @param getData a function that returns the list of data for a genotype. the size must match the number of allele in the dataByAllele map
     * @param filteringEngine
     * @param <T> the type of the data returned
     * @return the dataByAllele map filled in with the list of values for each allele
     */
    private static <T> LinkedHashMap<Allele, List<T>> combineDataByAllele(final LinkedHashMap<Allele, List<T>> dataByAllele, final VariantContext vc, Predicate<Genotype> preconditions, Function<Genotype, List<T>> getData, final Mutect2FilteringEngine filteringEngine) {

        // pull all the allele specific data out of each genotype (that passes preconditions) and add it to the list for the allele
        vc.getGenotypes().stream().filter(preconditions)
                .forEach(g -> {
                    Iterator<T> alleleDataIterator = getData.apply(g).iterator();
                    Iterator<List<T>> dataByAlleleIterator = dataByAllele.values().iterator();
                    while(alleleDataIterator.hasNext() && dataByAlleleIterator.hasNext())
                        dataByAlleleIterator.next().add(alleleDataIterator.next());

                });

        return dataByAllele;
    }


    /**
     * Returns a list of probabilities, one for each alternate allele which is the probability that the allele should be filtered out.
     * An empty list is returned when the filter is not/can not be evaluated.
     * @param vc
     * @param filteringEngine
     * @param referenceContext
     * @return The probability that each alternate allele should be filtered out. This list should NOT include data for the reference allele
     */
    @Override
    public List<Double> errorProbabilities(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        return requiredInfoAnnotations().stream().allMatch(vc::hasAttribute) ?
                calculateErrorProbabilityForAlleles(vc, filteringEngine, referenceContext)
                        .stream().map(prob -> Mutect2FilteringEngine.roundFinitePrecisionErrors(prob)).collect(Collectors.toList()) :
                Collections.<Double>emptyList();
    }

    // returning an empty list means filter is not evaluated
    protected abstract List<Double> calculateErrorProbabilityForAlleles(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext);
}
