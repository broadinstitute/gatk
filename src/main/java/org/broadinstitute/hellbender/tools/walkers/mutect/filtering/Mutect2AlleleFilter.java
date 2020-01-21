package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public abstract class Mutect2AlleleFilter extends Mutect2Filter {

    public static <T> LinkedHashMap<Allele, List<T>> getDataByAllele(final VariantContext vc, Predicate<Genotype> preconditions, Function<Genotype, List<T>> getData, final Mutect2FilteringEngine filteringEngine) {
        // create and initialize a map with all the alleles in the vc as keys and new, empty lists as values
        LinkedHashMap<Allele, List<T>> dataByAllele = vc.getAlleles().stream().collect(Collectors.toMap(Function.identity(), allele -> new ArrayList<>(), (a, b) -> a, () -> new LinkedHashMap<>()));
        return combineDataByAllele(dataByAllele, vc, preconditions, getData, filteringEngine);
    }

    public static <T> LinkedHashMap<Allele, List<T>> getAltDataByAllele(final VariantContext vc, Predicate<Genotype> preconditions, Function<Genotype, List<T>> getAltData, final Mutect2FilteringEngine filteringEngine) {
        // create and initialize a map with all the alleles in the vc as keys and new, empty lists as values
        LinkedHashMap<Allele, List<T>> dataByAllele = vc.getAlternateAlleles().stream().collect(Collectors.toMap(Function.identity(), allele -> new ArrayList<>(), (a, b) -> a, () -> new LinkedHashMap<>()));
        return combineDataByAllele(dataByAllele, vc, preconditions, getAltData, filteringEngine);
    }

    private static <T> LinkedHashMap<Allele, List<T>> combineDataByAllele(final LinkedHashMap<Allele, List<T>> dataByAllele, final VariantContext vc, Predicate<Genotype> preconditions, Function<Genotype, List<T>> getData, final Mutect2FilteringEngine filteringEngine) {

        // pull all the allele specific data out of each genotype (that passes preconditions) and add it to the list for the allele
        vc.getGenotypes().stream().filter(preconditions).filter(filteringEngine::isTumor)
                .forEach(g -> {
                    Iterator<T> alleleDataIterator = getData.apply(g).iterator();
                    Iterator<List<T>> dataByAlleleIterator = dataByAllele.values().iterator();
                    while(alleleDataIterator.hasNext() && dataByAlleleIterator.hasNext())
                        dataByAlleleIterator.next().add(alleleDataIterator.next());

                });

        return dataByAllele;
    }


    /**
     *
     * @param vc
     * @param filteringEngine
     * @param referenceContext
     * @return The probability that each alternate allele should be filtered out. This list should NOT include data for the reference allele
     */
    @Override
    public List<Double> errorProbabilities(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        return requiredAnnotations().stream().allMatch(vc::hasAttribute) ?
                calculateErrorProbabilityForAlleles(vc, filteringEngine, referenceContext)
                        .stream().map(prob -> Mutect2FilteringEngine.roundFinitePrecisionErrors(prob)).collect(Collectors.toList()) :
                Collections.<Double>emptyList();
    }

    // returning an empty list means filter is not evaluated
    protected abstract List<Double> calculateErrorProbabilityForAlleles(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext);
}
