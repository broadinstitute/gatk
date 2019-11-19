package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public abstract class Mutect2AlleleFilter<T> extends Mutect2Filter {

    public Map<Allele, Double> applyFilter(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        // create and initialize a map with all the alleles in the vc as keys and new, empty lists as values
        LinkedHashMap<Allele, List<T>> dataByAllele = vc.getAlleles().stream().collect(Collectors.toMap(Function.identity(), allele -> new ArrayList<>(), (a, b) -> a, () -> new LinkedHashMap<>()));

        // pull all the allele specific data out of each genotype (that passes preconditions) and add it to the list for the allele
        vc.getGenotypes().stream().filter(filteringEngine::isTumor)
                .filter(checkPreconditions())
                .forEach(g -> {
                    Iterator<T> alleleDataIterator = getData(g).iterator();
                    Iterator<List<T>> dataByAlleleIterator = dataByAllele.values().iterator();
                    while(alleleDataIterator.hasNext() && dataByAlleleIterator.hasNext())
                        dataByAlleleIterator.next().add(alleleDataIterator.next());

                });

        // construct output map with defaults
        LinkedHashMap<Allele, Double> probabilityByAllele = vc.getAlleles().stream().collect(Collectors.toMap(Function.identity(), allele -> this.getDefaultProbability(), (a, b) -> null, () -> new LinkedHashMap<>()));

        // now invoke the filter giving it the map with the data separated by allele
        calculateErrorProbabilityForAlleles(probabilityByAllele, dataByAllele, vc, filteringEngine, referenceContext);
        return probabilityByAllele;
    }

    /**
     * Subclasses should override if they want a different default probability.
     * Keep in mind that in the final output, NaN is used to determine when indicating that no probability was computed, and . will be the output for those alleles
     * @return the default probability to use
     */
    public Double getDefaultProbability() {
        return Double.NaN;
    }

    /**
     * All subclass filters should implement if they need to verify the data needed exists in the genotype
     * @return a predicate that will be applied to determine which genotypes will be part of the filter
     */
    public abstract Predicate<Genotype> checkPreconditions();

    /**
     * All subclass filters should implement this method to return the necessary data needed to apply the filter
     * @param g the genotype to pull the data from
     * @return A list of per-allele data for each allele in the variant context (the data in the genotype is ordered by the alleles returned from vc.getAlleles()
     */
    public abstract List<T> getData(Genotype g);

    public Map<Allele, Double> errorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        return requiredAnnotations().stream().allMatch(vc::hasAttribute) ?
                applyFilter(vc, filteringEngine, referenceContext) :
                // TODO make sure that somewhere the roundFinitePrecisionErrors is called when not a hard filter
//                        .entrySet().stream().collect(Collectors.toMap(Map.Entry::getKey, entry -> Mutect2FilteringEngine.roundFinitePrecisionErrors(entry.getValue()))) :
                Collections.<Allele, Double>emptyMap();
    }

    protected abstract void calculateErrorProbabilityForAlleles(LinkedHashMap<Allele, Double> results, final LinkedHashMap<Allele, List<T>> dataByAllele, final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext);
}
