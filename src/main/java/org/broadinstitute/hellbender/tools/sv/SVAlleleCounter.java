package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Simple allele counter for SVs. Supports multi-allelic variants, except for multi-allelic CNVs that lack genotype alleles.
 */
public class SVAlleleCounter {

    private final List<Allele> alleles;
    private final int[] counts;
    private final double[] frequencies;
    private final int number;

    public SVAlleleCounter(final List<Allele> alleles,
                           final List<Genotype> genotypes) {
        this.alleles = alleles;
        final Map<Allele, Long> alleleCountsMap = computeCounts(genotypes);
        number = alleleCountsMap.values().stream().mapToInt(Long::intValue).sum();
        counts = new int[alleles.size()];
        for (int i = 0; i < alleles.size(); i++) {
            counts[i] = alleleCountsMap.getOrDefault(alleles.get(i), 0L).intValue();
        }
        this.frequencies = computeFrequencies(counts, number);
    }

    public List<Allele> getAlleles() {
        return alleles;
    }

    public int[] getCounts() {
        return counts;
    }

    public int getNumber() {
        return number;
    }

    public double[] getFrequencies() { return frequencies; }

    /**
     * Counts unique alleles in the given set of genotypes.
     */
    private static Map<Allele, Long> computeCounts(final Collection<Genotype> genotypes) {
        return genotypes.stream().map(Genotype::getAlleles).flatMap(Collection::stream)
                .filter(a -> !(a == null || a.equals(Allele.NO_CALL)))
                .collect(Collectors.groupingBy(a -> a, Collectors.counting()));
    }

    /**
     * Compute allele frequencies (AF) based on counts.
     */
    private static double[] computeFrequencies(final int[] counts, final int number) {
        final double[] freq = new double[counts.length];
        if (number == 0) {
            Arrays.fill(freq, Double.NaN);
        } else {
            for (int i = 0; i < counts.length; i++) {
                freq[i] = counts[i] / (double) number;
            }
        }
        return freq;
    }
}
