package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Simple allele counter for SVs. Supports multi-allelic variants, except for multi-allelic CNVs that lack genotype alleles.
 */
public class SVAlleleCounter {

    private final List<Allele> alleles;
    private final int[] alleleCounts;
    private final int alleleNumber;

    protected SVAlleleCounter(final List<Allele> alleles,
                           final int[] alleleCounts,
                           final int alleleNumber) {
        this.alleles = alleles;
        this.alleleCounts = alleleCounts;
        this.alleleNumber = alleleNumber;
    }

    public List<Allele> getAlleles() {
        return alleles;
    }

    public int[] getAlleleCounts() {
        return alleleCounts;
    }

    public int getAlleleNumber() {
        return alleleNumber;
    }

    /**
     * Factory method for counting alleles. All calculations except AF are performed here.
     * @param alleles alternate alleles to count
     * @param genotypes genotypes with assigned alleles
     * @return result
     */
    public static SVAlleleCounter create(final List<Allele> alleles, final List<Genotype> genotypes) {
        final Map<Allele, Long> alleleCountsMap = SVCallRecordUtils.getAlleleCounts(genotypes);
        final int alleleNumber = alleleCountsMap.values().stream().mapToInt(Long::intValue).sum();

        final int[] alleleCounts = new int[alleles.size()];
        for (int i = 0; i < alleles.size(); i++) {
            alleleCounts[i] = alleleCountsMap.getOrDefault(alleles.get(i), 0L).intValue();
        }
        return new SVAlleleCounter(alleles, alleleCounts, alleleNumber);
    }

    /**
     * Get allele frequencies (AF). Avoid calling this repeatedly on the same object, as computations are performed
     * on the fly.
     */
    public double[] computeFrequencies() {
        final double[] freq = new double[alleleCounts.length];
        if (alleleNumber == 0) {
            Arrays.fill(freq, Double.NaN);
        } else {
            for (int i = 0; i < alleleCounts.length; i++) {
                freq[i] = alleleCounts[i] / (double) alleleNumber;
            }
        }
        return freq;
    }
}
