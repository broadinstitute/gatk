package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashMap;
import java.util.Map;

/**
 * A class containing utility methods used in the calculation of annotations related to cohort heterozygosity, e.g. InbreedingCoefficient and ExcessHet
 * Stores sample count to make sure we never have to iterate the genotypes more than once.
 */
public final class HeterozygosityCalculator {

    //round the likelihoods to return integer numbers of counts (as doubles)
    private static final boolean RETURN_ROUNDED = false;

    private int sampleCount = -1;
    private Map<Allele, Double> hetCounts;
    private Map<Allele, Double> alleleCounts;
    private final boolean returnRounded = RETURN_ROUNDED;
    private final VariantContext vc;

    /**
     * Create a new HeterozygosityUtils -- a new class should be instantiated for each VariantContext to store data for that VC
     */
    public HeterozygosityCalculator(final VariantContext vc) {
        this.vc = vc;
        doGenotypeCalculations();
    }

    /**
     * Get the count of heterozygotes in vc for a specific altAllele (both reference and non-reference hets, e.g. 1/2)
     */
    private void doGenotypeCalculations() {
        final GenotypesContext genotypes = vc.getGenotypes();
        if (genotypes == null || !vc.isVariant()) {
            return;
        }

        final int numAlleles = vc.getNAlleles();

        sampleCount = 0;
        if (hetCounts == null && alleleCounts == null) {
            hetCounts = new HashMap<>();
            alleleCounts = new HashMap<>();
            for (final Allele a : vc.getAlleles()) {
                if (a.isNonReference()) {
                    hetCounts.put(a, 0.0);
                }
                alleleCounts.put(a, 0.0);
            }

            int idxAB;

            //for each sample
            for (final Genotype g : genotypes) {
                if (g.isCalled() && g.hasLikelihoods() && g.getPloidy() == 2){  // only work for diploid samples
                    sampleCount++;
                } else {
                    continue;
                }

                int altIndex = 0;
                for(final Allele a : vc.getAlternateAlleles()) {
                    //for each alt allele index from 1 to N
                    altIndex++;

                        final double[] normalizedLikelihoods = MathUtils.normalizeFromLog10ToLinearSpace(g.getLikelihoods().getAsVector());
                        if (returnRounded) {
                            MathUtils.applyToArrayInPlace(normalizedLikelihoods, Math::round);
                        }

                        //iterate over the other alleles
                        for (int i = 0; i < numAlleles; i++) {
                            //only add homozygotes to alleleCounts, not hetCounts
                            if (i == altIndex) {
                                final double currentAlleleCounts = alleleCounts.get(a);
                                alleleCounts.put(a, currentAlleleCounts + 2*normalizedLikelihoods[GenotypeLikelihoods.calculatePLindex(altIndex,altIndex)]);
                                continue;
                            }
                            //pull out the heterozygote PL index, ensuring that the first allele index < second allele index
                            idxAB = GenotypeLikelihoods.calculatePLindex(Math.min(i,altIndex), Math.max(i,altIndex));
                            final double aHetCounts = hetCounts.get(a);
                            hetCounts.put(a, aHetCounts + normalizedLikelihoods[idxAB]);
                            final double currentAlleleCounts = alleleCounts.get(a);
                            //these are guaranteed to be hets
                            alleleCounts.put(a, currentAlleleCounts + normalizedLikelihoods[idxAB]);
                            final double refAlleleCounts = alleleCounts.get(vc.getReference());
                            alleleCounts.put(vc.getReference(), refAlleleCounts + normalizedLikelihoods[idxAB]);
                        }
                    //add in ref/ref likelihood
                    final double refAlleleCounts = alleleCounts.get(vc.getReference());
                    alleleCounts.put(vc.getReference(), refAlleleCounts + 2*normalizedLikelihoods[0]);
                }

            }
        }
    }

    /**
     * Get the count of heterozygotes in vc for a specific altAllele (both reference and non-reference hets, e.g. 1/2)
     * @param altAllele the alternate allele of interest
     * @return number of hets
     */
    public double getHetCount(final Allele altAllele) {
        Utils.nonNull(altAllele);
        return hetCounts.containsKey(altAllele)? hetCounts.get(altAllele) : 0;
    }

    public double getAlleleCount(final Allele allele) {
        Utils.nonNull(allele);
        return alleleCounts.containsKey(allele)? alleleCounts.get(allele) : 0;
    }

    public int getSampleCount() {
        return sampleCount;
    }
}
