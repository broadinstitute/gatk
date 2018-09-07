package org.broadinstitute.hellbender.utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

public final class GenotypeUtils {
    private GenotypeUtils(){}

    /**
     * Returns true of the genotype is a called diploid genotype with likelihoods.
     */
    public static boolean isDiploidWithLikelihoods(final Genotype g) {
        return Utils.nonNull(g).isCalled() && g.hasLikelihoods() && g.getPloidy() == 2;
    }

    /**
     * Returns a triple of ref/het/hom genotype "counts".
     *
     * The exact meaning of the count is dependent on the rounding behavior.
     * if {@code roundContributionFromEachGenotype}: the counts are discrete integer counts of the most probable genotype for each {@link Genotype}
     * else: they are the sum of the normalized likelihoods of each genotype and will not be integers
     *
     * Skips non-diploid genotypes.
     *
     *
     * @param vc the VariantContext that the {@link Genotype}s originated from, non-null
     * @param genotypes a GenotypesContext containing genotypes to count, these must be a subset of {@code vc.getGenotypes()}, non-null
     * @param roundContributionFromEachGenotype if this is true, the normalized likelihood from each genotype will be rounded before
     *                                          adding to the total count
     */
    public static GenotypeCounts computeDiploidGenotypeCounts(final VariantContext vc, final GenotypesContext genotypes,
                                                              final boolean roundContributionFromEachGenotype){
        Utils.nonNull(vc, "vc");
        Utils.nonNull(genotypes, "genotypes");
        final boolean doMultiallelicMapping = !vc.isBiallelic();

        int idxAA = 0, idxAB = 1, idxBB = 2;

        double refCount = 0;
        double hetCount = 0;
        double homCount = 0;

        for (final Genotype g : genotypes) {
            if (! isDiploidWithLikelihoods(g)){
                continue;
            }

            // Genotype::getLikelihoods returns a new array, so modification in-place is safe
            final double[] normalizedLikelihoods = MathUtils.normalizeFromLog10ToLinearSpace(g.getLikelihoods().getAsVector());


            if (doMultiallelicMapping) {
                if (g.isHetNonRef()) {
                    //all likelihoods go to homCount
                    homCount++;
                    continue;
                }

                //get alternate allele for each sample
                final Allele a1 = g.getAllele(0);
                final Allele a2 = g.getAllele(1);
                if (a2.isNonReference()) {
                    final int[] idxVector = vc.getGLIndicesOfAlternateAllele(a2);
                    idxAA = idxVector[0];
                    idxAB = idxVector[1];
                    idxBB = idxVector[2];
                }
                //I expect hets to be reference first, but there are no guarantees (e.g. phasing)
                else if (a1.isNonReference()) {
                    final int[] idxVector = vc.getGLIndicesOfAlternateAllele(a1);
                    idxAA = idxVector[0];
                    idxAB = idxVector[1];
                    idxBB = idxVector[2];
                }
            }

            if( roundContributionFromEachGenotype ) {
                refCount += MathUtils.fastRound(normalizedLikelihoods[idxAA]);
                hetCount += MathUtils.fastRound(normalizedLikelihoods[idxAB]);
                homCount += MathUtils.fastRound(normalizedLikelihoods[idxBB]);
            } else {
                refCount += normalizedLikelihoods[idxAA];
                hetCount += normalizedLikelihoods[idxAB];
                homCount += normalizedLikelihoods[idxBB];
            }
        }
        return new GenotypeCounts(refCount, hetCount, homCount);
    }
}
