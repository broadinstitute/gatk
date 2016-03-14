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
     * Returns a triple of ref/het/hom genotype counts. Skips non-diploid genotypes.
     */
    public static GenotypeCounts computeDiploidGenotypeCounts(final VariantContext vc, final GenotypesContext genotypes){
        Utils.nonNull(vc, "vc");
        Utils.nonNull(genotypes, "genotypes");
        final boolean doMultiallelicMapping = !vc.isBiallelic();

        int idxAA = 0, idxAB = 1, idxBB = 2;

        int refCount = 0;
        int hetCount = 0;
        int homCount = 0;

        for (final Genotype g : genotypes) {
            if (! isDiploidWithLikelihoods(g)){
                continue;
            }

            //Need to round the likelihoods to deal with small numerical deviations due to normalizing
            final double[] normalizedLikelihoodsUnrounded = MathUtils.normalizeFromLog10(g.getLikelihoods().getAsVector());
            final int[] normalizedLikelihoods = new int[normalizedLikelihoodsUnrounded.length];
            for (int i = 0; i < normalizedLikelihoodsUnrounded.length; i++) {
                normalizedLikelihoods[i] = (int) Math.round(normalizedLikelihoodsUnrounded[i]);
            }

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
                    final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a2);
                    idxAA = idxVector[0];
                    idxAB = idxVector[1];
                    idxBB = idxVector[2];
                }
                //I expect hets to be reference first, but there are no guarantees (e.g. phasing)
                else if (a1.isNonReference()) {
                    final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a1);
                    idxAA = idxVector[0];
                    idxAB = idxVector[1];
                    idxBB = idxVector[2];
                }
            }

            refCount += normalizedLikelihoods[idxAA];
            hetCount += normalizedLikelihoods[idxAB];
            homCount += normalizedLikelihoods[idxBB];
        }

        /*
         * Note: all that likelihood normalization etc may have accumulated some error.
         * We smooth it out my rounding the numbers to integers before the final computation.
         */
        refCount = Math.round(refCount);
        hetCount = Math.round(hetCount);
        homCount = Math.round(homCount);
        return new GenotypeCounts(refCount, hetCount, homCount);
    }
}
