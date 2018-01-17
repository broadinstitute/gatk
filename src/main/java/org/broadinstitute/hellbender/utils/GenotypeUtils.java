package org.broadinstitute.hellbender.utils;

import htsjdk.variant.variantcontext.*;
import picard.util.MathUtil;

public final class GenotypeUtils {
    private GenotypeUtils(){}

    /**
     * Returns true if the genotype is a diploid genotype with likelihoods.
     */
    public static boolean isDiploidWithLikelihoods(final Genotype g) {
        return Utils.nonNull(g).hasLikelihoods() && g.getPloidy() == 2;
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

        int idxAA = 0;
        int idxAB = 1;
        int idxBB = 2;

        double refCount = 0;
        double hetCount = 0;
        double varCount = 0;

        for (final Genotype g : genotypes) {
            if (! isDiploidWithLikelihoods(g)){
                continue;
            }

            // Genotype::getLikelihoods returns a new array, so modification in-place is safe
            final double[] normalizedLikelihoods = MathUtils.normalizeFromLog10ToLinearSpace(g.getLikelihoods().getAsVector());
            final double[] biallelicLikelihoods;

            //if there are multiple alts, use the biallelic PLs for the best alt
            if (vc.getAlternateAlleles().size() > 1 ) {
                //check for
                int maxInd = MathUtil.indexOfMax(normalizedLikelihoods);
                GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(maxInd);
                if (alleles.alleleIndex1 != 0 && alleles.alleleIndex2 != 0) {
                    //all likelihoods go to varCount because no ref allele is called
                    varCount++;
                    continue;
                }

                double maxLikelihood = normalizedLikelihoods[idxAB];
                int hetIndex = idxAB;
                int varIndex = idxBB;
                for (final Allele currAlt : vc.getAlternateAlleles()) {
                    final int[] idxVector = vc.getGLIndicesOfAlternateAllele(currAlt);
                    int tempIndex = idxVector[1];
                    if (normalizedLikelihoods[tempIndex] > maxLikelihood) {
                        maxLikelihood = normalizedLikelihoods[tempIndex];
                        hetIndex = tempIndex;
                        varIndex = idxVector[2];
                    }
                }
                biallelicLikelihoods = MathUtils.normalizeFromRealSpace(new double[] {normalizedLikelihoods[idxAA], normalizedLikelihoods[hetIndex], normalizedLikelihoods[varIndex]});
            }
            else {
                biallelicLikelihoods = normalizedLikelihoods;
            }

            double refLikelihood = biallelicLikelihoods[idxAA];
            double hetLikelihood = biallelicLikelihoods[idxAB];
            double varLikelihood = biallelicLikelihoods[idxBB];

            //NOTE: rounding is special cased for [0,0,X] and [X,0,0] PLs because likelihoods can come out as [0.5, 0.5, 0] and both counts round up
            if( roundContributionFromEachGenotype ) {
                refCount += MathUtils.fastRound(refLikelihood);
                if (refLikelihood != hetLikelihood) {   //if GQ = 0 (specifically [0,0,X] PLs) count as homRef and don't add to the other counts
                    hetCount += MathUtils.fastRound(hetLikelihood);
                }
                if (varLikelihood != hetLikelihood) {  //if GQ = 0 (specifically [X,0,0] PLs) count as het and don't count as variant
                    varCount += MathUtils.fastRound(varLikelihood); //need specific varCount (rather than complement of the others) for PL[0,0,0] case
                }
            } else {
                refCount += refLikelihood;
                hetCount += hetLikelihood;
                varCount += varLikelihood;
            }
        }
        return new GenotypeCounts(refCount, hetCount, varCount);
    }
}
