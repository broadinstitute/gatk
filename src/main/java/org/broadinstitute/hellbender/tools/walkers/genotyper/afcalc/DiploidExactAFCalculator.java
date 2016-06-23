package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

public abstract class DiploidExactAFCalculator extends ExactAFCalculator {

    protected abstract AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                                      final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker);

    @Override
    protected final GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(allelesToUse, "allelesToUse is null");
        return GATKVariantContextUtils.subsetAlleles(vc, allelesToUse, GATKVariantContextUtils.GenotypeAssignmentMethod.SET_TO_NO_CALL);
    }

    @Override
    protected final void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(likelihoodSums, "likelihoodSums is null");
        final List<double[]> GLs = getGLs(vc.getGenotypes(), true);
        for ( final double[] likelihoods : GLs ) {
            final int PLindexOfBestGL = MathUtils.maxElementIndex(likelihoods);
            if ( PLindexOfBestGL != PL_INDEX_OF_HOM_REF ) {
                final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindexOfBestGL);
                final int alleleLikelihoodIndex1 = alleles.alleleIndex1 - 1;
                final int alleleLikelihoodIndex2 = alleles.alleleIndex2 - 1;
                if ( alleles.alleleIndex1 != 0 ) {
                    likelihoodSums[alleleLikelihoodIndex1].sum += likelihoods[PLindexOfBestGL] - likelihoods[PL_INDEX_OF_HOM_REF];
                }
                // don't double-count it
                if ( alleles.alleleIndex2 != 0 && alleles.alleleIndex2 != alleles.alleleIndex1 ) {
                    likelihoodSums[alleleLikelihoodIndex2].sum += likelihoods[PLindexOfBestGL] - likelihoods[PL_INDEX_OF_HOM_REF];
                }
            }
        }
    }

    @Override
    public final GenotypesContext subsetAlleles(final VariantContext vc,
                                          final int defaultPloidy,
                                          final List<Allele> allelesToUse,
                                          final boolean assignGenotypes) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(allelesToUse, "allelesToUse is null");

        return allelesToUse.size() == 1
                ? GATKVariantContextUtils.subsetToRefOnly(vc, defaultPloidy)
                : GATKVariantContextUtils.subsetAlleles(vc, allelesToUse,
                     assignGenotypes ? GATKVariantContextUtils.GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN : GATKVariantContextUtils.GenotypeAssignmentMethod.SET_TO_NO_CALL);
    }
}
