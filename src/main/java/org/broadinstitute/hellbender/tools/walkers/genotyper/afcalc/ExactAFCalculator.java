package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Uses the Exact calculation of Heng Li
 */
abstract class ExactAFCalculator extends AFCalculator {

    protected static final int HOM_REF_INDEX = 0;  // AA likelihoods are always first
    protected static final int PL_INDEX_OF_HOM_REF = 0;

    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with higher likelihood are first.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_SUM_COMPARATOR = Comparator.<LikelihoodSum>comparingDouble(o->o.sum).reversed();

    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with higher likelihood are first but make sure that
     * NON_REF alleles are last.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_NON_REF_THEN_SUM_COMPARATOR = (o1, o2) -> {
        if (o1.allele == GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) {
            return 1;
        } else if (o2.allele == GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) {
            return -1;
        } else {
            return o1.compareTo(o2);
        }
    };
    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with lower alternative allele index are first regardless of
     * the likelihood sum.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_INDEX_COMPARATOR = Comparator.comparingInt(o->o.index);

    /**
     * Wrapper class that compares two likelihoods associated with two alleles
     */
    protected static final class LikelihoodSum implements Comparable<LikelihoodSum> {
        public double sum = 0.0;
        public final Allele allele;
        public final int index;

        public LikelihoodSum(final Allele allele, final int index) { this.allele = allele; this.index = index; }

        public int compareTo(final LikelihoodSum other) {
            final double diff = Double.compare(sum, other.sum);
            return ( diff < 0.0 ) ? 1 : (diff > 0.0 ) ? -1 : 0;
        }
    }

    /**
     * Unpack GenotypesContext into arraylist of doubel values
     * @param GLs            Input genotype context
     * @return               ArrayList of doubles corresponding to GL vectors
     */
    protected static List<double[]> getGLs(final GenotypesContext GLs, final boolean includeDummy) {
        final List<double[]> genotypeLikelihoods = new ArrayList<>(GLs.size() + 1);

        if ( includeDummy ) {
            genotypeLikelihoods.add(new double[]{0.0, 0.0, 0.0}); // dummy
        }
        for ( final Genotype sample : GLs.iterateInSampleNameOrder() ) {
            if ( sample.hasLikelihoods() ) {
                final double[] gls = sample.getLikelihoods().getAsVector();

                if ( MathUtils.sum(gls) < GATKVariantContextUtils.SUM_GL_THRESH_NOCALL ) {
                    genotypeLikelihoods.add(gls);
                }
            }
        }

        return genotypeLikelihoods;
    }

    @Override
    protected VariantContext reduceScope(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles) {
        Utils.nonNull(vc, "vc is null");
        // don't try to genotype too many alternate alleles
        final List<Allele> inputAltAlleles = vc.getAlternateAlleles();
        final List<Allele> outputAltAlleles = reduceScopeAlleles(vc, defaultPloidy, maximumAlternativeAlleles);

        // only if output allele has reduced from the input alt allele set size we should care.
        final int altAlleleReduction = inputAltAlleles.size() - outputAltAlleles.size();

        if (altAlleleReduction == 0) {
            return vc;
        }
        logger.warn("this tool is currently set to genotype at most " + maximumAlternativeAlleles
                + " alternate alleles in a given context, but the context at " + vc.getContig() + ":" + vc.getStart()
                + " has " + (vc.getAlternateAlleles().size())
                + " alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument");

        final List<Allele> alleles = new ArrayList<>(maximumAlternativeAlleles + 1);
        alleles.add(vc.getReference());
        alleles.addAll(reduceScopeAlleles(vc, defaultPloidy, maximumAlternativeAlleles));
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        builder.alleles(alleles);
        builder.genotypes(reduceScopeGenotypes(vc, defaultPloidy, alleles));
        if (altAlleleReduction < 0) {
            throw new IllegalStateException("unexpected: reduction increased the number of alt. alleles!: " + -altAlleleReduction + " " + vc + " " + builder.make());
        }
        return builder.make();
    }

    /**
     * Returns a the new set of alleles to use.
     * @param vc target variant context.
     * @param numAllelesToChoose number of alleles to keep.
     * @return the list of alternative allele to keep.
     */
    protected List<Allele> reduceScopeAlleles(final VariantContext vc, final int defaultPloidy, final int numAllelesToChoose) {
        Utils.nonNull(vc, "vc is null");

        // Look  for the <NON_REF> allele to exclude it from the pruning if present.
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();

        final int nonRefAltAlleleIndex = GATKVariantContextUtils.indexOfAltAllele(vc, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, false);
        final boolean nonRefAltAllelePresent = nonRefAltAlleleIndex >= 0;

        // <NON_REF> should not be considered in the downsizing, so we need to count it out when
        // considering if alt. allele downsizing is required.
        final int numProperOriginalAltAlleles = numOriginalAltAlleles - (nonRefAltAllelePresent ? 1 : 0);

        // Avoid pointless allele reduction:
        if (numAllelesToChoose >= numProperOriginalAltAlleles) {
            return vc.getAlternateAlleles();
        }

        final LikelihoodSum[] likelihoodSums = new LikelihoodSum[numOriginalAltAlleles];
        for ( int i = 0; i < numOriginalAltAlleles; i++ ) {
            final Allele allele = vc.getAlternateAllele(i);
            likelihoodSums[i] = new LikelihoodSum(allele,i);
        }

        // Calculate the allele likelihood sums.
        reduceScopeCalculateLikelihoodSums(vc, defaultPloidy, likelihoodSums);

        // sort them by probability mass and choose the best ones
        // Make sure that the <NON_REF> allele is last if present.
        Collections.sort(Arrays.asList(likelihoodSums), nonRefAltAllelePresent ? LIKELIHOOD_NON_REF_THEN_SUM_COMPARATOR : LIKELIHOOD_SUM_COMPARATOR);

        // We need to return the best likelihood alleles in the original alternative allele index order.
        // This heap will keep track of that index order.
        final PriorityQueue<LikelihoodSum> mostLikelyAllelesHeapByIndex = new PriorityQueue<>(numOriginalAltAlleles, LIKELIHOOD_INDEX_COMPARATOR);

        for ( int i = 0; i < numAllelesToChoose; i++ ) {
            mostLikelyAllelesHeapByIndex.add(likelihoodSums[i]);
        }

        // guaranteed no to have been added at this point thanks for checking on whether reduction was
        // needed in the first place.
        if (nonRefAltAllelePresent) {
            mostLikelyAllelesHeapByIndex.add(likelihoodSums[nonRefAltAlleleIndex]);
        }

        final List<Allele> orderedBestAlleles = new ArrayList<>(numAllelesToChoose);

        while (!mostLikelyAllelesHeapByIndex.isEmpty()) {
            orderedBestAlleles.add(mostLikelyAllelesHeapByIndex.remove().allele);
        }

        return orderedBestAlleles;
    }

    /**
     * Update the likelihood sums with using the variant context genotype likelihoods.
     * @param vc source variant context.
     * @param likelihoodSums where to update the likelihood sums.
     */
    protected abstract void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums);

    /**
     * Transforms the genotypes of the variant context according to the new subset of possible alleles.
     *
     * @param vc original variant-context.
     * @param allelesToUse possible alleles.
     * @return never {@code null}, the new set of genotype calls for the reduced scope.
     */
    protected abstract GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse);
}