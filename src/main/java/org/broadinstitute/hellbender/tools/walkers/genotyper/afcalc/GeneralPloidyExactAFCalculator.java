package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

public final class GeneralPloidyExactAFCalculator extends ExactAFCalculator {

    private static final boolean VERBOSE = false;
    private static final GenotypeLikelihoodCalculators GL_CALCS = new GenotypeLikelihoodCalculators();

    @Override
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy, final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(log10AlleleFrequencyPriors, "log10AlleleFrequencyPriors is null");
        Utils.nonNull(stateTracker, "stateTracker is null");
        combineSinglePools(vc.getGenotypes(), defaultPloidy, vc.getNAlleles(), log10AlleleFrequencyPriors);
        return getResultFromFinalState(vc, log10AlleleFrequencyPriors, stateTracker);
    }

    /**
     * Simple wrapper class to hold values of combined pool likelihoods.
     * For fast hashing and fast retrieval, there's a hash map that shadows main list.
     *
     */
    private static final class CombinedPoolLikelihoods {
        private final List<ExactACset> alleleCountSetList;
        private final Map<ExactACcounts,ExactACset> conformationMap;
        private double maxLikelihood;

        CombinedPoolLikelihoods() {
            // final int numElements = GenotypeLikelihoods.numLikelihoods();
            alleleCountSetList = new LinkedList<>();
            conformationMap = new LinkedHashMap<>();
            maxLikelihood = Double.NEGATIVE_INFINITY;
        }

        public void add(final ExactACset set) {
            alleleCountSetList.add(set);
            conformationMap.put(set.getACcounts(), set);
            final double likelihood = set.getLog10Likelihoods()[0];

            if (likelihood > maxLikelihood ) {
                maxLikelihood = likelihood;
            }
        }

        public boolean hasConformation(final int[] ac) {
            return conformationMap.containsKey(new ExactACcounts(ac));
        }

        public double getLikelihoodOfConformation(final int[] ac) {
            return conformationMap.get(new ExactACcounts(ac)).getLog10Likelihoods()[0];
        }

        public double getGLOfACZero() {
            return alleleCountSetList.get(0).getLog10Likelihoods()[0]; // AC 0 is always at beginning of list
        }

        public int getLength() {
            return alleleCountSetList.size();
        }
    }



    /**
     * Simple non-optimized version that combines GLs from several pools and produces global AF distribution.
     * @param GLs                              Inputs genotypes context with per-pool GLs
     * @param numAlleles                       Number of alternate alleles
     * @param log10AlleleFrequencyPriors       Frequency priors
     */
    @VisibleForTesting
    void combineSinglePools(final GenotypesContext GLs,
                                    final int defaultPloidy,
                                    final int numAlleles,
                                    final double[] log10AlleleFrequencyPriors) {

        // Combine each pool incrementally - likelihoods will be renormalized at each step

        // first element: zero ploidy, e.g. trivial degenerate distribution
        final int numAltAlleles = numAlleles - 1;
        final int[] zeroCounts = new int[numAlleles];
        final ExactACset set = new ExactACset(1, new ExactACcounts(zeroCounts));
        set.getLog10Likelihoods()[0] = 0.0;
        final StateTracker stateTracker = getStateTracker(false,numAltAlleles);
        int combinedPloidy = 0;
        CombinedPoolLikelihoods combinedPoolLikelihoods = new CombinedPoolLikelihoods();
        combinedPoolLikelihoods.add(set);

        for (final Genotype genotype : GLs.iterateInSampleNameOrder()) {
            // recover gls and check if they qualify.
            if (!genotype.hasPL()) {
                continue;
            }
            final double[] gls = genotype.getLikelihoods().getAsVector();
            if (!GATKVariantContextUtils.isInformative(gls)) {
                continue;
            }
            stateTracker.reset();
            final int declaredPloidy = genotype.getPloidy();
            final int ploidy = declaredPloidy < 1 ? defaultPloidy : declaredPloidy;
            // they do qualify so we proceed.
            combinedPoolLikelihoods = fastCombineMultiallelicPool(combinedPoolLikelihoods, gls,
                    combinedPloidy, ploidy, numAlleles, log10AlleleFrequencyPriors, stateTracker);
            combinedPloidy = ploidy + combinedPloidy; // total number of chromosomes in combinedLikelihoods
        }
        if (combinedPloidy == 0) {
            stateTracker.setLog10LikelihoodOfAFzero(0.0);
        }
    }

    private CombinedPoolLikelihoods fastCombineMultiallelicPool(final CombinedPoolLikelihoods originalPool,
                                                               final double[] newGL,
                                                               final int originalPloidy,
                                                               final int newGLPloidy,
                                                               final int numAlleles,
                                                               final double[] log10AlleleFrequencyPriors,
                                                               final StateTracker stateTracker) {
        final Deque<ExactACset> ACqueue = new LinkedList<>();
        // mapping of ExactACset indexes to the objects
        final Map<ExactACcounts, ExactACset> indexesToACset = new LinkedHashMap<>();
        final CombinedPoolLikelihoods newPool = new CombinedPoolLikelihoods();

        // add AC=0 to the queue
        final int[] zeroCounts = new int[numAlleles];
        final int newPloidy = originalPloidy + newGLPloidy;
        zeroCounts[0] = newPloidy;

        final ExactACset zeroSet = new ExactACset(1, new ExactACcounts(zeroCounts));

        ACqueue.add(zeroSet);
        indexesToACset.put(zeroSet.getACcounts(), zeroSet);

        // keep processing while we have AC conformations that need to be calculated
        while ( !ACqueue.isEmpty() ) {
            // compute log10Likelihoods
            final ExactACset ACset = ACqueue.remove();

            calculateACConformationAndUpdateQueue(ACset, newPool, originalPool, newGL, log10AlleleFrequencyPriors, originalPloidy, newGLPloidy, ACqueue, indexesToACset, stateTracker);

            // clean up memory
            indexesToACset.remove(ACset.getACcounts());
            if ( VERBOSE ) {
                System.out.printf(" *** removing used set=%s%n", ACset.getACcounts());
            }

        }
        return newPool;
    }

    /**
     *
     * @param set                       ExactACset holding conformation to be computed
     * @param newPool                   New pool likelihood holder
     * @param originalPool              Original likelihood holder
     * @param newGL                     New pool GL vector to combine
     * @param log10AlleleFrequencyPriors Prior object
     * @param originalPloidy             Total ploidy of original combined pool
     * @param newGLPloidy                Ploidy of GL vector
     * @param ACqueue                    Queue of conformations to compute
     * @param indexesToACset             AC indices of objects in queue
     * @return                           max log likelihood
     */
    private double calculateACConformationAndUpdateQueue(final ExactACset set,
                                                         final CombinedPoolLikelihoods newPool,
                                                         final CombinedPoolLikelihoods originalPool,
                                                         final double[] newGL,
                                                         final double[] log10AlleleFrequencyPriors,
                                                         final int originalPloidy,
                                                         final int newGLPloidy,
                                                         final Deque<ExactACset> ACqueue,
                                                         final Map<ExactACcounts, ExactACset> indexesToACset,
                                                         final StateTracker stateTracker) {

        // compute likelihood in "set" of new set based on original likelihoods
        final int numAlleles = set.getACcounts().getCounts().length;
        final int newPloidy = set.getACsum();
        final double log10LofK = computeLofK(set, originalPool, newGL, log10AlleleFrequencyPriors, numAlleles, originalPloidy, newGLPloidy, stateTracker);


        // add to new pool
        if (!Double.isInfinite(log10LofK)) {
            newPool.add(set);
        }

        if ( stateTracker.abort(log10LofK, set.getACcounts(), true, true) ) {
            return log10LofK;
        }

        // iterate over higher frequencies if possible
        // by convention, ACcounts contained in set have full vector of possible pool ac counts including ref count.
        // so, if first element is zero, it automatically means we have no wiggle since we're in a corner of the conformation space
        final int ACwiggle = set.getACcounts().getCounts()[0];
        if ( ACwiggle == 0 ) // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
        {
            return log10LofK;
        }


        // add conformations for other cases
        for ( int allele = 1; allele < numAlleles; allele++ ) {
            final int[] ACcountsClone = set.getACcounts().getCounts().clone();
            ACcountsClone[allele]++;
            // is this a valid conformation?
            final int altSum = (int)MathUtils.sum(ACcountsClone) - ACcountsClone[0];
            ACcountsClone[0] = newPloidy - altSum;
            if (ACcountsClone[0] < 0) {
                continue;
            }

            updateACset(ACcountsClone, ACqueue, indexesToACset);
        }


        return log10LofK;
    }

    /**
     * Compute likelihood of a particular AC conformation and update AFresult object
     * @param set                     Set of AC counts to compute
     * @param firstGLs                  Original pool likelihoods before combining
     * @param secondGL                  New GL vector with additional pool
     * @param log10AlleleFrequencyPriors     Allele frequency priors
     * @param numAlleles                Number of alleles (including ref)
     * @param ploidy1                   Ploidy of original pool (combined)
     * @param ploidy2                   Ploidy of new pool
     * @return                          log-likelihood of requested conformation
     */
    private double computeLofK(final ExactACset set,
                               final CombinedPoolLikelihoods firstGLs,
                               final double[] secondGL,
                               final double[] log10AlleleFrequencyPriors,
                               final int numAlleles, final int ploidy1, final int ploidy2, final StateTracker stateTracker) {

        final int newPloidy = ploidy1 + ploidy2;

        // sanity check
        int totalAltK = set.getACsum();
        if (newPloidy != totalAltK) {
            throw new GATKException("BUG: inconsistent sizes of set.getACsum and passed ploidy values");
        }

        totalAltK -= set.getACcounts().getCounts()[0];
        // totalAltK has sum of alt alleles of conformation now


        // special case for k = 0 over all k
        if ( totalAltK == 0 ) {   // all-ref case
            final double log10Lof0 = firstGLs.getGLOfACZero() + secondGL[HOM_REF_INDEX];
            set.getLog10Likelihoods()[0] = log10Lof0;
            stateTracker.setLog10LikelihoodOfAFzero(log10Lof0);
            stateTracker.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);
            return log10Lof0;

        }   else {
            // initialize result with denominator
            // ExactACset holds by convention the conformation of all alleles, and the sum of all allele count is just the ploidy.
            // To compute n!/k1!k2!k3!... we need to compute first n!/(k2!k3!...) and then further divide by k1! where k1=ploidy-sum_k_i

            final int[] currentCount = set.getACcounts().getCounts();
            final double denom =  -MathUtils.log10MultinomialCoefficient(newPloidy, currentCount);

            final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(ploidy2, numAlleles);
            for (int PLIndex = 0; PLIndex < glCalc.genotypeCount(); PLIndex++) {
                final GenotypeAlleleCounts alleleCounts = glCalc.genotypeAlleleCountsAt(PLIndex);

                final int[] acCount2 = alleleCounts.alleleCountsByIndex(numAlleles - 1);
                final int[] acCount1 = MathUtils.vectorDiff(currentCount, acCount2);
                // see if conformation is valid and if original pool had this conformation
                // for conformation to be valid, all elements of g2 have to be <= elements of current AC set
                if (isValidConformation(acCount1,ploidy1) && firstGLs.hasConformation(acCount1)) {
                    final double gl2 = secondGL[PLIndex];
                    if (!Double.isInfinite(gl2)) {
                        final double firstGL = firstGLs.getLikelihoodOfConformation(acCount1);
                        final double num1 = MathUtils.log10MultinomialCoefficient(ploidy1, acCount1);
                        final double num2 = MathUtils.log10MultinomialCoefficient(ploidy2, acCount2);
                        final double sum = firstGL + gl2 + num1 + num2;
                        set.getLog10Likelihoods()[0] = MathUtils.approximateLog10SumLog10(set.getLog10Likelihoods()[0], sum);
                    }
                }
            }
            set.getLog10Likelihoods()[0] += denom;
        }

        double log10LofK = set.getLog10Likelihoods()[0];

        // update the MLE if necessary
        final int[] altCounts = Arrays.copyOfRange(set.getACcounts().getCounts(), 1, set.getACcounts().getCounts().length);
        // TODO -- GUILLERMO THIS CODE MAY PRODUCE POSITIVE LIKELIHOODS OR -INFINITY
        stateTracker.updateMLEifNeeded(Math.max(log10LofK, -Double.MAX_VALUE), altCounts);

        // apply the priors over each alternate allele
        for (final int ACcount : altCounts ) {
            if ( ACcount > 0 ) {
                log10LofK += log10AlleleFrequencyPriors[ACcount];
            }
        }
        // TODO -- GUILLERMO THIS CODE MAY PRODUCE POSITIVE LIKELIHOODS OR -INFINITY
        stateTracker.updateMAPifNeeded(Math.max(log10LofK, -Double.MAX_VALUE), altCounts);

        return log10LofK;
    }

    //Small helper routine - is a particular AC conformation vector valid? ie are all elements non-negative and sum to ploidy?
    private static boolean isValidConformation(final int[] alleleCounts, final int ploidy) {
        return (Arrays.stream(alleleCounts).noneMatch(ac -> ac < 0) && MathUtils.sum(alleleCounts) == ploidy);
    }


    private static void updateACset(final int[] newSetCounts,
                                   final Deque<ExactACset> ACqueue,
                                   final Map<ExactACcounts, ExactACset> indexesToACset) {

        final ExactACcounts index = new ExactACcounts(newSetCounts);
        if ( !indexesToACset.containsKey(index) ) {
            final ExactACset newSet = new ExactACset(1, index);
            indexesToACset.put(index, newSet);
            ACqueue.add(newSet);
        }

    }
}
