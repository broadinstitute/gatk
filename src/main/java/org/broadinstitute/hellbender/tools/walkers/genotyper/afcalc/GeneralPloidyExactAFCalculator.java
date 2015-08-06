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

    static final int MAX_LENGTH_FOR_POOL_PL_LOGGING = 100; // if PL vectors longer than this # of elements, don't log them

    private static final boolean VERBOSE = false;

    @Override
    protected GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(allelesToUse, "allelesToUse is null");
        return subsetAlleles(vc,defaultPloidy,allelesToUse,false);
    }

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
            conformationMap = new HashMap<>();
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

    @Override
    protected void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(likelihoodSums, "likelihoodSums is null");

        final int numOriginalAltAlleles = likelihoodSums.length;
        final GenotypesContext genotypes = vc.getGenotypes();
        for ( final Genotype genotype : genotypes.iterateInSampleNameOrder() ) {
            if (!genotype.hasPL()) {
                continue;
            }
            final double[] gls = genotype.getLikelihoods().getAsVector();
            if (MathUtils.sum(gls) >= GATKVariantContextUtils.SUM_GL_THRESH_NOCALL) {
                continue;
            }

            final int PLindexOfBestGL = MathUtils.maxElementIndex(gls);

            final double bestToHomRefDiffGL = PLindexOfBestGL == PL_INDEX_OF_HOM_REF ? 0.0 : gls[PLindexOfBestGL] - gls[PL_INDEX_OF_HOM_REF];
            final int declaredPloidy = genotype.getPloidy();
            final int ploidy = declaredPloidy <= 0 ? defaultPloidy : declaredPloidy;

            final int[] acCount = getAlleleCountFromPLIndex(1 + numOriginalAltAlleles, ploidy, PLindexOfBestGL);
            // by convention, first count coming from getAlleleCountFromPLIndex comes from reference allele
            for (int k=1; k < acCount.length;k++) {
                if (acCount[k] > 0) {
                    likelihoodSums[k - 1].sum += acCount[k] * bestToHomRefDiffGL;
                }
            }
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
            if (MathUtils.sum(gls) >= GATKVariantContextUtils.SUM_GL_THRESH_NOCALL) {
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
        final Map<ExactACcounts, ExactACset> indexesToACset = new HashMap<>();
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

            // for current conformation, get all possible ways to break vector K into two components G1 and G2
            final SumIterator innerIterator = new SumIterator(numAlleles,ploidy2);
            set.getLog10Likelihoods()[0] = Double.NEGATIVE_INFINITY;
            while (innerIterator.hasNext()) {
                // check if breaking current conformation into g1 and g2 is feasible.
                final int[] acCount2 = innerIterator.getCurrentVector();
                final int[] acCount1 = MathUtils.vectorDiff(currentCount, acCount2);
                final int idx2 = innerIterator.getLinearIndex();
                // see if conformation is valid and if original pool had this conformation
                // for conformation to be valid, all elements of g2 have to be <= elements of current AC set
                if (isValidConformation(acCount1,ploidy1) && firstGLs.hasConformation(acCount1)) {
                    final double gl2 = secondGL[idx2];
                    if (!Double.isInfinite(gl2)) {
                        final double firstGL = firstGLs.getLikelihoodOfConformation(acCount1);
                        final double num1 = MathUtils.log10MultinomialCoefficient(ploidy1, acCount1);
                        final double num2 = MathUtils.log10MultinomialCoefficient(ploidy2, acCount2);
                        final double sum = firstGL + gl2 + num1 + num2;

                        set.getLog10Likelihoods()[0] = MathUtils.approximateLog10SumLog10(set.getLog10Likelihoods()[0], sum);
                    }
                }
                innerIterator.next();
            }

            set.getLog10Likelihoods()[0] += denom;
        }

        double log10LofK = set.getLog10Likelihoods()[0];

        // update the MLE if necessary
        final int altCounts[] = Arrays.copyOfRange(set.getACcounts().getCounts(), 1, set.getACcounts().getCounts().length);
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

    /**
     * Small helper routine - is a particular AC conformation vector valid? ie are all elements non-negative and sum to ploidy?
     * @param set                            AC conformation vector
     * @param ploidy                         Ploidy of set
     * @return                               Valid conformation
     */
    private static boolean isValidConformation(final int[] set, final int ploidy) {
        int sum=0;
        for (final int ac: set) {
            if (ac < 0) {
                return false;
            }
            sum += ac;
        }

        return (sum == ploidy);
    }

    /**
     * From a given variant context, extract a given subset of alleles, and update genotype context accordingly,
     * including updating the PL's, and assign genotypes accordingly
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param defaultPloidy                     ploidy to assume in case that {@code vc} does not contain that information
     *                                          for a sample.
     * @param allelesToUse                      alleles to subset
     * @param assignGenotypes                   true: assign hard genotypes, false: leave as no-call
     * @return                                  GenotypesContext with new PLs
     */
    public GenotypesContext subsetAlleles(final VariantContext vc,
                                          final int defaultPloidy,
                                          final List<Allele> allelesToUse,
                                          final boolean assignGenotypes) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(allelesToUse, "allelesToUse is null");

        // the genotypes with PLs
        final GenotypesContext oldGTs = vc.getGenotypes();

        // samples
        final List<String> sampleIndices = oldGTs.getSampleNamesOrderedByName();

        // the new genotypes to create
        final GenotypesContext newGTs = GenotypesContext.create();

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;


        // create the new genotypes
        for ( int k = 0; k < oldGTs.size(); k++ ) {
            final Genotype g = oldGTs.get(sampleIndices.get(k));
            final int declaredPloidy = g.getPloidy();
            final int ploidy = declaredPloidy <= 0 ? defaultPloidy : declaredPloidy;
            if ( !g.hasLikelihoods() ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), GATKVariantContextUtils.noCallAlleles(ploidy)));
                continue;
            }

            // create the new likelihoods array from the alleles we are allowed to use
            final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
            double[] newLikelihoods;

            // Optimization: if # of new alt alleles = 0 (pure ref call), keep original likelihoods so we skip normalization
            // and subsetting
            if ( numOriginalAltAlleles == numNewAltAlleles || numNewAltAlleles == 0) {
                newLikelihoods = originalLikelihoods;
            } else {
                newLikelihoods = subsetToAlleles(originalLikelihoods, ploidy, vc.getAlleles(), allelesToUse);

                // might need to re-normalize
                newLikelihoods = MathUtils.normalizeFromLog10(newLikelihoods, false, true);
            }

            // if there is no mass on the (new) likelihoods, then just no-call the sample
            if ( MathUtils.sum(newLikelihoods) > GATKVariantContextUtils.SUM_GL_THRESH_NOCALL ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), GATKVariantContextUtils.noCallAlleles(ploidy)));
            }
            else {
                final GenotypeBuilder gb = new GenotypeBuilder(g);

                if ( numNewAltAlleles == 0 ) {
                    gb.noPL();
                } else {
                    gb.PL(newLikelihoods);
                }

                // if we weren't asked to assign a genotype, then just no-call the sample
                if ( !assignGenotypes || MathUtils.sum(newLikelihoods) > GATKVariantContextUtils.SUM_GL_THRESH_NOCALL ) {
                    gb.alleles(GATKVariantContextUtils.noCallAlleles(ploidy));
                } else {
                    assignGenotype(gb, newLikelihoods, allelesToUse, ploidy);
                }
                newGTs.add(gb.make());
            }
        }

        return newGTs;

    }

    /**
     * Assign genotypes (GTs) to the samples in the Variant Context greedily based on the PLs
     *
     * @param newLikelihoods       the PL array
     * @param allelesToUse         the list of alleles to choose from (corresponding to the PLs)
     * @param numChromosomes        Number of chromosomes per pool
     */
    private static void assignGenotype(final GenotypeBuilder gb,
                                       final double[] newLikelihoods,
                                       final List<Allele> allelesToUse,
                                       final int numChromosomes) {
        final int numNewAltAlleles = allelesToUse.size() - 1;

        // find the genotype with maximum likelihoods
        final int PLindex = numNewAltAlleles == 0 ? 0 : MathUtils.maxElementIndex(newLikelihoods);
        final GenotypeLikelihoodCalculator calculator = new GenotypeLikelihoodCalculators().getInstance(numChromosomes, allelesToUse.size());
        final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(PLindex);

        gb.alleles(alleleCounts.asAlleleList(allelesToUse));

        // remove PLs if necessary
        if (newLikelihoods.length > MAX_LENGTH_FOR_POOL_PL_LOGGING) {
            gb.noPL();
        }

        if ( numNewAltAlleles > 0 ) {
            gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods));
        }
    }

    private static final int MAX_NUM_ALLELES_TO_CACHE = 20;
    private static final int MAX_NUM_SAMPLES_PER_POOL = 1000;


    /**
     * Crucial inner class that handles addressing elements of pool likelihoods. We store likelihoods as a map
     * of form int[] -> double (to be more precise, IntArrayWrapper -> Double).
     * For a given ploidy (chromosome count) and number of alleles, we need a form to iterate deterministically
     * across all possible allele conformations.
     * Problem equivalent to listing in determistic order all possible ways in which N integers will sum to P,
     * where N is number of alleles and P is number of chromosomes.
     * There's an option to list all integers so that sum will be UP to P.
     * For example, with P=2,N=2, restrictSumTo = 2 iterator will produce
     * [2 0 ] [1 1] [ 0 2]
     *
     *
     */
    private static final class SumIterator {
        private int[] currentState;
        private final int[] finalState;
        private final int restrictSumTo;
        private final int dim;
        private boolean hasNext;
        private int linearIndex;
        private int currentSum;

        /**
         * Default constructor. Typical use case: restrictSumTo = -1 if there's no sum restriction, or will generate int[]
         * vectors so that all add to this value.
         *
         * @param finalState                    End state - typically we should set value to (P,P,P,...)
         * @param restrictSumTo                 See above
         */
        public SumIterator(final int[] finalState,final int restrictSumTo) {
            this.finalState = finalState;
            this.dim = finalState.length;
            this.restrictSumTo = restrictSumTo;
            currentState = new int[dim];
            reset();

        }

        /**
         * Shortcut constructor for common use case: iterator will produce
         * all vectors of length numAlleles whose sum = numChromosomes
         * @param numAlleles              Number of alleles
         * @param numChromosomes          Ploidy
         */
        public SumIterator(final int numAlleles, final int numChromosomes) {
            this(getInitialStateVector(numAlleles, numChromosomes), numChromosomes);
        }


        private static int[] getInitialStateVector(final int nAlleles, final int numChromosomes) {
            final int[] initialState = new int[nAlleles];
            Arrays.fill(initialState, numChromosomes);
            return initialState;
        }

        public void next() {
            final int initialDim = (restrictSumTo > 0)?1:0;
            hasNext = next(finalState, initialDim);
            if (hasNext) {
                linearIndex++;
            }
        }

        private boolean next(final int[] finalState, final int initialDim) {
            boolean hasNextState = false;
            for (int currentDim=initialDim; currentDim < finalState.length; currentDim++) {
                final int x = currentState[currentDim]+1;

                if (x > finalState[currentDim] || (currentSum >= restrictSumTo && initialDim > 0)) {
                    // update vector sum, and reset position
                    currentSum -= currentState[currentDim];
                    currentState[currentDim] = 0;
                    if (currentDim >= dim-1) {
                        hasNextState = false;
                        break;
                    }
                }
                else {
                    currentState[currentDim] = x;
                    hasNextState = true;
                    currentSum++;
                    break;
                }
            }
            if (initialDim > 0) {
                currentState[0] = restrictSumTo - currentSum;
            }
            return hasNextState;
        }

        public void reset() {
            Arrays.fill(currentState, 0);
            if (restrictSumTo > 0) {
                currentState[0] = restrictSumTo;
            }
            hasNext = true;
            linearIndex = 0;
            currentSum = 0;
        }
        public int[] getCurrentVector() {
            return currentState;
        }

        public int getLinearIndex() {
            return linearIndex;
        }

        public boolean hasNext() {
            return hasNext;
        }
    }


    /**
     * Given set of alleles with corresponding vector of likelihoods, subset to a new set of alleles
     *
     * @param oldLikelihoods        Vector of PL's corresponding to original alleles
     * @param numChromosomes        Ploidy (number of chromosomes describing PL's)
     * @param originalAlleles       List of original alleles
     * @param allelesToSubset       Alleles to subset
     * @return                      Vector of new PL's, ordered accorrding to SumIterator's ordering
     */
    private static double[] subsetToAlleles(final double[] oldLikelihoods, final int numChromosomes,
                                           final List<Allele> originalAlleles, final List<Allele> allelesToSubset) {

        final int newPLSize = getNumLikelihoodElements(allelesToSubset.size(), numChromosomes);
        final double[] newPLs = new double[newPLSize];


        int idx = 0;
        // First fill boolean array stating whether each original allele is present in new mapping
        final boolean [] allelePresent = new boolean[originalAlleles.size()];
        for ( final Allele allele : originalAlleles ) {
            allelePresent[idx++] = allelesToSubset.contains(allele);
        }


        // compute mapping from old idx to new idx
        // This might be needed in case new allele set is not ordered in the same way as old set
        // Example. Original alleles: {T*,C,G,A}. New alleles: {G,C}. Permutation key = [2,1]

        final int[] permutationKey = new int[allelesToSubset.size()];
        for (int k=0; k < allelesToSubset.size(); k++)
        // for each allele to subset, find corresponding index in original allele list
        {
            permutationKey[k] = originalAlleles.indexOf(allelesToSubset.get(k));
        }


        final SumIterator iterator = new SumIterator(originalAlleles.size(),numChromosomes);

        while (iterator.hasNext()) {
            // for each entry in logPL table, associated originally with allele count stored in vec[],
            // see if this allele count conformation will be present in new logPL table.
            // For entry to be present, elements in dimensions not present in requested allele list have to have count = 0
            final int[] pVec = iterator.getCurrentVector();
            final double pl = oldLikelihoods[iterator.getLinearIndex()];

            boolean keyPresent = true;
            for (int k=0; k < allelePresent.length; k++) {
                if (pVec[k] > 0 && !allelePresent[k]) {
                    keyPresent = false;
                }
            }

            if (keyPresent) {// skip to next entry in logPLs if this conformation is not present in subset

                final int[] newCount = new int[allelesToSubset.size()];

                // map from old allele mapping count to new allele mapping
                // In pseudo-Matlab notation: newCount = vec[permutationKey] for permutationKey vector
                for (idx = 0; idx < newCount.length; idx++) {
                    newCount[idx] = pVec[permutationKey[idx]];
                }

                // get corresponding index from new count
                final int outputIdx = getLinearIndex(newCount, allelesToSubset.size(), numChromosomes);
                newPLs[outputIdx] = pl;
            }
            iterator.next();
        }

        return  newPLs;
    }

    private static int getLinearIndex(final int[] vectorIdx, final int numAlleles, final int ploidy) {

        if (ploidy <= 0) {
            return 0;
        }

        int linearIdx = 0;
        int cumSum = ploidy;
        for (int k=numAlleles-1;k>=1; k--) {
            final int idx = vectorIdx[k];
            // how many blocks are before current position
            if (idx == 0) {
                continue;
            }
            for (int p=0; p < idx; p++) {
                linearIdx += getNumLikelihoodElements(k, cumSum - p);
            }

            cumSum -= idx;
        }

        return linearIdx;

    }

    /**
     * Given a scalar index, what's the alelle count conformation corresponding to it?
     * @param nAlleles                    Number of alleles
     * @param numChromosomes              Ploidy
     * @param PLindex                     Index to query
     * @return                            Allele count conformation, according to iteration order from SumIterator
     */
    private static int[] getAlleleCountFromPLIndex(final int nAlleles, final int numChromosomes, final int PLindex) {

        final GenotypeLikelihoodCalculator calculator = new GenotypeLikelihoodCalculators().getInstance(numChromosomes, nAlleles);
        final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(PLindex);
        return alleleCounts.alleleCountsByIndex(nAlleles - 1);
    }

    /*
    * a cache of the PL ivector sizes as a function of # of alleles and pool sizes
    */
    @VisibleForTesting
    static int getNumLikelihoodElements(final int numAlleles, final int ploidy) {
        return GLVECTORSIZES[numAlleles][ploidy];
    }

    //Note: this is shared state but it's not modified at runtime
    private static final int[][] GLVECTORSIZES = fillGLVectorSizeCache(MAX_NUM_ALLELES_TO_CACHE, 2*MAX_NUM_SAMPLES_PER_POOL);

    private static int[][] fillGLVectorSizeCache(final int maxAlleles, final int maxPloidy) {

        final int[][] cache = new int[maxAlleles][maxPloidy];
        for (int numAlleles=1; numAlleles < maxAlleles; numAlleles++) {
            for (int ploidy=0; ploidy < maxPloidy; ploidy++) {
                if (numAlleles == 1) {
                    cache[numAlleles][ploidy] = 1;
                } else if (ploidy == 1) {
                    cache[numAlleles][ploidy] = numAlleles;
                } else {
                    int acc =0;
                    for (int k=0; k <= ploidy; k++ ) {
                        acc += cache[numAlleles - 1][ploidy - k];
                    }

                    cache[numAlleles][ploidy] = acc;
                }
            }
        }
        return cache;
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
