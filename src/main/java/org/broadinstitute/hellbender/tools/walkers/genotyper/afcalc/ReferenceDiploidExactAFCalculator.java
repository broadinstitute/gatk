package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Reference implementation of multi-allelic EXACT model.  Extremely slow for many alternate alleles.
 */
public final class ReferenceDiploidExactAFCalculator extends ExactAFCalculator {

    private static final double LOG10_OF_2 = MathUtils.log10(2);

    protected ReferenceDiploidExactAFCalculator() {
    }

    @Override
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                                      final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(log10AlleleFrequencyPriors, "log10AlleleFrequencyPriors is null");
        Utils.nonNull(stateTracker, "stateTracker is null");
        final int numAlternateAlleles = vc.getNAlleles() - 1;

        final List<double[]> genotypeLikelihoods = getGLs(vc.getGenotypes(), true, vc.hasAllele(Allele.NON_REF_ALLELE));
        final int numSamples = genotypeLikelihoods.size()-1;
        final int numChr = 2*numSamples;

        // queue of AC conformations to process
        final Deque<ExactACset> ACqueue = new LinkedList<>();

        // mapping of ExactACset indexes to the objects
        final Map<ExactACcounts, ExactACset> indexesToACset = new LinkedHashMap<>(numChr+1);

        // add AC=0 to the queue
        final int[] zeroCounts = new int[numAlternateAlleles];
        final ExactACset zeroSet = new ExactACset(numSamples+1, new ExactACcounts(zeroCounts));
        ACqueue.add(zeroSet);
        indexesToACset.put(zeroSet.getACcounts(), zeroSet);

        while ( !ACqueue.isEmpty() ) {

            // compute log10Likelihoods
            final ExactACset set = ACqueue.remove();

            calculateAlleleCountConformation(set, genotypeLikelihoods, numChr, ACqueue, indexesToACset, log10AlleleFrequencyPriors,stateTracker);

            // clean up memory
            indexesToACset.remove(set.getACcounts());
        }

        return getResultFromFinalState(vc, log10AlleleFrequencyPriors, stateTracker);
    }


    private double calculateAlleleCountConformation(final ExactACset set,
                                                    final List<double[]> genotypeLikelihoods,
                                                    final int numChr,
                                                    final Deque<ExactACset> ACqueue,
                                                    final Map<ExactACcounts, ExactACset> indexesToACset,
                                                    final double[] log10AlleleFrequencyPriors,
                                                    final StateTracker stateTracker) {

        // compute the log10Likelihoods
        computeLofK(set, genotypeLikelihoods, log10AlleleFrequencyPriors, stateTracker);

        final double log10LofK = set.getLog10Likelihoods()[set.getLog10Likelihoods().length-1];

        // can we abort early because the log10Likelihoods are so small?
        if ( stateTracker.abort(log10LofK, set.getACcounts(), true, false) ) {
            return log10LofK;
        }

        // iterate over higher frequencies if possible
        final int ACwiggle = numChr - set.getACsum();
        if ( ACwiggle == 0 ){ // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
            return log10LofK;
        }

        final int numAltAlleles = set.getACcounts().getCounts().length;

        // add conformations for the k+1 case
        for ( int allele = 0; allele < numAltAlleles; allele++ ) {
            final int[] ACcountsClone = set.getACcounts().getCounts().clone();
            ACcountsClone[allele]++;
            // to get to this conformation, a sample would need to be AB (remember that ref=0)
            final int PLindex = GenotypeLikelihoods.calculatePLindex(0, allele + 1);
            updateACset(ACcountsClone, numChr, set, PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
        }

        // add conformations for the k+2 case if it makes sense; note that the 2 new alleles may be the same or different
        if ( ACwiggle > 1 ) {
            final List<DependentSet> differentAlleles = new ArrayList<>(numAltAlleles * numAltAlleles);
            final List<DependentSet> sameAlleles = new ArrayList<>(numAltAlleles);

            for ( int allele_i = 0; allele_i < numAltAlleles; allele_i++ ) {
                for ( int allele_j = allele_i; allele_j < numAltAlleles; allele_j++ ) {
                    final int[] ACcountsClone = set.getACcounts().getCounts().clone();
                    ACcountsClone[allele_i]++;
                    ACcountsClone[allele_j]++;

                    // to get to this conformation, a sample would need to be BB or BC (remember that ref=0, so add one to the index)
                    final int PLindex = GenotypeLikelihoods.calculatePLindex(allele_i + 1, allele_j + 1);
                    if ( allele_i == allele_j ) {
                        sameAlleles.add(new DependentSet(ACcountsClone, PLindex));
                    } else {
                        differentAlleles.add(new DependentSet(ACcountsClone, PLindex));
                    }
                }
            }

            // IMPORTANT: we must first add the cases where the 2 new alleles are different so that the queue maintains its ordering
            for ( final DependentSet dependent : differentAlleles ) {
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
            }
            for ( final DependentSet dependent : sameAlleles ) {
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
            }
        }

        return log10LofK;
    }

    private static void computeLofK(final ExactACset set,
                                    final List<double[]> genotypeLikelihoods,
                                    final double[] log10AlleleFrequencyPriors,
                                    final StateTracker stateTracker) {

        final double[] setLog10Likelihoods = set.getLog10Likelihoods();
        setLog10Likelihoods[0] = 0.0; // the zero case
        final int totalK = set.getACsum();

        // special case for k = 0 over all k
        if ( totalK == 0 ) {
            for (int j = 1, n = setLog10Likelihoods.length; j < n; j++ ) {
                setLog10Likelihoods[j] = setLog10Likelihoods[j - 1] + genotypeLikelihoods.get(j)[HOM_REF_INDEX];
            }

            final double log10Lof0 = setLog10Likelihoods[setLog10Likelihoods.length-1];
            stateTracker.setLog10LikelihoodOfAFzero(log10Lof0);
            stateTracker.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);
            return;
        }

        // if we got here, then k > 0 for at least one k.
        // the non-AA possible conformations were already dealt with by pushes from dependent sets;
        // now deal with the AA case (which depends on previous cells in this column) and then update the L(j,k) value
        for (int j = 1, n = setLog10Likelihoods.length; j < n; j++ ) {
            if ( totalK < 2*j-1 ) {
                final double[] gl = genotypeLikelihoods.get(j);
                final double conformationValue = MathUtils.log10(2*j-totalK) + MathUtils.log10(2*j-totalK-1) + setLog10Likelihoods[j-1] + gl[HOM_REF_INDEX];
                setLog10Likelihoods[j] = MathUtils.approximateLog10SumLog10(setLog10Likelihoods[j], conformationValue);
            }

            final double logDenominator = MathUtils.log10(2*j) + MathUtils.log10(2*j-1);
            setLog10Likelihoods[j] = setLog10Likelihoods[j] - logDenominator;
        }

        double log10LofK = setLog10Likelihoods[setLog10Likelihoods.length-1];

        // update the MLE if necessary
        stateTracker.updateMLEifNeeded(log10LofK, set.getACcounts().getCounts());

        // apply the priors over each alternate allele
        for ( final int ACcount : set.getACcounts().getCounts() ) {
            if ( ACcount > 0 ) {
                log10LofK += log10AlleleFrequencyPriors[ACcount];
            }
        }

        stateTracker.updateMAPifNeeded(log10LofK, set.getACcounts().getCounts());
    }

    private static final class DependentSet {
        public final int[] ACcounts;
        public final int PLindex;

        DependentSet(final int[] ACcounts, final int PLindex) {
            this.ACcounts = ACcounts;
            this.PLindex = PLindex;
        }
    }


    // adds the ExactACset represented by the ACcounts to the ACqueue if not already there (creating it if needed) and
    // also pushes its value to the given callingSetIndex.
    private static void updateACset(final int[] newSetCounts,
                                    final int numChr,
                                    final ExactACset dependentSet,
                                    final int PLsetIndex,
                                    final Queue<ExactACset> ACqueue,
                                    final Map<ExactACcounts, ExactACset> indexesToACset,
                                    final List<double[]> genotypeLikelihoods) {
        final ExactACcounts index = new ExactACcounts(newSetCounts);
        if ( !indexesToACset.containsKey(index) ) {
            final ExactACset set = new ExactACset(numChr/2 +1, index);
            indexesToACset.put(index, set);
            ACqueue.add(set);
        }

        // push data from the dependency to the new set
        pushData(indexesToACset.get(index), dependentSet, PLsetIndex, genotypeLikelihoods);
    }

    private static void pushData(final ExactACset targetSet,
                                 final ExactACset dependentSet,
                                 final int PLsetIndex,
                                 final List<double[]> genotypeLikelihoods) {
        final int totalK = targetSet.getACsum();

        final double[] targetSetLog10Likelihoods = targetSet.getLog10Likelihoods();
        final double[] dependentSetLog10Likelihoods = dependentSet.getLog10Likelihoods();
        final int[] counts = targetSet.getACcounts().getCounts();

        for ( int j = 1, n = targetSetLog10Likelihoods.length; j < n; j++ ) {
            if (2 * j >= totalK) { // skip impossible conformations
                final double[] gl = genotypeLikelihoods.get(j);
                final double conformationValue =
                        determineCoefficient(PLsetIndex, j, counts, totalK) + dependentSetLog10Likelihoods[j-1] + gl[PLsetIndex];
                targetSetLog10Likelihoods[j] = MathUtils.approximateLog10SumLog10(targetSetLog10Likelihoods[j], conformationValue);
            }
        }
    }

    private static double determineCoefficient(final int PLindex, final int j, final int[] ACcounts, final int totalK) {
        // the closed form representation generalized for multiple alleles is as follows:
        // AA: (2j - totalK) * (2j - totalK - 1)
        // AB: 2k_b * (2j - totalK)
        // AC: 2k_c * (2j - totalK)
        // BB: k_b * (k_b - 1)
        // BC: 2 * k_b * k_c
        // CC: k_c * (k_c - 1)

        // find the 2 alleles that are represented by this PL index
        final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindex);

        // *** note that throughout this method we subtract one from the alleleIndex because ACcounts ***
        // *** doesn't consider the reference allele whereas the GenotypeLikelihoods PL cache does.   ***

        // the AX het case
        if ( alleles.alleleIndex1 == 0 ) {
            return MathUtils.log10(2 * ACcounts[alleles.alleleIndex2 - 1]) + MathUtils.log10(2 * j - totalK);
        }

        final int k_i = ACcounts[alleles.alleleIndex1-1];

        // the hom var case (e.g. BB, CC, DD)
        final double coeff;
        if ( alleles.alleleIndex1 == alleles.alleleIndex2 ) {
            coeff = MathUtils.log10(k_i) + MathUtils.log10(k_i - 1);
        } else {        // the het non-ref case (e.g. BC, BD, CD)
            final int k_j = ACcounts[alleles.alleleIndex2-1];
            coeff = LOG10_OF_2 + MathUtils.log10(k_i) + MathUtils.log10(k_j);
        }

        return coeff;
    }
}
