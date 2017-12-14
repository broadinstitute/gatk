package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Computes the conditional bi-allelic exact results
 *
 * Suppose vc contains 2 alt allele: A* with C and T.  This function first computes:
 *
 * (1) P(D | AF_c > 0 && AF_t == *) [i.e., T can be anything]
 *
 * it then computes the conditional probability on AF_c == 0:
 *
 * (2) P(D | AF_t > 0 && AF_c == 0)
 *
 * Thinking about this visually, we have the following likelihood matrix where each cell is
 * the P(D | AF_c == i && AF_t == j):
 *
 *     0 AF_c > 0
 *    -----------------
 * 0  |  |
 *    |--|-------------
 * a  |  |
 * f  |  |
 * _  |  |
 * t  |  |
 * >  |  |
 * 0  |  |
 *
 * What we really want to know how
 *
 * (3) P(D | AF_c == 0 & AF_t == 0)
 *
 * compares with
 *
 * (4) P(D | AF_c > 0 || AF_t > 0)
 *
 * This is effectively asking for the value in the upper left vs. the sum of all cells.
 *
 * This class implements the conditional likelihoods summation for any number of alt
 * alleles, where each alt allele has its EXACT probability of segregating calculated by
 * reducing each alt B into the case XB and computing P(D | AF_b > 0 ) as follows:
 *
 * Suppose we have for a A/B/C site the following GLs:
 *
 * AA AB BB AC BC CC
 *
 * and we want to get the bi-allelic GLs for X/B, where X is everything not B
 *
 * XX = AA + AC + CC (since X = A or C)
 * XB = AB + BC
 * BB = BB
 *
 * After each allele has its probability calculated we compute the joint posterior
 * as P(D | AF_* == 0) = prod_i P (D | AF_i == 0), after applying the theta^i
 * prior for the ith least likely allele.
 */
 public final class IndependentAllelesDiploidExactAFCalculator extends ExactAFCalculator {

    private static final int[] BIALLELIC_NON_INFORMATIVE_PLS = {0,0,0};
    private static final List<Allele> BIALLELIC_NOCALL = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
    private static final double PHRED_2_LOG10_COEFF = -.1;

    /**
     * Sorts AFCalcResults by their posteriors of AF > 0, so the
     */
    private static final Comparator<AFCalculationResult> compareAFCalcResultsByPNonRef = Comparator.<AFCalculationResult>comparingDouble(o -> o.getLog10PosteriorOfAFGT0()).reversed();

    /**
     * The AFCalc model we are using to do the bi-allelic computation
     */
    private final AFCalculator biAlleleExactModel;
    private final GenotypeLikelihoodCalculators calculators;

    IndependentAllelesDiploidExactAFCalculator() {
        biAlleleExactModel = new ReferenceDiploidExactAFCalculator();
        calculators = new GenotypeLikelihoodCalculators();
    }

    @Override
    public AFCalculationResult computeLog10PNonRef(final VariantContext vc,
                                                   final int defaultPloidy,
                                                   final double[] log10AlleleFrequencyPriors,
                                                   final StateTracker stateTracker) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(log10AlleleFrequencyPriors, "log10AlleleFrequencyPriors is null");
        Utils.nonNull(stateTracker, "stateTracker is null");

        final List<AFCalculationResult> independentResultTrackers = computeAlleleIndependentExact(vc, defaultPloidy, log10AlleleFrequencyPriors);

        if ( independentResultTrackers.isEmpty() ) {
            throw new IllegalStateException("Independent alleles model returned an empty list of results at VC " + vc);
        }

        if ( independentResultTrackers.size() == 1 ) {
            // fast path for the very common bi-allelic use case
            return independentResultTrackers.get(0);
        } else {
            final AFCalculationResult combinedAltAllelesResult = combineAltAlleleIndependentExact(vc,defaultPloidy,log10AlleleFrequencyPriors);
            // we are a multi-allelic, so we need to actually combine the results
            final List<AFCalculationResult> withMultiAllelicPriors = applyMultiAllelicPriors(independentResultTrackers);
            return combineIndependentPNonRefs(vc, withMultiAllelicPriors, combinedAltAllelesResult);
        }
    }

    private AFCalculationResult combineAltAlleleIndependentExact(final VariantContext vc, final int defaultPloidy, final double[] log10AlleleFrequencyPriors) {
        final VariantContext combinedAltAllelesVariantContext = makeCombinedAltAllelesVariantContext(vc);
        return biAlleleExactModel.getLog10PNonRef(combinedAltAllelesVariantContext, defaultPloidy, vc.getNAlleles() - 1, log10AlleleFrequencyPriors);
    }

    private VariantContext makeCombinedAltAllelesVariantContext(final VariantContext vc) {
        final int nAltAlleles = vc.getNAlleles() - 1;

        if ( nAltAlleles == 1 ) {
            return vc;
        }
        final VariantContextBuilder vcb = new VariantContextBuilder(vc);
        final Allele reference = vcb.getAlleles().get(0);
        vcb.alleles(Arrays.asList(reference, Allele.NON_REF_ALLELE));
        final int genotypeCount = calculators.genotypeCount(2, vc.getNAlleles());
        final double[] hetLikelihoods = new double[vc.getNAlleles() - 1];
        final double[] homAltLikelihoods = new double[genotypeCount - hetLikelihoods.length - 1];
        final double[] newLikelihoods = new double[3];
        final List<Genotype> newGenotypes = new ArrayList<>(vc.getNSamples());
        for (final Genotype oldGenotype : vc.getGenotypes()) {
            final GenotypeBuilder gb = new GenotypeBuilder(oldGenotype);
            final List<Allele> oldAlleles = oldGenotype.getAlleles();
            if (oldAlleles != null) {
                final List<Allele> newAlleles = new ArrayList<>(oldAlleles.size());
                for (int i = 0; i < oldAlleles.size(); i++) {
                    final Allele oldAllele = oldAlleles.get(i);
                    if (oldAllele.isReference()) {
                        newAlleles.add(reference);
                    } else if (oldAllele.isNoCall()) {
                        newAlleles.add(Allele.NO_CALL);
                    } else {
                        newAlleles.add(Allele.NON_REF_ALLELE);
                    }
                }
                gb.alleles(newAlleles);
            }
            if (oldGenotype.isNonInformative()) {
                gb.PL(BIALLELIC_NON_INFORMATIVE_PLS);
            } else if (combineAltAlleleLikelihoods(oldGenotype, genotypeCount, newLikelihoods, hetLikelihoods, homAltLikelihoods)) {
                gb.PL(newLikelihoods);
            }

            newGenotypes.add(gb.make());
        }
        return vcb.genotypes(newGenotypes).make();
    }

    /**
     * Compute the conditional exact AFCalcResult for each allele in vc independently, returning
     * the result of each, in order of the alt alleles in VC
     *
     * @param vc the VariantContext we want to analyze, with at least 1 alt allele
     * @param log10AlleleFrequencyPriors the priors
     * @return a list of the AFCalcResults for each bi-allelic sub context of vc
     */
    private List<AFCalculationResult> computeAlleleIndependentExact(final VariantContext vc, final int defaultPloidy,
                                                                    final double[] log10AlleleFrequencyPriors) {
        final List<AFCalculationResult> results = new LinkedList<>();

        for ( final VariantContext subvc : makeAlleleConditionalContexts(vc) ) {
            final AFCalculationResult resultTracker = biAlleleExactModel.getLog10PNonRef(subvc, defaultPloidy, vc.getNAlleles() - 1, log10AlleleFrequencyPriors);
            results.add(resultTracker);
        }

        return results;
    }

    /**
     * Returns the bi-allelic variant context for each alt allele in vc with bi-allelic likelihoods, in order
     *
     * @param vc the variant context to split.  Must have n.alt.alleles > 1
     * @return a bi-allelic variant context for each alt allele in vc
     */
    @VisibleForTesting
    static List<VariantContext> makeAlleleConditionalContexts(final VariantContext vc) {
        final int nAltAlleles = vc.getNAlleles() - 1;

        if ( nAltAlleles == 1 ) {
            // fast path for bi-allelic case.
            return Collections.singletonList(vc);
        } else {
            // go through the work of ripping up the VC into its biallelic components
            final List<VariantContext> vcs = new LinkedList<>();

            for ( int altI = 0; altI < nAltAlleles; altI++ ) {
                vcs.add(biallelicCombinedGLs(vc, altI + 1));
            }

            return vcs;
        }
    }

    /**
     * Create a single bi-allelic variant context from rootVC with alt allele with index altAlleleIndex
     *
     * @param rootVC the root (potentially multi-allelic) variant context
     * @param altAlleleIndex index of the alt allele, from 0 == first alt allele
     * @return a bi-allelic variant context based on rootVC
     */
    private static VariantContext biallelicCombinedGLs(final VariantContext rootVC, final int altAlleleIndex) {
        if ( rootVC.isBiallelic() ) {
            return rootVC;
        }
        final int nAlts = rootVC.getNAlleles() - 1;
        final List<Genotype> biallelicGenotypes = new ArrayList<>(rootVC.getNSamples());
        for ( final Genotype g : rootVC.getGenotypes() ) {
            biallelicGenotypes.add(combineGLsPrecise(g, altAlleleIndex, nAlts));
        }

        final VariantContextBuilder vcb = new VariantContextBuilder(rootVC);
        final Allele altAllele = rootVC.getAlternateAllele(altAlleleIndex - 1);
        vcb.alleles(Arrays.asList(rootVC.getReference(), altAllele));
        vcb.genotypes(biallelicGenotypes);
        return vcb.make();
    }

    /**
     * Returns a new Genotype with the PLs of the multi-allelic original reduced to a bi-allelic case.
     *
     * <p>Uses the log-sum-exp trick in order to work well with very low PLs</p>
     *
     * <p>This is handled in the following way:</p>
     *
     * <p>Suppose we have for a A/B/C site the following GLs:</p>
     *
     * <p>AA AB BB AC BC CC</p>
     *
     * <p>and we want to get the bi-allelic GLs for X/B, where X is everything not B</p>
     *
     * <p>XX = AA + AC + CC (since X = A or C)<br/>
     * XB = AB + BC                           <br/>
     * BB = BB                                     <br/>
     * </p>
     * <p>
     *     This implementation use the log sum trick in order to avoid numeric inestability.
     * </p>
     *
     * @param original the original multi-allelic genotype
     * @param altIndex the index of the alt allele we wish to keep in the bialleic case -- with ref == 0
     * @param nAlts the total number of alt alleles
     * @return a new biallelic genotype with appropriate PLs
     */
    @VisibleForTesting
    static Genotype combineGLsPrecise(final Genotype original, final int altIndex, final int nAlts) {

        if ( original.isNonInformative() ) {
            return new GenotypeBuilder(original).PL(BIALLELIC_NON_INFORMATIVE_PLS).alleles(BIALLELIC_NOCALL).make();
        }

        if ( altIndex < 1 || altIndex > nAlts ) {
            throw new IllegalStateException("altIndex must be between 1 and nAlts " + nAlts);
        }

        final int[] pls = original.getPL();

        final int nAlleles = nAlts + 1;

        final int plCount = pls.length;

        double BB = 0;
        final double[] XBvalues = new double[nAlleles - 1];
        final double[] XXvalues = new double[plCount - nAlleles];

        int xbOffset = 0;
        int xxOffset = 0;
        for ( int index = 0; index < plCount; index++ ) {
            final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair pair = GenotypeLikelihoods.getAllelePair(index);
            final int i = pair.alleleIndex1;
            final int j = pair.alleleIndex2;
            if (i == j) {
              if (i == altIndex) {
                  BB = PHRED_2_LOG10_COEFF * pls[index];
              } else {
                  XXvalues[xxOffset++] = PHRED_2_LOG10_COEFF * pls[index];
              }
            } else if (i == altIndex || j == altIndex) {
                XBvalues[xbOffset++] = PHRED_2_LOG10_COEFF * pls[index];
            } else {
                XXvalues[xxOffset++] = PHRED_2_LOG10_COEFF * pls[index];
            }
        }

        final double XB = MathUtils.log10SumLog10(XBvalues);
        final double XX = MathUtils.log10SumLog10(XXvalues);

        final double[] GLs = { XX, XB, BB};
        return new GenotypeBuilder(original).PL(GLs).alleles(BIALLELIC_NOCALL).make();
    }

    @VisibleForTesting
    static List<AFCalculationResult> applyMultiAllelicPriors(final List<AFCalculationResult> conditionalPNonRefResults) {
        final List<AFCalculationResult> sorted = new ArrayList<>(conditionalPNonRefResults);

        // sort the results, so the most likely allele is first
        Collections.sort(sorted, compareAFCalcResultsByPNonRef);

        final double lastPosteriorGt0 = sorted.get(0).getLog10PosteriorOfAFGT0();
        final double log10SingleAllelePriorOfAFGt0 = conditionalPNonRefResults.get(0).getLog10PriorOfAFGT0();

        for ( int i = 0; i < sorted.size(); i++ ) {
            if ( sorted.get(i).getLog10PosteriorOfAFGT0() > lastPosteriorGt0 ) {
                throw new IllegalStateException("pNonRefResults not sorted: lastPosteriorGt0 " + lastPosteriorGt0 + " but current is " + sorted.get(i).getLog10PosteriorOfAFGT0());
            }

                final double log10PriorAFGt0 = (i + 1) * log10SingleAllelePriorOfAFGt0;
            final double log10PriorAFEq0 = Math.log10(1 - Math.pow(10, log10PriorAFGt0));
            final double[] thetaTONPriors = new double[] { log10PriorAFEq0, log10PriorAFGt0 };

            // bind pNonRef for allele to the posterior value of the AF > 0 with the new adjusted prior
            sorted.set(i, sorted.get(i).copyWithNewPriors(MathUtils.normalizeLog10(thetaTONPriors)));
        }

        return sorted;
    }

    /**
     * Take the independent estimates of pNonRef for each alt allele and combine them into a single result
     *
     * Given n independent calculations for each of n alternate alleles create a single
     * combined AFCalcResult with:
     *
     * priors for AF == 0 equal to theta^N for the nth least likely allele
     * posteriors that reflect the combined chance that any alleles are segregating and corresponding
     * likelihoods
     * combined MLEs in the order of the alt alleles in vc
     *
     * @param sortedResultsWithThetaNPriors the pNonRef result for each allele independently
     */
    private static AFCalculationResult combineIndependentPNonRefs(final VariantContext vc,
                                                                  final List<AFCalculationResult> sortedResultsWithThetaNPriors,
                                                                  final AFCalculationResult combinedAltAllelesResult) {


        final int nAltAlleles = sortedResultsWithThetaNPriors.size();
        final int[] alleleCountsOfMLE = new int[nAltAlleles];
        final Map<Allele, Double> log10pRefByAllele = new LinkedHashMap<>(nAltAlleles);

        // the sum of the log10 posteriors for AF == 0 and AF > 0 to determine joint probs

        for ( final AFCalculationResult sortedResultWithThetaNPriors : sortedResultsWithThetaNPriors ) {
            final Allele altAllele = sortedResultWithThetaNPriors.getAllelesUsedInGenotyping().get(1);
            final int altI = vc.getAlleles().indexOf(altAllele) - 1;

            // MLE of altI allele is simply the MLE of this allele in altAlleles
            alleleCountsOfMLE[altI] = sortedResultWithThetaNPriors.getAlleleCountAtMLE(altAllele);

            // bind pNonRef for allele to the posterior value of the AF > 0 with the new adjusted prior
            log10pRefByAllele.put(altAllele, sortedResultWithThetaNPriors.getLog10PosteriorOfAFEq0());
        }

        return new AFCalculationResult(alleleCountsOfMLE, vc.getAlleles(),
                // necessary to ensure all values < 0
                MathUtils.normalizeLog10(new double[] { combinedAltAllelesResult.getLog10LikelihoodOfAFEq0(), combinedAltAllelesResult.getLog10LikelihoodOfAFGT0() }),
                // priors incorporate multiple alt alleles, must be normalized
                MathUtils.normalizeLog10(new double[] { combinedAltAllelesResult.getLog10PriorOfAFEq0(), combinedAltAllelesResult.getLog10PriorOfAFGT0() }),
                log10pRefByAllele);
    }

    private static boolean combineAltAlleleLikelihoods(final Genotype g, final int plMaxIndex, final double[] dest,
                                                       final double[] hetLikelihoods, final double[] homAltLikelihoods) {

        final int[] pls = g.getPL();
        if (pls == null) {
            return false;
        }
        int hetNextIndex = 0;
        int homAltNextIndex = 0;
        for (int plIndex = 1; plIndex < plMaxIndex; plIndex++) {
            final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(plIndex);
            if (alleles.alleleIndex1 == 0 || alleles.alleleIndex2 == 0) {
                hetLikelihoods[hetNextIndex++] = pls[plIndex] * PHRED_2_LOG10_COEFF;
            } else {
                homAltLikelihoods[homAltNextIndex++] = pls[plIndex] * PHRED_2_LOG10_COEFF;
            }
        }
        dest[0] = pls[0] * PHRED_2_LOG10_COEFF;
        dest[1] = MathUtils.approximateLog10SumLog10(hetLikelihoods);
        dest[2] = MathUtils.approximateLog10SumLog10(homAltLikelihoods);
        return true;
    }
}
