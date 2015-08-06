package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.fragments.FragmentCollection;
import org.broadinstitute.hellbender.utils.fragments.FragmentUtils;
import org.broadinstitute.hellbender.utils.genotyper.DiploidGenotype;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.List;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

/**
 * Stable, error checking version of the Bayesian genotyper.  Useful for calculating the likelihoods, priors,
 * and posteriors given a pile of bases and quality scores
 *
 * Suppose we have bases b1, b2, ..., bN with qualities scores q1, q2, ..., qN.  This object
 * calculates:
 *
 * P(G | D) = P(G) * P(D | G)
 *
 * where
 *
 * P(D | G) = sum_i log10 P(bi | G)
 *
 * and
 *
 * P(bi | G) = 1 - P(error | q1) if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for homozygous genotypes and for heterozygous genotypes:
 *
 * P(bi | G) = 1 - P(error | q1) / 2 + P(error | q1) / 6 if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for each of the 10 unique diploid genotypes AA, AC, AG, .., TT
 *
 * Everything is stored as arrays indexed by DiploidGenotype.ordinal() values in log10 space.
 *
 * The priors contain the relative probabilities of each genotype, and must be provided at object creation.
 * From then on, you can call any of the add() routines to update the likelihoods and posteriors in the above
 * model.
 */
public class DiploidSNPGenotypeLikelihoods {

    public final static double DEFAULT_PCR_ERROR_RATE = FragmentUtils.DEFAULT_PCR_ERROR_RATE;

    protected final static int FIXED_PLOIDY = 2;
    protected final static int MAX_PLOIDY = FIXED_PLOIDY + 1;
    protected final static double ploidyAdjustment = log10(FIXED_PLOIDY);
    protected final static double log10_3 = log10(3.0);

    protected boolean VERBOSE = false;

    //
    // The fundamental data arrays associated with a Genotype Likelihoods object
    //
    protected double[] log10Likelihoods = null;

    // TODO: don't calculate this each time through
    protected double log10_PCR_error_3;
    protected double log10_1_minus_PCR_error;

    /**
     * Create a new GenotypeLikelhoods object with given PCR error rate for each diploid genotype
     *
     * @param PCR_error_rate  the PCR error rate
     */
    public DiploidSNPGenotypeLikelihoods(final double PCR_error_rate) {
        log10_PCR_error_3 = log10(PCR_error_rate) - log10_3;
        log10_1_minus_PCR_error = log10(1.0 - PCR_error_rate);
        setToZero();
    }

    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        final DiploidSNPGenotypeLikelihoods c = (DiploidSNPGenotypeLikelihoods)super.clone();
        c.log10Likelihoods = log10Likelihoods.clone();
        return c;
    }

    protected void setToZero() {
        log10Likelihoods = genotypeZeros.clone();                 // likelihoods are all zeros
    }

    /**
     * Returns an array of log10 likelihoods for each genotype, indexed by DiploidGenotype.ordinal values()
     * @return likelihoods array
     */
    public double[] getLikelihoods() {
        return log10Likelihoods;
    }

    // -------------------------------------------------------------------------------------
    //
    // add() routines.  These are the workhorse routines for calculating the overall genotype
    // likelihoods given observed bases and reads.  Includes high-level operators all the
    // way down to single base and qual functions.
    //
    // -------------------------------------------------------------------------------------

    /**
     * Updates likelihoods and posteriors to reflect the additional observations contained within the
     * read-based pileup up by calling add(observedBase, qualityScore) for each base / qual in the
     * pileup
     *
     * @param pileup                    read pileup
     * @param ignoreBadBases            should we ignore bad bases?
     * @param capBaseQualsAtMappingQual should we cap a base's quality by its read's mapping quality?
     * @param minBaseQual               the minimum base quality at which to consider a base valid
     * @return the number of good bases found in the pileup
     */
    public int add(final ReadPileup pileup, final boolean ignoreBadBases, final boolean capBaseQualsAtMappingQual, final int minBaseQual) {
        int n = 0;

        // for each fragment, add to the likelihoods
        final FragmentCollection<PileupElement> fpile = pileup.toFragments();

        for ( final PileupElement p : fpile.getSingletonReads() ) {
            n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);
        }

        for ( final List<PileupElement> overlappingPair : fpile.getOverlappingPairs() ) {
            n += add(overlappingPair, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);
        }

        return n;
    }

    public int add(final PileupElement elt, final boolean ignoreBadBases, final boolean capBaseQualsAtMappingQual, final int minBaseQual) {
        final byte obsBase = elt.getBase();
        final byte qual = qualToUse(elt, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);
        if ( qual == 0 ) {
            return 0;
        }

        return add(obsBase, qual, (byte)0, (byte)0, 1);
    }

    public int add(final List<PileupElement> overlappingPair, final boolean ignoreBadBases, final boolean capBaseQualsAtMappingQual, final int minBaseQual) {
        final PileupElement p1 = overlappingPair.get(0);
        final PileupElement p2 = overlappingPair.get(1);

        final byte observedBase1 = p1.getBase();
        final byte qualityScore1 = qualToUse(p1, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);
        final byte observedBase2 = p2.getBase();
        final byte qualityScore2 = qualToUse(p2, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);

        if ( qualityScore1 == 0 ) {
            if ( qualityScore2 == 0 ) // abort early if we didn't see any good bases
            {
                return 0;
            } else {
                return add(observedBase2, qualityScore2, (byte)0, (byte)0);
            }
        } else {
            return add(observedBase1, qualityScore1, observedBase2, qualityScore2);
        }
    }

    /**
     *
     * @param obsBase1    first observed base
     * @param qual1       base qual of first observed base
     * @param obsBase2    second observed base
     * @param qual2       base qual of second observed base; can be 0, indicating no second base was observed for this fragment
     * @param nObs        the number of times this quad of values was seen.  Generally 1, but reduced reads can have nObs > 1 for synthetic reads
     * @return 0 if the base is bad, 1 otherwise
     */
    private int add(final byte obsBase1, final byte qual1, final byte obsBase2, final byte qual2, final int nObs) {
        // TODO-- Right now we assume that there are at most 2 reads per fragment.  This assumption is fine
        // TODO--   given the current state of next-gen sequencing, but may need to be fixed in the future.
        // TODO--   However, when that happens, we'll need to be a lot smarter about the caching we do here.

        // Just look up the cached result if it's available, or compute and store it
        final DiploidSNPGenotypeLikelihoods gl;
        if ( ! inCache(obsBase1, qual1, obsBase2, qual2, FIXED_PLOIDY) ) {
            gl = calculateCachedGenotypeLikelihoods(obsBase1, qual1, obsBase2, qual2, FIXED_PLOIDY);
        } else {
            gl = getCachedGenotypeLikelihoods(obsBase1, qual1, obsBase2, qual2, FIXED_PLOIDY);
        }

        // for bad bases, there are no likelihoods
        if ( gl == null ) {
            return 0;
        }

        final double[] likelihoods = gl.getLikelihoods();

        for ( final DiploidGenotype g : DiploidGenotype.values() ) {
            final double likelihood = likelihoods[g.ordinal()];
            log10Likelihoods[g.ordinal()] += likelihood * nObs;
        }

        return 1;
    }

    private int add(final byte obsBase1, final byte qual1, final byte obsBase2, final byte qual2) {
        return add(obsBase1, qual1, obsBase2, qual2, 1);
    }

    // -------------------------------------------------------------------------------------
    //
    // Dealing with the cache routines
    //
    // -------------------------------------------------------------------------------------

    static DiploidSNPGenotypeLikelihoods[][][][][] CACHE = new DiploidSNPGenotypeLikelihoods[BaseUtils.BASES.length][QualityUtils.MAX_SAM_QUAL_SCORE +1][BaseUtils.BASES.length+1][QualityUtils.MAX_SAM_QUAL_SCORE +1][MAX_PLOIDY];

    protected boolean inCache(final byte observedBase1, final byte qualityScore1, final byte observedBase2, final byte qualityScore2, final int ploidy) {
        return getCache(CACHE, observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy) != null;
    }

    protected DiploidSNPGenotypeLikelihoods getCachedGenotypeLikelihoods(final byte observedBase1, final byte qualityScore1, final byte observedBase2, final byte qualityScore2, final int ploidy) {
        final DiploidSNPGenotypeLikelihoods gl = getCache(CACHE, observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy);
        if ( gl == null ) {
            throw new RuntimeException(String.format("BUG: trying to fetch an unset cached genotype likelihood at base1=%c, qual1=%d, base2=%c, qual2=%d, ploidy=%d",
                    observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy));
        }
        return gl;
    }

    protected DiploidSNPGenotypeLikelihoods calculateCachedGenotypeLikelihoods(final byte observedBase1, final byte qualityScore1, final byte observedBase2, final byte qualityScore2, final int ploidy) {
        final DiploidSNPGenotypeLikelihoods gl = calculateGenotypeLikelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2);
        setCache(CACHE, observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy, gl);
        return gl;
    }

    protected void setCache( final DiploidSNPGenotypeLikelihoods[][][][][] cache,
                             final byte observedBase1, final byte qualityScore1, final byte observedBase2, final byte qualityScore2, final int ploidy,
                             final DiploidSNPGenotypeLikelihoods val ) {
        final int i = BaseUtils.simpleBaseToBaseIndex(observedBase1);
        final int j = qualityScore1;
        final int k = qualityScore2 != 0 ? BaseUtils.simpleBaseToBaseIndex(observedBase2) : BaseUtils.BASES.length;
        final int l = qualityScore2;
        final int m = ploidy;

        cache[i][j][k][l][m] = val;
    }

    protected DiploidSNPGenotypeLikelihoods getCache(final DiploidSNPGenotypeLikelihoods[][][][][] cache,
                                            final byte observedBase1, final byte qualityScore1, final byte observedBase2, final byte qualityScore2, final int ploidy) {
        final int i = BaseUtils.simpleBaseToBaseIndex(observedBase1);
        final int j = qualityScore1;
        final int k = qualityScore2 != 0 ? BaseUtils.simpleBaseToBaseIndex(observedBase2) : BaseUtils.BASES.length;
        final int l = qualityScore2;
        final int m = ploidy;
        return cache[i][j][k][l][m];
    }

    protected DiploidSNPGenotypeLikelihoods calculateGenotypeLikelihoods(final byte observedBase1, final byte qualityScore1, final byte observedBase2, final byte qualityScore2) {
        final double[] log10FourBaseLikelihoods = computeLog10Likelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2);

        try {

            final DiploidSNPGenotypeLikelihoods gl = (DiploidSNPGenotypeLikelihoods)this.clone();
            gl.setToZero();

            // we need to adjust for ploidy.  We take the raw p(obs | chrom) / ploidy, which is -log10(ploidy) in log space
            for ( final DiploidGenotype g : DiploidGenotype.values() ) {

                // todo assumes ploidy is 2 -- should be generalized.  Obviously the below code can be turned into a loop
                double p_base = 0.0;
                p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base1)] - ploidyAdjustment);
                p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base2)] - ploidyAdjustment);

                final double likelihood = log10(p_base);
                gl.log10Likelihoods[g.ordinal()] += likelihood;
            }

            if ( VERBOSE ) {
                for ( final DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
                System.out.println();
            for ( final DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]); }
                System.out.println();
            }

            return gl;

         } catch ( final CloneNotSupportedException e ) {
             throw new RuntimeException(e);
         }
    }

    /**
     * Updates likelihoods and posteriors to reflect an additional observation of observedBase with
     * qualityScore.
     *
     * @param observedBase1  the base observed on the 1st read of the fragment
     * @param qualityScore1  the qual of the base on the 1st read of the fragment, or zero if NA
     * @param observedBase2  the base observed on the 2nd read of the fragment
     * @param qualityScore2  the qual of the base on the 2nd read of the fragment, or zero if NA
     * @return likelihoods for this observation or null if the base was not considered good enough to add to the likelihoods (Q0 or 'N', for example)
     */
    protected double[] computeLog10Likelihoods(final byte observedBase1, final byte qualityScore1, final byte observedBase2, final byte qualityScore2) {
        final double[] log10FourBaseLikelihoods = baseZeros.clone();

        for ( final byte trueBase : BaseUtils.BASES ) {
            double likelihood = 0.0;

            for ( final byte fragmentBase : BaseUtils.BASES ) {
                double log10FragmentLikelihood = (trueBase == fragmentBase ? log10_1_minus_PCR_error : log10_PCR_error_3);
                if ( qualityScore1 != 0 ) {
                    log10FragmentLikelihood += log10PofObservingBaseGivenChromosome(observedBase1, fragmentBase, qualityScore1);
                }
                if ( qualityScore2 != 0 ) {
                    log10FragmentLikelihood += log10PofObservingBaseGivenChromosome(observedBase2, fragmentBase, qualityScore2);
                }

                //if ( VERBOSE ) {
                //    System.out.printf("  L(%c | b=%s, Q=%d) = %f / %f%n",
                //            observedBase, trueBase, qualityScore, pow(10,likelihood) * 100, likelihood);
                //}

                likelihood += pow(10, log10FragmentLikelihood);
            }

            log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(trueBase)] = log10(likelihood);
        }

        return log10FourBaseLikelihoods;
    }

    /**
     *
     * @param observedBase observed base
     * @param chromBase    target base
     * @param qual         base quality
     * @return log10 likelihood
     */
    protected double log10PofObservingBaseGivenChromosome(final byte observedBase, final byte chromBase, final byte qual) {

        final double logP;

        if ( observedBase == chromBase ) {
            // the base is consistent with the chromosome -- it's 1 - e
            //logP = oneMinusData[qual];
            final double e = pow(10, (qual / -10.0));
            logP = log10(1.0 - e);
        } else {
            // the base is inconsistent with the chromosome -- it's e * P(chromBase | observedBase is an error)
            logP = qual / -10.0 + (-log10_3);
        }

        //System.out.printf("%c %c %d => %f%n", observedBase, chromBase, qual, logP);
        return logP;
    }

    /**
     * Helper function that returns the phred-scaled base quality score we should use for calculating
     * likelihoods for a pileup element.  May return 0 to indicate that the observation is bad, and may
     * cap the quality score by the mapping quality of the read itself.
     *
     * @param p                           Pileup element
     * @param ignoreBadBases              Should we ignore bad bases?
     * @param capBaseQualsAtMappingQual   Should we cap the base qualities at the mapping quality of the read?
     * @param minBaseQual                 Minimum allowed base quality
     * @return the actual base quality to use
     */
    private static byte qualToUse(final PileupElement p, final boolean ignoreBadBases, final boolean capBaseQualsAtMappingQual, final int minBaseQual) {
        if ( ignoreBadBases && !BaseUtils.isRegularBase( p.getBase() ) ) {
            return 0;
        }

        byte qual = p.getQual();

        if ( qual > SAMUtils.MAX_PHRED_SCORE ) {
            throw new UserException.MisencodedBAM(p.getRead(), "we encountered an extremely high quality score (" + (int) qual + ")");
        }
        if ( capBaseQualsAtMappingQual ) {
            qual = (byte) Math.min(0xff & qual, p.getMappingQual());
        }
        if ( (int)qual < minBaseQual ) {
            qual = (byte) 0;
        }

        return qual;
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // helper routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Return a string representation of this object in a moderately usable form
     *
     * @return string representation
     */
    public String toString() {
        double sum = 0;
        final StringBuilder s = new StringBuilder();
        for (final DiploidGenotype g : DiploidGenotype.values()) {
            s.append(String.format("%s %.10f ", g, log10Likelihoods[g.ordinal()]));
			sum += Math.pow(10, log10Likelihoods[g.ordinal()]);
        }
		s.append(String.format(" %f", sum));
        return s.toString();
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // Validation routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    public boolean validate() {
        return validate(true);
    }

    public boolean validate(final boolean throwException) {
        try {
            for ( final DiploidGenotype g : DiploidGenotype.values() ) {
                String bad = null;

                final int i = g.ordinal();
                if ( ! MathUtils.wellFormedDouble(log10Likelihoods[i]) || ! MathUtils.isNegativeOrZero(log10Likelihoods[i]) ) {
                    bad = String.format("Likelihood %f is badly formed", log10Likelihoods[i]);
                }

                if ( bad != null ) {
                    throw new IllegalStateException(String.format("At %s: %s", g.toString(), bad));
                }
            }
        } catch ( final IllegalStateException e ) {
            if ( throwException ) {
                throw new RuntimeException(e);
            } else {
                return false;
            }
        }

        return true;
    }

    //
    // Constant static data
    //
    private final static double[] genotypeZeros = new double[DiploidGenotype.values().length];
    private final static double[] baseZeros = new double[BaseUtils.BASES.length];

    static {
        for ( final DiploidGenotype g : DiploidGenotype.values() ) {
            genotypeZeros[g.ordinal()] = 0.0;
        }
        for ( final byte base : BaseUtils.BASES ) {
            baseZeros[BaseUtils.simpleBaseToBaseIndex(base)] = 0.0;
        }
    }
}
