package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Stores the most likely and second most likely alleles, along with a threshold
 * for assuming computing that a read is informative.
 *
 * If the difference between the most-likely allele and the next-most-likely allele is < INFORMATIVE_LIKELIHOOD_THRESHOLD
 * then the most likely allele is set to "no call", and isInformative will return false.  This constant can be
 * overridden simply by using one of the version of these calls that accepts informative threshold as an argument.
 *
 * For convenience, there are functions called getAlleleIfInformative that return either the most likely allele, or
 * NO_CALL if two or more alleles have likelihoods within INFORMATIVE_LIKELIHOOD_THRESHOLD of one another.
 *
 * By default empty allele maps will return NO_CALL, and allele maps with a single entry will return the
 * corresponding key
 */
public final class MostLikelyAllele {
    public static final double INFORMATIVE_LIKELIHOOD_THRESHOLD = 0.2;

    public static final MostLikelyAllele NO_CALL = new MostLikelyAllele(Allele.NO_CALL, null, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);

    private final Allele mostLikely;
    private final Allele secondMostLikely;
    private final double log10LikelihoodOfMostLikely;
    private final double log10LikelihoodOfSecondBest;

    /**
     * Create a new MostLikelyAllele
     *
     * If there's a meaningful most likely allele, allele should be a real allele.  If none can be determined,
     * mostLikely should be a NO_CALL allele.
     *
     * @param mostLikely the most likely allele
     * @param secondMostLikely the most likely allele after mostLikely (may be null)
     * @param log10LikelihoodOfMostLikely the log10 likelihood of the most likely allele      (or {@link Double.NEGATIVE_INFINITY} if none is available)
     * @param log10LikelihoodOfSecondBest the log10 likelihood of the next most likely allele (or {@link Double.NEGATIVE_INFINITY} if none is available)
     */
    public MostLikelyAllele(final Allele mostLikely, final Allele secondMostLikely, final double log10LikelihoodOfMostLikely, final double log10LikelihoodOfSecondBest) {
        Utils.nonNull( mostLikely, "mostLikely allele cannot be null");
        if ( mostLikely.equals(secondMostLikely)){
            throw new IllegalArgumentException("most likely allele and second most likely allele should be different");
        };
        if ( log10LikelihoodOfMostLikely != Double.NEGATIVE_INFINITY && ! MathUtils.goodLog10Probability(log10LikelihoodOfMostLikely) ) {
            throw new IllegalArgumentException("log10LikelihoodOfMostLikely must be either -Infinity or a good log10 prob but got " + log10LikelihoodOfMostLikely);
        }
        if ( log10LikelihoodOfSecondBest != Double.NEGATIVE_INFINITY && ! MathUtils.goodLog10Probability(log10LikelihoodOfSecondBest) ) {
            throw new IllegalArgumentException("log10LikelihoodOfSecondBest must be either -Infinity or a good log10 prob but got " + log10LikelihoodOfSecondBest);
        }
        if ( log10LikelihoodOfMostLikely < log10LikelihoodOfSecondBest ) {
            throw new IllegalArgumentException("log10LikelihoodOfMostLikely must be <= log10LikelihoodOfSecondBest but got " + log10LikelihoodOfMostLikely + " vs 2nd " + log10LikelihoodOfSecondBest);
        }

        this.mostLikely = mostLikely;
        this.secondMostLikely = secondMostLikely;
        this.log10LikelihoodOfMostLikely = log10LikelihoodOfMostLikely;
        this.log10LikelihoodOfSecondBest = log10LikelihoodOfSecondBest;
    }

    /**
     * Create a new MostLikelyAllele with only 1 allele.
     *
     * If there's a meaningful most likely allele, allele should be a real allele.  If none can be determined,
     * mostLikely should be a NO_CALL allele.
     *
     * @param mostLikely the most likely allele
     * @param log10LikelihoodOfMostLikely the log10 likelihood of the most likely allele      (or {@link Double.NEGATIVE_INFINITY} if none is available)
     */
    public MostLikelyAllele(final Allele mostLikely, final double log10LikelihoodOfMostLikely) {
        this(mostLikely, null, log10LikelihoodOfMostLikely, Double.NEGATIVE_INFINITY);
    }

    /**
     * Retruns the most likely allele. Never null.
     */
    public Allele getMostLikelyAllele() {
        return mostLikely;
    }

    /**
     * Retruns the second most likely allele or null if there is none.
     */
    public Allele getSecondMostLikelyAllele() {
        return secondMostLikely;
    }

    /**
     * Retruns the log10 likelihood of the most likely allele or {@link Double.NEGATIVE_INFINITY} if none is available.
     */
    public double getLog10LikelihoodOfMostLikely() {
        return log10LikelihoodOfMostLikely;
    }

    /**
     * Retruns the log10 likelihood of the second most likely allele or {@link Double.NEGATIVE_INFINITY} if none is available.
     */
    public double getLog10LikelihoodOfSecondBest() {
        return log10LikelihoodOfSecondBest;
    }

    /**
     * @see #isInformative(double) with threshold of INFORMATIVE_LIKELIHOOD_THRESHOLD
     */
    public boolean isInformative() {
        return isInformative(INFORMATIVE_LIKELIHOOD_THRESHOLD);
    }

    /**
     * Was this allele selected from an object that was specifically informative about the allele?
     *
     * The calculation that implements this is whether the likelihood of the most likely allele is larger
     * than the second most likely by at least the log10ThresholdForInformative
     *
     * @return true if so, false if not
     */
    public boolean isInformative(final double log10ThresholdForInformative) {
        return log10LikelihoodOfMostLikely - log10LikelihoodOfSecondBest > log10ThresholdForInformative;
    }

    /**
     * @see #getAlleleIfInformative(double) with threshold of INFORMATIVE_LIKELIHOOD_THRESHOLD
     */
    public Allele getAlleleIfInformative() {
        return getAlleleIfInformative(INFORMATIVE_LIKELIHOOD_THRESHOLD);
    }

    /**
     * Get the most likely allele if isInformative(log10ThresholdForInformative) is true, or NO_CALL otherwise
     *
     * @param log10ThresholdForInformative a log10 threshold to determine if the most likely allele was informative
     * @return a non-null allele
     */
    public Allele getAlleleIfInformative(final double log10ThresholdForInformative) {
        return isInformative(log10ThresholdForInformative) ? mostLikely : Allele.NO_CALL;
    }
}
