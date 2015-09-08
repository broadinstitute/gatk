package org.broadinstitute.hellbender.utils;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.Serializable;
import java.util.Comparator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

public final class MannWhitneyU {

    private static final NormalDistribution STANDARD_NORMAL = new NormalDistribution(0.0,1.0);
    private static final NormalDistribution APACHE_NORMAL = new NormalDistribution(0.0,1.0,1e-2);
    private static final double LNSQRT2PI = Math.log(Math.sqrt(2.0 * Math.PI));

    private final SortedSet<Pair<Number,USet>> observations;
    private int sizeSet1;
    private int sizeSet2;
    private final ExactMode exactMode;

    public MannWhitneyU(final ExactMode mode, final boolean dither) {
        if ( dither ) {
            observations = new TreeSet<>(new DitheringComparator());
        } else {
            observations = new TreeSet<>(new NumberedPairComparator());
        }
        sizeSet1 = 0;
        sizeSet2 = 0;
        exactMode = mode;
    }

    public MannWhitneyU(final boolean dither) {
        this(ExactMode.POINT, dither);
    }

    /**
     * Add an observation into the observation tree
     * @param n: the observation (a number)
     * @param set: whether the observation comes from set 1 or set 2
     */
    public void add(final Number n, final USet set) {
        observations.add(new MutablePair<>(n, set));
        if ( set == USet.SET1 ) {
            ++sizeSet1;
        } else {
            ++sizeSet2;
        }
    }

    /**
     * Runs the one-sided test under the hypothesis that the data in set "lessThanOther" stochastically
     * dominates the other set
     * @param lessThanOther - either Set1 or Set2
     * @return - u-based z-approximation, and p-value associated with the test (p-value is exact for small n,m)
     */
    public Pair<Double,Double> runOneSidedTest(final USet lessThanOther) {
        final long u = calculateOneSidedU(observations, lessThanOther);
        final int n = lessThanOther == USet.SET1 ? sizeSet1 : sizeSet2;
        final int m = lessThanOther == USet.SET1 ? sizeSet2 : sizeSet1;
        if ( n == 0 || m == 0 ) {
            // test is uninformative as one or both sets have no observations
            return new MutablePair<>(Double.NaN, Double.NaN);
        }

        // the null hypothesis is that {N} is stochastically less than {M}, so U has counted
        // occurrences of {M}s before {N}s. We would expect that this should be less than (n*m+1)/2 under
        // the null hypothesis, so we want to integrate from K=0 to K=U for cumulative cases. Always.
        return calculateP(n, m, u, false, exactMode);
    }


    /**
     * Runs the one-sided test under the hypothesis that the data in set "vals2" stochastically
     * dominates the vals1 set
     * @return - u-based z-approximation, and p-value associated with the test (p-value is exact for small n,m)
     */
    public static Pair<Double,Double> runOneSidedTest(final boolean dithering, final List<? extends Number> vals1, final List<? extends Number> vals2) {
        Utils.nonNull(vals1);
        Utils.nonNull(vals2);
        final MannWhitneyU mannWhitneyU = new MannWhitneyU(dithering);
        for (final Number qual : vals1) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET1);
        }
        for (final Number qual : vals2) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET2);
        }
        return mannWhitneyU.runOneSidedTest(MannWhitneyU.USet.SET1);
    }

    /**
     * Given a u statistic, calculate the p-value associated with it, dispatching to approximations where appropriate
     * @param n - The number of entries in the stochastically smaller (dominant) set
     * @param m - The number of entries in the stochastically larger (dominated) set
     * @param u - the Mann-Whitney U value
     * @param twoSided - is the test twosided
     * @return the (possibly approximate) p-value associated with the MWU test, and the (possibly approximate) z-value associated with it
     * todo -- there must be an approximation for small m and large n
     */
    private static Pair<Double,Double> calculateP(final int n, final int m, final long u, final boolean twoSided, final ExactMode exactMode) {
        final Pair<Double,Double> zandP;
        if ( n > 8 && m > 8 ) {
            // large m and n - normal approx
            zandP = calculatePNormalApproximation(n,m,u, twoSided);
        } else if ( n > 5 && m > 7 ) {
            // large m, small n - sum uniform approx
            // todo -- find the appropriate regimes where this approximation is actually better enough to merit slowness
            // pval = calculatePUniformApproximation(n,m,u);
            zandP = calculatePNormalApproximation(n, m, u, twoSided);
        } else if ( n > 8 || m > 8 ) {
            zandP = calculatePFromTable(n, m, u, twoSided);
        } else {
            // small m and n - full approx
            zandP = calculatePRecursively(n,m,u,twoSided,exactMode);
        }

        return zandP;
    }

    public static Pair<Double,Double> calculatePFromTable(final int n, final int m, final long u, final boolean twoSided) {
        // todo -- actually use a table for:
        // todo      - n large, m small
        return calculatePNormalApproximation(n, m, u, twoSided);
    }

    /**
     * Uses a normal approximation to the U statistic in order to return a cdf p-value. See Mann, Whitney [1947]
     * @param n - The number of entries in the stochastically smaller (dominant) set
     * @param m - The number of entries in the stochastically larger (dominated) set
     * @param u - the Mann-Whitney U value
     * @param twoSided - whether the test should be two sided
     * @return p-value associated with the normal approximation
     */
    public static Pair<Double,Double> calculatePNormalApproximation(final int n, final int m, final long u, final boolean twoSided) {
        final double z = getZApprox(n,m,u);
        if ( twoSided ) {
            return new MutablePair<>(z,2.0*(z < 0 ? STANDARD_NORMAL.cumulativeProbability(z) : 1.0-STANDARD_NORMAL.cumulativeProbability(z)));
        } else {
            return new MutablePair<>(z,STANDARD_NORMAL.cumulativeProbability(z));
        }
    }

    /**
     * Calculates the Z-score approximation of the u-statistic
     * @param n - The number of entries in the stochastically smaller (dominant) set
     * @param m - The number of entries in the stochastically larger (dominated) set
     * @param u - the Mann-Whitney U value
     * @return the asymptotic z-approximation corresponding to the MWU p-value for n < m
     */
    private static double getZApprox(final int n, final int m, final long u) {
        final double mean = ( ((long)m)*n+1.0)/2;
        final double var = (((long) n)*m*(n+m+1.0))/12;
        final double z = ( u - mean )/ Math.sqrt(var);
        return z;
    }

    /**
     * Calculates the U-statistic associated with the one-sided hypothesis that "dominator" stochastically dominates
     * the other U-set. Note that if S1 dominates S2, we want to count the occurrences of points in S2 coming before points in S1.
     * @param observed - the observed data points, tagged by each set
     * @param dominator - the set that is hypothesized to be stochastically dominating
     * @return the u-statistic associated with the hypothesis that dominator stochastically dominates the other set
     */
    public static long calculateOneSidedU(final SortedSet<Pair<Number,USet>> observed, final USet dominator) {
        long otherBeforeDominator = 0l;
        int otherSeenSoFar = 0;
        for ( final Pair<Number,USet> dataPoint : observed ) {
            if ( dataPoint.getRight() != dominator ) {
                ++otherSeenSoFar;
            } else {
                otherBeforeDominator += otherSeenSoFar;
            }
        }

        return otherBeforeDominator;
    }

    /**
     * The Mann-Whitney U statistic follows a recursive equation (that enumerates the proportion of possible
     * binary strings of "n" zeros, and "m" ones, where a one precedes a zero "u" times). This accessor
     * calls into that recursive calculation.
     * @param n: number of set-one entries (hypothesis: set one is stochastically less than set two)
     * @param m: number of set-two entries
     * @param u: number of set-two entries that precede set-one entries (e.g. 0,1,0,1,0 -> 3 )
     * @param twoSided: whether the test is two sided or not. The recursive formula is symmetric, multiply by two for two-sidedness.
     * @param  mode: whether the mode is a point probability, or a cumulative distribution
     * @return the probability under the hypothesis that all sequences are equally likely of finding a set-two entry preceding a set-one entry "u" times.
     */
    public static Pair<Double,Double> calculatePRecursively(final int n, final int m, final long u, final boolean twoSided, final ExactMode mode) {
        if ( m > 8 && n > 5 ) { throw new GATKException(String.format("Please use the appropriate (normal or sum of uniform) approximation. Values n: %d, m: %d", n, m)); }
        final double p = mode == ExactMode.POINT ? cpr(n,m,u) : cumulativeCPR(n,m,u);
        //p *= twoSided ? 2.0 : 1.0;
        final double z;
        if ( mode == ExactMode.CUMULATIVE ) {
            z = APACHE_NORMAL.inverseCumulativeProbability(p);
        } else {
            final double sd = Math.sqrt((1.0 + 1.0 / (1 + n + m)) * (n * m) * (1.0 + n + m) / 12); // biased variance empirically better fit to distribution then asymptotic variance
            //System.out.printf("SD is %f and Max is %f and prob is %f%n",sd,1.0/Math.sqrt(sd*sd*2.0*Math.PI),p);
            if ( p > 1.0/ Math.sqrt(sd * sd * 2.0 * Math.PI) ) { // possible for p-value to be outside the range of the normal. Happens at the mean, so z is 0.
                z = 0.0;
            } else {
                if ( u >= n*m/2 ) {
                    z = Math.sqrt(-2.0 * (Math.log(sd) + Math.log(p) + LNSQRT2PI));
                } else {
                    z = -Math.sqrt(-2.0 * (Math.log(sd) + Math.log(p) + LNSQRT2PI));
                }
            }
        }

        return new MutablePair<>(z,(twoSided ? 2.0*p : p));
    }

    /**
     * : just a shorter name for calculatePRecursively. See Mann, Whitney, [1947]
     * @param n: number of set-1 entries
     * @param m: number of set-2 entries
     * @param u: number of times a set-2 entry as preceded a set-1 entry
     * @return recursive p-value
     */
    private static double cpr(final int n, final int m, final long u) {
        if ( u < 0 ) {
            return 0.0;
        }
        if ( m == 0 || n == 0 ) {
            // there are entries in set 1 or set 2, so no set-2 entry can precede a set-1 entry; thus u must be zero.
            // note that this exists only for edification, as when we reach this point, the coefficient on this term is zero anyway
            return ( u == 0 ) ? 1.0 : 0.0;
        }


        return (((double)n)/(n+m))*cpr(n-1,m,u-m) + (((double)m)/(n+m))*cpr(n,m-1,u);
    }

    private static double cumulativeCPR(final int n, final int m, final long u ) {
        // from above:
        // the null hypothesis is that {N} is stochastically less than {M}, so U has counted
        // occurrences of {M}s before {N}s. We would expect that this should be less than (n*m+1)/2 under
        // the null hypothesis, so we want to integrate from K=0 to K=U for cumulative cases. Always.
        double p = 0.0;
        // optimization using symmetry, use the least amount of sums possible
        final long uSym = ( u <= n*m/2 ) ? u : ((long)n)*m-u;
        for ( long uu = 0; uu < uSym; uu++ ) {
            p += cpr(n,m,uu);
        }
        // correct by 1.0-p if the optimization above was used (e.g. 1-right tail = left tail)
        return (u <= n*m/2) ? p : 1.0-p;
    }

    /**
     * hook into the data tree, for testing purposes only
     * @return  observations
     */
    @VisibleForTesting
    SortedSet<Pair<Number,USet>> getObservations() {
        return observations;
    }

    /**
     * A comparator class which uses dithering on tie-breaking to ensure that the internal treeset drops no values
     * and to ensure that rank ties are broken at random.
     */
    private static class DitheringComparator implements Comparator<Pair<Number,USet>>, Serializable {
        private static final long serialVersionUID = 0L;

        @Override
        public int compare(final Pair<Number,USet> left, final Pair<Number,USet> right) {
            final double comp = Double.compare(left.getLeft().doubleValue(), right.getLeft().doubleValue());
            if ( comp > 0 ) { return 1; }
            if ( comp < 0 ) { return -1; }
            return Utils.getRandomGenerator().nextBoolean() ? -1 : 1;
        }
    }

    /**
     * A comparator that reaches into the pair and compares numbers without tie-braking.
     */
    private static class NumberedPairComparator implements Comparator<Pair<Number,USet>>, Serializable {
        private static final long serialVersionUID = 0L;

        @Override
        public int compare(final Pair<Number,USet> left, final Pair<Number,USet> right ) {
            return Double.compare(left.getLeft().doubleValue(), right.getLeft().doubleValue());
        }
    }

    public enum USet { SET1, SET2 }
    public enum ExactMode { POINT, CUMULATIVE }

}
