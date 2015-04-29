package org.broadinstitute.hellbender.utils;

import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.Serializable;
import java.util.Comparator;
import java.util.TreeSet;

public class MannWhitneyU {

    private static NormalDistribution STANDARD_NORMAL = new NormalDistribution(0.0,1.0);
    private static NormalDistribution APACHE_NORMAL = new NormalDistribution(0.0,1.0,1e-2);
    private static double LNSQRT2PI = Math.log(Math.sqrt(2.0 * Math.PI));

    private TreeSet<Pair<Number,USet>> observations;
    private int sizeSet1;
    private int sizeSet2;
    private ExactMode exactMode;

    public MannWhitneyU(ExactMode mode, boolean dither) {
        if ( dither )
            observations = new TreeSet<>(new DitheringComparator());
        else
            observations = new TreeSet<>(new NumberedPairComparator());
        sizeSet1 = 0;
        sizeSet2 = 0;
        exactMode = mode;
    }

    public MannWhitneyU() {
        this(ExactMode.POINT,true);
    }

    public MannWhitneyU(boolean dither) {
        this(ExactMode.POINT,dither);
    }

    /**
     * Add an observation into the observation tree
     * @param n: the observation (a number)
     * @param set: whether the observation comes from set 1 or set 2
     */
    public void add(Number n, USet set) {
        observations.add(new MutablePair<>(n,set));
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
    public Pair<Double,Double> runOneSidedTest(USet lessThanOther) {
        long u = calculateOneSidedU(observations, lessThanOther);
        int n = lessThanOther == USet.SET1 ? sizeSet1 : sizeSet2;
        int m = lessThanOther == USet.SET1 ? sizeSet2 : sizeSet1;
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
     * Runs the standard two-sided test,
     * returns the u-based z-approximate and p values.
     * @return a pair holding the u and p-value.
     */
    public Pair<Double,Double> runTwoSidedTest() {
        Pair<Long,USet> uPair = calculateTwoSidedU(observations);
        long u = uPair.getLeft();
        int n = uPair.getRight() == USet.SET1 ? sizeSet1 : sizeSet2;
        int m = uPair.getRight() == USet.SET1 ? sizeSet2 : sizeSet1;
        if ( n == 0 || m == 0 ) {
            // test is uninformative as one or both sets have no observations
            return new MutablePair<>(Double.NaN, Double.NaN);
        }
        return calculateP(n, m, u, true, exactMode);
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
    protected static Pair<Double,Double> calculateP(int n, int m, long u, boolean twoSided, ExactMode exactMode) {
        Pair<Double,Double> zandP;
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

    public static Pair<Double,Double> calculatePFromTable(int n, int m, long u, boolean twoSided) {
        // todo -- actually use a table for:
        // todo      - n large, m small
        return calculatePNormalApproximation(n,m,u, twoSided);
    }

    /**
     * Uses a normal approximation to the U statistic in order to return a cdf p-value. See Mann, Whitney [1947]
     * @param n - The number of entries in the stochastically smaller (dominant) set
     * @param m - The number of entries in the stochastically larger (dominated) set
     * @param u - the Mann-Whitney U value
     * @param twoSided - whether the test should be two sided
     * @return p-value associated with the normal approximation
     */
    public static Pair<Double,Double> calculatePNormalApproximation(int n,int m,long u, boolean twoSided) {
        double z = getZApprox(n,m,u);
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
    private static double getZApprox(int n, int m, long u) {
        double mean = ( ((long)m)*n+1.0)/2;
        double var = (((long) n)*m*(n+m+1.0))/12;
        double z = ( u - mean )/ Math.sqrt(var);
        return z;
    }

    /**
     * Uses a sum-of-uniform-0-1 random variable approximation to the U statistic in order to return an approximate
     * p-value. See Buckle, Kraft, van Eeden [1969] (approx) and Billingsly [1995] or Stephens, MA [1966, biometrika] (sum of uniform CDF)
     * @param n - The number of entries in the stochastically smaller (dominant) set
     * @param m - The number of entries in the stochastically larger (dominated) set
     * @param u - mann-whitney u value
     * @return p-value according to sum of uniform approx
     * todo -- this is currently not called due to not having a good characterization of where it is significantly more accurate than the
     * todo -- normal approxmation (e.g. enough to merit the runtime hit)
     */
    public static double calculatePUniformApproximation(int n, int m, long u) {
        long R = u + (n*(n+1))/2;
        double a = Math.sqrt(m * (n + m + 1));
        double b = (n/2.0)*(1- Math.sqrt((n + m + 1) / m));
        double z = b + ((double)R)/a;
        if ( z < 0 ) { return 1.0; }
        else if ( z > n ) { return 0.0; }
        else {
            if ( z > ((double) n) /2 ) {
                return 1.0-1/(CombinatoricsUtils.factorialDouble(n))*uniformSumHelper(z, (int) Math.floor(z), n, 0);
            } else {
                return 1/(CombinatoricsUtils.factorialDouble(n))*uniformSumHelper(z, (int) Math.floor(z), n, 0);
            }
        }
    }

    /**
     * Helper function for the sum of n uniform random variables
     * @param z - value at which to compute the (un-normalized) cdf
     * @param m - a cutoff integer (defined by m <= z < m + 1)
     * @param n - the number of uniform random variables
     * @param k - holder variable for the recursion (alternatively, the index of the term in the sequence)
     * @return the (un-normalized) cdf for the sum of n random variables
     */
    private static double uniformSumHelper(double z, int m, int n, int k) {
        if ( k > m ) { return 0; }
        int coef = (k % 2 == 0) ? 1 : -1;
        return coef*CombinatoricsUtils.binomialCoefficient(n, k)* Math.pow(z - k, n) + uniformSumHelper(z,m,n,k+1);
    }

    /**
     * Calculates the U-statistic associated with a two-sided test (e.g. the RV from which one set is drawn
     * stochastically dominates the RV from which the other set is drawn); two-sidedness is accounted for
     * later on simply by multiplying the p-value by 2.
     *
     * Recall: If X stochastically dominates Y, the test is for occurrences of Y before X, so the lower value of u is chosen
     * @param observed - the observed data
     * @return the minimum of the U counts (set1 dominates 2, set 2 dominates 1)
     */
    public static Pair<Long,USet> calculateTwoSidedU(TreeSet<Pair<Number,USet>> observed) {
        int set1SeenSoFar = 0;
        int set2SeenSoFar = 0;
        long uSet1DomSet2 = 0;
        long uSet2DomSet1 = 0;
        USet previous = null;
        for ( Pair<Number,USet> dataPoint : observed ) {

            if ( dataPoint.getRight() == USet.SET1 ) {
                ++set1SeenSoFar;
            } else {
                ++set2SeenSoFar;
            }

            if ( previous != null ) {
                if ( dataPoint.getRight() == USet.SET1 ) {
                    uSet2DomSet1 += set2SeenSoFar;
                } else {
                    uSet1DomSet2 += set1SeenSoFar;
                }
            }

            previous = dataPoint.getRight();
        }

        return uSet1DomSet2 < uSet2DomSet1 ? new MutablePair<>(uSet1DomSet2, USet.SET1) : new MutablePair<>(uSet2DomSet1, USet.SET2);
    }

    /**
     * Calculates the U-statistic associated with the one-sided hypothesis that "dominator" stochastically dominates
     * the other U-set. Note that if S1 dominates S2, we want to count the occurrences of points in S2 coming before points in S1.
     * @param observed - the observed data points, tagged by each set
     * @param dominator - the set that is hypothesized to be stochastically dominating
     * @return the u-statistic associated with the hypothesis that dominator stochastically dominates the other set
     */
    public static long calculateOneSidedU(TreeSet<Pair<Number,USet>> observed,USet dominator) {
        long otherBeforeDominator = 0l;
        int otherSeenSoFar = 0;
        for ( Pair<Number,USet> dataPoint : observed ) {
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
    public static Pair<Double,Double> calculatePRecursively(int n, int m, long u, boolean twoSided, ExactMode mode) {
        if ( m > 8 && n > 5 ) { throw new GATKException(String.format("Please use the appropriate (normal or sum of uniform) approximation. Values n: %d, m: %d", n, m)); }
        double p = mode == ExactMode.POINT ? cpr(n,m,u) : cumulativeCPR(n,m,u);
        //p *= twoSided ? 2.0 : 1.0;
        double z;
        if ( mode == ExactMode.CUMULATIVE ) {
            z = APACHE_NORMAL.inverseCumulativeProbability(p);
        } else {
            double sd = Math.sqrt((1.0 + 1.0 / (1 + n + m)) * (n * m) * (1.0 + n + m) / 12); // biased variance empirically better fit to distribution then asymptotic variance
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
     * Hook into CPR with sufficient warning (for testing purposes)
     * calls into that recursive calculation.
     * @param n: number of set-one entries (hypothesis: set one is stochastically less than set two)
     * @param m: number of set-two entries
     * @param u: number of set-two entries that precede set-one entries (e.g. 0,1,0,1,0 -> 3 )
     * @return same as cpr
     */
    protected static double calculatePRecursivelyDoNotCheckValuesEvenThoughItIsSlow(int n, int m, long u) {
        return cpr(n,m,u);
    }

    /**
     * : just a shorter name for calculatePRecursively. See Mann, Whitney, [1947]
     * @param n: number of set-1 entries
     * @param m: number of set-2 entries
     * @param u: number of times a set-2 entry as preceded a set-1 entry
     * @return recursive p-value
     */
    private static double cpr(int n, int m, long u) {
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

    private static double cumulativeCPR(int n, int m, long u ) {
        // from above:
        // the null hypothesis is that {N} is stochastically less than {M}, so U has counted
        // occurrences of {M}s before {N}s. We would expect that this should be less than (n*m+1)/2 under
        // the null hypothesis, so we want to integrate from K=0 to K=U for cumulative cases. Always.
        double p = 0.0;
        // optimization using symmetry, use the least amount of sums possible
        long uSym = ( u <= n*m/2 ) ? u : ((long)n)*m-u;
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
    protected TreeSet<Pair<Number,USet>> getObservations() {
        return observations;
    }

    /**
     * hook into the set sizes, for testing purposes only
     * @return size set 1, size set 2
     */
    protected Pair<Integer,Integer> getSetSizes() {
        return new MutablePair<>(sizeSet1,sizeSet2);
    }

    /**
     * A comparator class which uses dithering on tie-breaking to ensure that the internal treeset drops no values
     * and to ensure that rank ties are broken at random.
     */
    private static class DitheringComparator implements Comparator<Pair<Number,USet>>, Serializable {
        private static final long serialVersionUID = 0L;

        @Override
        public int compare(Pair<Number,USet> left, Pair<Number,USet> right) {
            double comp = Double.compare(left.getLeft().doubleValue(), right.getLeft().doubleValue());
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
        public int compare(Pair<Number,USet> left, Pair<Number,USet> right ) {
            return Double.compare(left.getLeft().doubleValue(), right.getLeft().doubleValue());
        }
    }

    public enum USet { SET1, SET2 }
    public enum ExactMode { POINT, CUMULATIVE }

}
