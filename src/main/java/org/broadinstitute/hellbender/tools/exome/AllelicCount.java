package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;
import java.util.TreeMap;

/**
 * Reference and alternate allele counts at a SNP site specified by an interval.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicCount implements Locatable {
    private final SimpleInterval interval;
    private final int refReadCount, altReadCount;

    // constants for Brent optimization of minor allele fraction
    private static final MaxEval BRENT_MAX_EVAL = new MaxEval(100);
    private static final double MINOR_ALLELE_FRACTION_RELATIVE_TOLERANCE = 0.00001;
    private static final double MINOR_ALLELE_FRACTION_ABSOLUTE_TOLERANCE = 0.00001;
    private static final BrentOptimizer OPTIMIZER =
            new BrentOptimizer(MINOR_ALLELE_FRACTION_RELATIVE_TOLERANCE, MINOR_ALLELE_FRACTION_ABSOLUTE_TOLERANCE);

    public AllelicCount(final SimpleInterval interval, final int refReadCount, final int altReadCount) {
        Utils.nonNull(interval, "Can't construct AllelicCount with null interval.");
        if (refReadCount < 0 || altReadCount < 0) {
            throw new IllegalArgumentException("Can't construct AllelicCount with negative read counts.");
        }
        if (refReadCount + altReadCount == 0) {
            throw new IllegalArgumentException("Can't construct AllelicCount with zero total counts.");
        }
        this.interval = interval;
        this.refReadCount = refReadCount;
        this.altReadCount = altReadCount;
    }

    @Override
    public String getContig() {return interval.getContig(); }

    @Override
    public int getStart() {return interval.getStart(); }

    @Override
    public int getEnd() {return interval.getEnd(); }

    public Interval getInterval() {
        return new Interval(interval.getContig(), interval.getStart(), interval.getEnd());
    }

    public int getRefReadCount() {  return refReadCount;        }

    public int getAltReadCount() {  return altReadCount;        }

    /**
     * Returns the maximum likelihood estimate of the alternate-allele fraction.
     * @return      alternate-allele fraction
     */
    public double estimateAltAlleleFraction() {
        return (double) altReadCount / (refReadCount + altReadCount);
    }

    /**
     * Returns the maximum likelihood estimate of the minor-allele fraction.
     *
     * Ignoring allelic bias the likelihood of a alt reads and r ref reads is proportional to f^a(1-f)^r + f^r(1-f)^a,
     * where f in [0,1/2] is the minor allele fraction and NOT the alt allele fraction.  The two terms derive from the
     * a priori equally likely cases that the alt and ref alleles are the minor allele, respectively.
     * This likelihood is unimodal on [0,1/2] so numerical max-finding is straightforward.
     *
     * @return      maximum likelihood estimate of the minor allele fraction
     */
    public double estimateMinorAlleleFraction() {
        return MinorAlleleFractionCache.get(altReadCount, refReadCount);
    }

    /**
     * Returns a TargetCoverage with coverage given by minor allele fraction (plus epsilon, to avoid issues caused by
     * taking logs).
     * @param name  target name
     * @return      TargetCoverage with coverage given by minor allele fraction
     */
    public TargetCoverage toMinorAlleleFractionTargetCoverage(final String name) {
        return new TargetCoverage(name, new SimpleInterval(interval), estimateMinorAlleleFraction() + Double.MIN_VALUE);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof AllelicCount)) {
            return false;
        }

        final AllelicCount count = (AllelicCount) o;
        return interval.equals(count.interval) && refReadCount == count.refReadCount
                && altReadCount == count.altReadCount;
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + refReadCount;
        result = 31 * result + altReadCount;
        return result;
    }

    /**
     * A helper class to maintain a cache of computed maximum likelihood estimates of minor allele fraction
     */
    private static final class MinorAlleleFractionCache {

        private static final Map<Pair<Integer, Integer>, Double> cache = new TreeMap<>();

        /**
         * get cached value or cache a new value as necessary
         * @param a alt read count
         * @param r ref read count
         * @return maximum likelihood estimate of minor allele fraction
         */
        public static double get(final int a, final int r) {
            final Pair<Integer, Integer> pair = new ImmutablePair<>(Math.min(a, r), Math.max(a, r));    //by symmetry
            if (cache.containsKey(pair)) {
                return cache.get(pair);
            } else {
                final double minorAlleleFraction = calculateMinorAlleleFraction(a, r);
                cache.put(pair, minorAlleleFraction);
                return minorAlleleFraction;
            }
        }

        private static double calculateMinorAlleleFraction(final int a, final int r) {
            final double altFraction = (double) a / (a + r);
            final double initialEstimate = altFraction < 0.5 ? altFraction : 1 - altFraction;
            final SearchInterval searchInterval = new SearchInterval(0.0, 0.5, initialEstimate);

            // work in log space to avoid underflow
            final UnivariateObjectiveFunction objective = new UnivariateObjectiveFunction(f -> {
                final double logf = Math.log(f);
                final double logOneMinusf = Math.log(1-f);
                final double altMinorLogLikelihood = a * logf + r * logOneMinusf;
                final double refMinorLogLikelihood = r * logf + a * logOneMinusf;
                return GATKProtectedMathUtils.naturalLogSumExp(altMinorLogLikelihood, refMinorLogLikelihood);
            });

            return OPTIMIZER.optimize(objective, GoalType.MAXIMIZE, searchInterval, BRENT_MAX_EVAL).getPoint();
        }
    }
}
