package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.IntegerDistribution;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.Random;

/** Histogram of observations on a compact set of non-negative integer values. */
@DefaultSerializer(IntHistogram.Serializer.class)
public final class IntHistogram {
    // the final slot (i.e., counts[counts.length-1]) accumulates the count of all observations greater than
    // the maximum tracked value -- it's the "overflow" bin
    private final long[] counts;
    private long totalObservations;

    public IntHistogram( final int maxTrackedValue ) {
        Utils.validateArg(maxTrackedValue >= 0,"maxTrackedValue must be non-negative");
        counts = new long[maxTrackedValue + 2];
        totalObservations = 0;
    }

    private IntHistogram( final long[] counts, final long totalObservations) {
        this.counts = counts;
        this.totalObservations = totalObservations;
    }

    private IntHistogram( final Kryo kryo, final Input input ) {
        final int len = input.readInt();
        counts = new long[len];
        long total = 0L;
        for ( int idx = 0; idx != len; ++idx ) {
            final long val = input.readLong();
            final double idxSum = val * idx;
            total += val;
            counts[idx] = val;
        }
        totalObservations = total;
    }

    private void serialize( final Kryo kryo, final Output output ) {
        output.writeInt(counts.length);
        for ( final long val : counts ) {
            output.writeLong(val);
        }
    }

    public void addObservation( final int observedValue ) {
        Utils.validateArg(observedValue >= 0,"observedValue must be non-negative");
        counts[observedValue >= counts.length ? counts.length-1 : observedValue] += 1;
        totalObservations += 1;
    }

    public void addObservations( final int observedValue, final long nObservations ) {
        Utils.validateArg(observedValue >= 0,"observedValue must be non-negative");
        Utils.validateArg(nObservations >= 0,"nObservations must be non-negative");
        counts[observedValue >= counts.length ? counts.length-1 : observedValue] += nObservations;
        totalObservations += nObservations;
    }

    public void addObservations( final IntHistogram intHistogram ) {
        final long[] thoseCounts = intHistogram.counts;
        Utils.validateArg(counts.length == thoseCounts.length,
                "The supplied histogram doesn't have the right shape.");
        for ( int idx = 0; idx != counts.length; ++idx ) {
            counts[idx] += thoseCounts[idx];
        }
        totalObservations += intHistogram.totalObservations;
    }

    public int getMaximumTrackedValue() { return counts.length - 2; }

    public long getTotalObservations() { return totalObservations; }

    public long getNObservations( final int observedValue ) {
        Utils.validateArg(observedValue >= 0 && observedValue < counts.length,
                "observedValue must be non-negative, and no more than 1 greater than the maximum tracked observed value");
        return counts[observedValue];
    }

    public void clear() {
        Arrays.fill(counts, 0L);
        totalObservations = 0L;
    }

    /** Returns a histogram with the smallest tracked maximum such that the fraction of untracked observations does not
     *  exceed 1/trackedMaximum.  If the fraction of untracked observations in the histogram already exceeds
     *  1/trackedMaximum without any further trimming, then this histogram is returned rather than a new one.
     */
    public IntHistogram trim() {
        long nUntrackedObservations = counts[counts.length - 1];
        int nBins = counts.length;

        if ( nBins * nUntrackedObservations >= totalObservations ) {
            return this;
        }

        do {
            nUntrackedObservations += counts[nBins - 2];
            nBins -= 1;
        } while ( nBins * nUntrackedObservations < totalObservations );

        nBins += 1;
        nUntrackedObservations -= counts[nBins - 2];
        final long[] newCounts = Arrays.copyOfRange(counts, 0, nBins);
        newCounts[nBins - 1] = nUntrackedObservations;
        return new IntHistogram(newCounts, totalObservations);
    }

    public CDF getCDF() { return new CDF(this); }

    public final static class Serializer extends com.esotericsoftware.kryo.Serializer<IntHistogram> {
        @Override
        public void write( final Kryo kryo, final Output output, final IntHistogram histogram ) {
            histogram.serialize(kryo, output);
        }

        @Override
        public IntHistogram read( final Kryo kryo, final Input input, final Class<IntHistogram> klass ) {
            return new IntHistogram(kryo, input);
        }
    }

    /**
     * Composes an empirical distribution based on the current contents fo the histogram.
     * <p>
     *     Later changes in the histogram won't affect the returned distribution.
     * </p>
     * <p>
     *     You can indicate an arbitrary "smoothing" count number which added to the
     *     observed frequency of every value. This way non observed values won't necessarily
     *     have a probability of zero.
     * </p>
     * <p>
     *     The supported space is enclosed in <code>[0,{@link Integer#MAX_VALUE max_int}]</code>.
     *     So the probalility of negative values is 0 and the probability of 0 and positive values is
     *     always at least <code>smoothing / num. of observations</code>.
     * </p>
     * <p>
     *     Despite that the probability of very large values not tracked by the histogram is not zero,
     *     the cumulative probability is said to reach 1.0 at the largest tracked number.
     *     As a consequence the distribution is actually not a proper one (unknown normalization constant).
     * </p>
     * <p>
     *     However since we expect that just a small fraction of the observations will fall outside the tracked range
     *     we assume that is a proper distribution as far as the calculation of the mean, variance, density and
     *     cumulative density is concern.
     * </p>
     *
     * @param smoothing zero or a positive number of counts to each value.
     * @return never {@code null}.
     */
    public AbstractIntegerDistribution empiricalDistribution(final int smoothing) {
        ParamUtils.isPositiveOrZero(smoothing, "the smoothing must be zero or positive");
        final long[] counts = Arrays.copyOfRange(this.counts, 0, this.counts.length - 1);
        final double[] cumulativeCounts = new double[counts.length];
        double sum = 0;
        double sqSum = 0;
        for (int i = 0; i < counts.length; i++) {
            final long newCount = (counts[i] += smoothing);
            sum += newCount * i;
            sqSum += i * newCount * i;
        }
        cumulativeCounts[0] = counts[0];
        for (int i = 1; i < counts.length; i++) {
            cumulativeCounts[i] = counts[i] + cumulativeCounts[i - 1];
        }
        final double totalCounts = cumulativeCounts[counts.length - 1];
        final double inverseTotalCounts = 1.0 / totalCounts;
        final double mean = sum / totalCounts;
        final double variance = sqSum / totalCounts - mean * mean;
        final int seed = Arrays.hashCode(counts);
        final RandomGenerator rdnGen = new JDKRandomGenerator();
        rdnGen.setSeed(seed);
        return new AbstractIntegerDistribution(rdnGen) {

            private static final long serialVersionUID = -1L;

            @Override
            public double probability(int x) {
                if (x < 0) {
                    return 0.0;
                } else if (x >= counts.length) {
                    return smoothing * inverseTotalCounts;
                } else {
                    return counts[x] * inverseTotalCounts;
                }
            }

            @Override
            public double cumulativeProbability(int x) {
                if (x < 0) {
                    return 0;
                } else if (x >= counts.length) {
                    return 1.0;
                } else {
                    return cumulativeCounts[x] * inverseTotalCounts;
                }
            }

            @Override
            public double getNumericalMean() {
                return mean;
            }

            @Override
            public double getNumericalVariance() {
                return variance;
            }

            @Override
            public int getSupportLowerBound() {
                return 0;
            }

            @Override
            public int getSupportUpperBound() {
                return Integer.MAX_VALUE;
            }

            @Override
            public boolean isSupportConnected() {
                return true;
            }
        };
    }

    @DefaultSerializer(CDF.Serializer.class)
    public final static class CDF {
        final float[] cdfFractions;
        final long nCounts;

        public CDF( final IntHistogram intHistogram ) {
            final int arrayLen = intHistogram.getMaximumTrackedValue() + 2;
            cdfFractions = new float[arrayLen];
            this.nCounts = intHistogram.getTotalObservations();
            long sum = 0;
            for ( int observedValue = 0; observedValue != arrayLen; ++observedValue ) {
                sum += intHistogram.getNObservations(observedValue);
                cdfFractions[observedValue] = (float)sum / nCounts;
            }
        }

        public CDF(final float[] cdfFractions, final long nCounts) {
            this.cdfFractions = cdfFractions;
            this.nCounts = nCounts;
        }

        private CDF( final Kryo kryo, final Input input ) {
            final int len = input.readInt();
            cdfFractions = new float[len];
            for ( int idx = 0; idx != len; ++idx ) {
                cdfFractions[idx] = input.readFloat();
            }
            nCounts = input.readLong();
        }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeInt(cdfFractions.length);
            for ( final float val : cdfFractions ) {
                output.writeFloat(val);
            }
            output.writeLong(nCounts);
        }

        public int size() { return cdfFractions.length; }
        public float getFraction( final int idx ) { return cdfFractions[idx]; }
        public long getTotalObservations() { return nCounts; }

        public IntHistogram createEmptyHistogram() {
            return new IntHistogram(cdfFractions.length - 2);
        }

        /** Using the Kolmogorov-Smirnov statistic for two samples:
         *  Is the specified histogram significantly different from the CDF?
         *  This is what we use in production code: It minimizes the number of comparisons, and quits
         *  as soon as it finds a significant K-S stat (rather than finding the max). */
        public boolean isDifferentByKSStatistic( final IntHistogram sampleHistogram, final float significance ) {
            final long[] sampleCounts = sampleHistogram.counts;

            Utils.validateArg(significance > 0f && significance < .2f,
                    "The significance must be specifed as a probability of a chance occurrence (a number like .01f).");
            Utils.validateArg(sampleCounts.length == cdfFractions.length,
                    "The supplied histogram doesn't have the right size.");

            final long mCounts = sampleHistogram.getTotalObservations();
            Utils.validateArg(mCounts > 0, "The supplied histogram is empty.");

            final double cOfAlpha = Math.sqrt(-.5 * Math.log(significance / 2.));
            // minimum K-S statistic for the given significance
            final float significantDiff = (float)(cOfAlpha * Math.sqrt((double)(nCounts + mCounts) / nCounts / mCounts));

            long sum = 0L; // sum from 0 to idx of sample histogram counts
            int idx = 0; // current index into sample histogram
            float curVal = 0.f; // sum/mCounts as a float (between 0 and 1)
            boolean prevZero = false; // whether the value at the previous index was 0 or not
            final long trackedTotal = mCounts - sampleCounts[sampleCounts.length - 1];
            while ( sum < trackedTotal ) { // i.e., until we've processed all non-zero counts in the sample histogram
                final long val = sampleCounts[idx];
                if ( val == 0 ) prevZero = true;
                else {
                    // an extremum can occur at the end of a run of zero counts
                    if ( prevZero && Math.abs(cdfFractions[idx - 1] - curVal) >= significantDiff ) return true;
                    prevZero = false;
                    sum += val;
                    curVal = (float)sum / mCounts;
                    // an extremum can occur at any non-zero count
                    if ( Math.abs(cdfFractions[idx] - curVal) >= significantDiff ) return true;
                }
                idx += 1;
            }
            return false;
        }

        public int median() { return internalPopStat(.5f); }
        public int leftMedianDeviation( final int median ) { return median - internalPopStat(.25f); }
        public int rightMedianDeviation( final int median ) {
            return internalPopStat(.75f) - median;
        }

        /** An observed value for which the specified fraction of the sample population has that value or less. */
        public int popStat( final float popFraction ) {
            Utils.validateArg(popFraction >= 0.f && popFraction <= 1.f, "popFraction must be between 0 and 1");
            return internalPopStat(popFraction);
        }

        private int internalPopStat( final float popFraction ) {
            final int idx = Arrays.binarySearch(cdfFractions, popFraction);
            return idx < 0 ? -idx - 1 : idx;
        }

        public final static class Serializer extends com.esotericsoftware.kryo.Serializer<CDF> {
            @Override
            public void write( final Kryo kryo, final Output output, final CDF cdf ) {
                cdf.serialize(kryo, output);
            }

            @Override
            public CDF read( final Kryo kryo, final Input input, final Class<CDF> klass ) {
                return new CDF(kryo, input);
            }
        }
    }
}
