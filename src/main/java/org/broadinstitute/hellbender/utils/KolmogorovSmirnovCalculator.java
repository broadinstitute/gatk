package org.broadinstitute.hellbender.utils;

import java.util.Arrays;

/**
 * Computes a two-sample Kolmogorov-Smirnov statistic, and tests it for statistical significance at a specified level.
 * Build this with a reference histogram of counts and a significance level,
 * and then test other histograms to see whether they're likely to be drawn from the same population or not.
 */
public class KolmogorovSmirnovCalculator {
    private final double cOfAlpha;
    private final long nCounts;
    private final float[] cdf;

    public KolmogorovSmirnovCalculator(final long[] refHisto, final float alpha ) {
        if ( alpha <= 0.f || alpha >= 1. ) {
            throw new IllegalArgumentException("Alpha should be a small positive number, like maybe .05 -- "+
                    alpha+" is not appropriate.");
        }
        cOfAlpha = Math.sqrt(-.5*Math.log(alpha/2.));
        nCounts = Arrays.stream(refHisto).sum();
        if ( nCounts < 1 ) {
            throw new IllegalArgumentException("Empty reference histogram.");
        }
        cdf = new float[refHisto.length];
        long sum = 0L;
        for ( int idx = 0; idx != refHisto.length; ++idx ) {
            final long val = refHisto[idx];
            if ( val < 0 ) {
                throw new IllegalArgumentException("All values in the reference histogram must be non-negative.");
            }
            sum += refHisto[idx];
            cdf[idx] = (float)sum/nCounts;
        }
    }

    public boolean isDifferent( final long[] sampleHisto ) {
        if ( sampleHisto.length != cdf.length ) {
            throw new IllegalArgumentException("Sample histogram and reference histogram have different lengths.");
        }
        // calculate mCounts, the sum of counts in the sample histogram
        long mCounts = 0L;
        for ( int idx = 0; idx != sampleHisto.length; ++idx ) {
            final long val = sampleHisto[idx];
            if ( val < 0 ) {
                throw new IllegalArgumentException("All values in the sample histogram must be non-negative.");
            }
            mCounts += val;
        }
        if ( mCounts < 1 ) {
            throw new IllegalArgumentException("Empty sample histogram.");
        }

        // K-S statistic for the given alpha
        final float significantDiff = (float)(cOfAlpha*Math.sqrt((double)(nCounts +mCounts)/nCounts/mCounts));

        long sum = 0L; // sum from 0 to idx of sample histogram counts
        int idx = 0; // current index into sample histogram
        float curVal = 0.f; // sum/mCounts as a float (between 0 and 1)
        boolean prevZero = false; // whether the value at the previous index was 0 or not
        while ( sum < mCounts ) { // i.e., until we've processed all non-zero counts in the sample histogram
            final long val = sampleHisto[idx];
            if ( val == 0 ) prevZero = true;
            else {
                // an extremum can occur at the end of a run of zero counts
                if ( prevZero && Math.abs(cdf[idx-1] - curVal) >= significantDiff ) return true;
                prevZero = false;
                sum += val;
                curVal = (float)sum/mCounts;
                // an extremum can occur at any non-zero count
                if ( Math.abs(cdf[idx] - curVal) >= significantDiff ) return true;
            }
            idx += 1; // this could run off the end, but the while condition will fail before we attempt to use it
        }
        return false;
    }
}
