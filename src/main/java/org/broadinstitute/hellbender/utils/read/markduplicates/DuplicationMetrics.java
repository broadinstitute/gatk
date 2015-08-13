package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.util.Histogram;

import java.io.Serializable;

/**
 * Metrics that are calculated during the process of marking duplicates
 * within a stream of SAMRecords.
 */
@SuppressWarnings("serial")
public final class DuplicationMetrics extends MetricBase implements Serializable {
    /** The library on which the duplicate marking was performed. */
    public String LIBRARY;

    /**
     * The number of mapped reads examined which did not have a mapped mate pair,
     * either because the read is unpaired, or the read is paired to an unmapped mate.
     */
    public long UNPAIRED_READS_EXAMINED;

    /** The number of mapped read pairs examined. */
    public long READ_PAIRS_EXAMINED;

    /** The total number of unmapped reads examined. */
    public long UNMAPPED_READS;

    /** The number of fragments that were marked as duplicates. */
    public long UNPAIRED_READ_DUPLICATES;

    /** The number of read pairs that were marked as duplicates. */
    public long READ_PAIR_DUPLICATES;

    /**
     * The number of read pairs duplicates that were caused by optical duplication.
     * Value is always < READ_PAIR_DUPLICATES, which counts all duplicates regardless of source.
     */
    public long READ_PAIR_OPTICAL_DUPLICATES;

    /** The percentage of mapped sequence that is marked as duplicate. */
    public Double PERCENT_DUPLICATION;

    /** The estimated number of unique molecules in the library based on PE duplication. */
    public Long ESTIMATED_LIBRARY_SIZE;

    /**
     * Fills in the ESTIMATED_LIBRARY_SIZE based on the paired read data examined where
     * possible and the PERCENT_DUPLICATION.
     */
    public void calculateDerivedMetrics() {
        this.ESTIMATED_LIBRARY_SIZE = estimateLibrarySize(this.READ_PAIRS_EXAMINED - this.READ_PAIR_OPTICAL_DUPLICATES,
                this.READ_PAIRS_EXAMINED - this.READ_PAIR_DUPLICATES);

        PERCENT_DUPLICATION = (UNPAIRED_READ_DUPLICATES + READ_PAIR_DUPLICATES *2) /(double) (UNPAIRED_READS_EXAMINED + READ_PAIRS_EXAMINED *2);
    }

    /**
     * Estimates the size of a library based on the number of paired end molecules observed
     * and the number of unique pairs ovserved.
     *
     * Based on the Lander-Waterman equation that states:
     *     C/X = 1 - exp( -N/X )
     * where
     *     X = number of distinct molecules in library
     *     N = number of read pairs
     *     C = number of distinct fragments observed in read pairs
     */
    public static Long estimateLibrarySize(final long readPairs, final long uniqueReadPairs) {
        final long readPairDuplicates = readPairs - uniqueReadPairs;

        if (readPairs > 0 && readPairDuplicates > 0) {
            long n = readPairs;
            long c = uniqueReadPairs;

            double m = 1.0, M = 100.0;

            if (c >= n || f(m*c, c, n) < 0) {
                throw new IllegalStateException("Invalid values for pairs and unique pairs: "
                        + n + ", " + c);

            }

            while( f(M*c, c, n) >= 0 ) M *= 10.0;

            for (int i=0; i<40; i++ ) {
                double r = (m+M)/2.0;
                double u = f( r * c, c, n );
                if ( u == 0 ) break;
                else if ( u > 0 ) m = r;
                else if ( u < 0 ) M = r;
            }

            return (long) (c * (m+M)/2.0);
        }
        else {
            return null;
        }
    }

    /** Method that is used in the computation of estimated library size. */
    private static double f(double x, double c, double n) {
        return c/x - 1 + Math.exp(-n/x);
    }

    /**
     * Estimates the ROI (return on investment) that one would see if a library was sequenced to
     * x higher coverage than the observed coverage.
     *
     * @param estimatedLibrarySize the estimated number of molecules in the library
     * @param x the multiple of sequencing to be simulated (i.e. how many X sequencing)
     * @param pairs the number of pairs observed in the actual sequencing
     * @param uniquePairs the number of unique pairs observed in the actual sequencing
     * @return a number z <= x that estimates if you had pairs*x as your sequencing then you
     *         would observe uniquePairs*z unique pairs.
     */
    public static double estimateRoi(long estimatedLibrarySize, double x, long pairs, long uniquePairs) {
        return estimatedLibrarySize * ( 1 - Math.exp(-(x*pairs)/estimatedLibrarySize) ) / uniquePairs;
    }

    /**
     * Calculates a histogram using the estimateRoi method to estimate the effective yield
     * doing x sequencing for x=1..10.
     */
    public Histogram<Double> calculateRoiHistogram() {
        if (ESTIMATED_LIBRARY_SIZE == null) {
            try {
                calculateDerivedMetrics();
                if (ESTIMATED_LIBRARY_SIZE == null) return null;
            }
            catch (IllegalStateException ise) { return null; }
        }

        long uniquePairs = READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES;
        Histogram<Double> histo = new Histogram<>();

        for (double x=1; x<=100; x+=1) {
            histo.increment(x, estimateRoi(ESTIMATED_LIBRARY_SIZE, x, READ_PAIRS_EXAMINED, uniquePairs));
        }

        return histo;
    }

    public DuplicationMetrics copy() {
        final DuplicationMetrics copy = new DuplicationMetrics();
        copy.LIBRARY = this.LIBRARY;
        copy.UNPAIRED_READS_EXAMINED = this.UNPAIRED_READS_EXAMINED;
        copy.READ_PAIRS_EXAMINED = this.READ_PAIRS_EXAMINED;
        copy.UNMAPPED_READS = this.UNMAPPED_READS;
        copy.UNPAIRED_READ_DUPLICATES = this.UNPAIRED_READ_DUPLICATES;
        copy.READ_PAIR_DUPLICATES = this.READ_PAIR_DUPLICATES;
        copy.READ_PAIR_OPTICAL_DUPLICATES = this.READ_PAIR_OPTICAL_DUPLICATES;
        copy.PERCENT_DUPLICATION = this.PERCENT_DUPLICATION;
        copy.ESTIMATED_LIBRARY_SIZE = this.ESTIMATED_LIBRARY_SIZE;

        return copy;
    }

}
