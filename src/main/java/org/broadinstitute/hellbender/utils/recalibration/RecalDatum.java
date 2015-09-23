package org.broadinstitute.hellbender.utils.recalibration;

import htsjdk.samtools.SAMUtils;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.io.Serializable;

/**
 * An individual piece of recalibration data. Each bin counts up the number of observations and the number
 * of reference mismatches seen for that combination of covariates.
 */
public final class RecalDatum implements Serializable {
    public static final byte MAX_RECALIBRATED_Q_SCORE = SAMUtils.MAX_PHRED_SCORE;
    private static final double UNINITIALIZED = -1.0;
    private static final long serialVersionUID = 1L;
    private static final double MULTIPLIER = 1000.0;

    /**
     * estimated reported quality score based on combined data's individual q-reporteds and number of observations
     */
    private double estimatedQReported;

    /**
     * the empirical quality for datums that have been collapsed together (by read group and reported quality, for example)
     */
    private double empiricalQuality;

    /**
     * number of bases seen in total
     */
    private long numObservations;

    /**
     * number of bases seen that didn't match the reference
     * (actually sum of the error weights - so not necessarily a whole number)
     * Stored with an internal multiplier to keep it closer to the floating-point sweet spot and avoid numerical error
     * (see https://github.com/broadinstitute/hellbender/wiki/Numerical-errors ).
     */
    private double numMismatches;

    /**
     * used when calculating empirical qualities to avoid division by zero
     */
    private static final int SMOOTHING_CONSTANT = 1;

    //---------------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Create a new RecalDatum with given observation and mismatch counts, and an reported quality
     *
     * @param _numObservations    observations
     * @param _numMismatches      mismatches
     * @param reportedQuality     Qreported
     */
    public RecalDatum(final long _numObservations, final double _numMismatches, final byte reportedQuality) {
        if ( _numObservations < 0 ) throw new IllegalArgumentException("numObservations < 0");
        if ( _numMismatches < 0.0 ) throw new IllegalArgumentException("numMismatches < 0");
        if ( reportedQuality < 0 ) throw new IllegalArgumentException("reportedQuality < 0");

        numObservations = _numObservations;
        numMismatches = (_numMismatches*MULTIPLIER);
        estimatedQReported = reportedQuality;
        empiricalQuality = UNINITIALIZED;
    }

    /**
     * Copy copy into this recal datum, overwriting all of this objects data
     * @param copy  RecalDatum to copy
     */
    public RecalDatum(final RecalDatum copy) {
        this.numObservations = copy.numObservations;
        this.numMismatches = copy.numMismatches;
        this.estimatedQReported = copy.estimatedQReported;
        this.empiricalQuality = copy.empiricalQuality;
    }

    /**
     * Add in all of the data from other into this object, updating the reported quality from the expected
     * error rate implied by the two reported qualities
     *
     * @param other  RecalDatum to combine
     */
    public void combine(final RecalDatum other) {
        final double sumErrors = this.calcExpectedErrors() + other.calcExpectedErrors();
        increment(other.getNumObservations(), other.getNumMismatches());
        estimatedQReported = -10 * Math.log10(sumErrors / getNumObservations());
        empiricalQuality = UNINITIALIZED;
    }

    public void setEstimatedQReported(final double estimatedQReported) {
        if ( estimatedQReported < 0 ) throw new IllegalArgumentException("estimatedQReported < 0");
        if ( Double.isInfinite(estimatedQReported) ) throw new IllegalArgumentException("estimatedQReported is infinite");
        if ( Double.isNaN(estimatedQReported) ) throw new IllegalArgumentException("estimatedQReported is NaN");

        this.estimatedQReported = estimatedQReported;
        empiricalQuality = UNINITIALIZED;
    }

    public final double getEstimatedQReported() {
        return estimatedQReported;
    }
    public final byte getEstimatedQReportedAsByte() {
        return (byte)(int)(Math.round(getEstimatedQReported()));
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // Empirical quality score -- derived from the num mismatches and observations
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Returns the error rate (in real space) of this interval, or 0 if there are no observations
     * @return the empirical error rate ~= N errors / N obs
     */
    public double getEmpiricalErrorRate() {
        if ( numObservations == 0 )
            return 0.0;
        else {
            // cache the value so we don't call log over and over again
            final double doubleMismatches = (numMismatches/MULTIPLIER) + SMOOTHING_CONSTANT;
            // smoothing is one error and one non-error observation, for example
            final double doubleObservations = numObservations + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT;
            return doubleMismatches / doubleObservations;
        }
    }

    public void setEmpiricalQuality(final double empiricalQuality) {
        if ( empiricalQuality < 0 ) throw new IllegalArgumentException("empiricalQuality < 0");
        if ( Double.isInfinite(empiricalQuality) ) throw new IllegalArgumentException("empiricalQuality is infinite");
        if ( Double.isNaN(empiricalQuality) ) throw new IllegalArgumentException("empiricalQuality is NaN");

        this.empiricalQuality = empiricalQuality;
    }

    public final double getEmpiricalQuality() {
        return getEmpiricalQuality(getEstimatedQReported());
    }

    public final double getEmpiricalQuality(final double conditionalPrior) {
        if (empiricalQuality == UNINITIALIZED) {
            calcEmpiricalQuality(conditionalPrior);
        }
        return empiricalQuality;
    }

    public final byte getEmpiricalQualityAsByte() {
        return (byte)(Math.round(getEmpiricalQuality()));
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // toString methods
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public String toString() {
        return String.format("%d,%.2f,%.2f", getNumObservations(), getNumMismatches(), getEmpiricalQuality());
    }

    public String stringForCSV() {
        return String.format("%s,%.2f,%.2f", toString(), getEstimatedQReported(), getEmpiricalQuality() - getEstimatedQReported());
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // increment methods
    //
    //---------------------------------------------------------------------------------------------------------------

    public final long getNumObservations() {
        return numObservations;
    }

    public final void setNumObservations(final long numObservations) {
        if ( numObservations < 0 ) throw new IllegalArgumentException("numObservations < 0");
        this.numObservations = numObservations;
        empiricalQuality = UNINITIALIZED;
    }

    public final double getNumMismatches() {
        return numMismatches/MULTIPLIER;
    }

    public final void setNumMismatches(final double numMismatches) {
        if ( numMismatches < 0 ) throw new IllegalArgumentException("numMismatches < 0");
        this.numMismatches = (numMismatches*MULTIPLIER);
        empiricalQuality = UNINITIALIZED;
    }

    public final void incrementNumObservations(final long by) {
        numObservations += by;
        empiricalQuality = UNINITIALIZED;
    }

    public final void incrementNumMismatches(final double by) {
        numMismatches += (by*MULTIPLIER);
        empiricalQuality = UNINITIALIZED;
    }

    public final void increment(final long incObservations, final double incMismatches) {
        numObservations += incObservations;
        numMismatches += (incMismatches*MULTIPLIER);
        empiricalQuality = UNINITIALIZED;
    }

    public final void increment(final boolean isError) {
        increment(1, isError ? 1.0 : 0.0);
    }

    // -------------------------------------------------------------------------------------
    //
    // Private implementation helper functions
    //
    // -------------------------------------------------------------------------------------

    /**
     * calculate the expected number of errors given the estimated Q reported and the number of observations
     * in this datum.
     *
     * @return a positive (potentially fractional) estimate of the number of errors
     */
    private double calcExpectedErrors() {
        return getNumObservations() * QualityUtils.qualToErrorProb(estimatedQReported);
    }

    /**
     * Calculate and cache the empirical quality score from mismatches and observations (expensive operation)
     */
    private void calcEmpiricalQuality(final double conditionalPrior) {

        // smoothing is one error and one non-error observation
        final long mismatches = (long)(getNumMismatches() + 0.5) + SMOOTHING_CONSTANT;
        final long observations = getNumObservations() + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT;

        final double empiricalQual = RecalDatum.bayesianEstimateOfEmpiricalQuality(observations, mismatches, conditionalPrior);

        empiricalQuality = Math.min(empiricalQual, (double) MAX_RECALIBRATED_Q_SCORE);
    }

    //static final boolean DEBUG = false;
    private static final double RESOLUTION_BINS_PER_QUAL = 1.0;

    public static double bayesianEstimateOfEmpiricalQuality(final long nObservations, final long nErrors, final double QReported) {

        final int numBins = (QualityUtils.MAX_REASONABLE_Q_SCORE + 1) * (int)RESOLUTION_BINS_PER_QUAL;

        final double[] logPosteriors = new double[numBins];

        for ( int bin = 0; bin < numBins; bin++ ) {

            final double QEmpOfBin = bin / RESOLUTION_BINS_PER_QUAL;

            logPosteriors[bin] = logQempPrior(QEmpOfBin, QReported) + logQempLikelihood(QEmpOfBin, nObservations, nErrors);
        }

        final double[] normalizedPosteriors = MathUtils.normalizeFromLog(logPosteriors);
        final int MLEbin = MathUtils.maxElementIndex(normalizedPosteriors);

        final double Qemp = MLEbin / RESOLUTION_BINS_PER_QUAL;
        return Qemp;
    }

    /**
     * Quals above this value should be capped down to this value (because they are too high)
     * in the base quality score recalibrator
     */
    public static final byte MAX_GATK_USABLE_Q_SCORE = 40;
    private static final double[] logQempPriorCache = new double[MAX_GATK_USABLE_Q_SCORE + 1];
    static {
        // f(x) = a*exp(-((x - b)^2 / (2*c^2)))
        // Note that a is the height of the curve's peak, b is the position of the center of the peak, and c controls the width of the "bell".
        final double GF_a = 0.9;
        final double GF_b = 0.0;
        final double GF_c = 0.5;   // with these parameters, deltas can shift at most ~20 Q points

        final Gaussian gaussian = new Gaussian(GF_a, GF_b, GF_c);
        for ( int i = 0; i <= MAX_GATK_USABLE_Q_SCORE; i++ ) {
            double logPrior = Math.log(gaussian.value((double) i));
            if ( Double.isInfinite(logPrior) )
                logPrior = -Double.MAX_VALUE;
            logQempPriorCache[i] = logPrior;
        }
    }

    protected static double logQempPrior(final double Qempirical, final double Qreported) {
        final int difference = Math.min(Math.abs((int) (Qempirical - Qreported)), MAX_GATK_USABLE_Q_SCORE);
        return logQempPriorCache[difference];
    }

    private static final long MAX_NUMBER_OF_OBSERVATIONS = Integer.MAX_VALUE - 1;

    protected static double logQempLikelihood(final double Qempirical, long nObservations, long nErrors) {
        if ( nObservations == 0 )
            return 0.0;

        // the binomial code requires ints as input (because it does caching).  This should theoretically be fine because
        // there is plenty of precision in 2^31 observations, but we need to make sure that we don't have overflow
        // before casting down to an int.
        if ( nObservations > MAX_NUMBER_OF_OBSERVATIONS ) {
            // we need to decrease nErrors by the same fraction that we are decreasing nObservations
            final double fraction = (double)MAX_NUMBER_OF_OBSERVATIONS / (double)nObservations;
            nErrors = Math.round((double) nErrors * fraction);
            nObservations = MAX_NUMBER_OF_OBSERVATIONS;
        }

        // this is just a straight binomial PDF
        double logProb = MathUtils.logBinomialProbability((int) nObservations, (int) nErrors, QualityUtils.qualToLogErrorProb(Qempirical));
        if ( Double.isInfinite(logProb) || Double.isNaN(logProb) )
            logProb = -Double.MAX_VALUE;

        return logProb;
    }
}
