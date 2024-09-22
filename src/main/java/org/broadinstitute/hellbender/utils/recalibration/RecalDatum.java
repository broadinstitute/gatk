package org.broadinstitute.hellbender.utils.recalibration;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.io.Serializable;

/**
 * An individual piece of recalibration data. Each bin counts up the number of observations and the number
 * of reference mismatches seen for that combination of covariates.
 */
public final class RecalDatum implements Serializable {
    public static final byte MAX_RECALIBRATED_Q_SCORE = SAMUtils.MAX_PHRED_SCORE;
    private static final int UNINITIALIZED_EMPIRICAL_QUALITY = -1;
    private static final long serialVersionUID = 1L;
    private static final double MULTIPLIER = 100000.0;  //See discussion in numMismatches about what the multiplier is.

    /** estimated reported quality score based on combined data's individual q-reporteds and number of observations
     */ // tsato: rename to reportedQuality (forget the estimating when combining; that's just implementation detail for the edge case)
    private double estimatedQReported; // tsato: estimated by Illumina? Should be renamed to OriginalQuality
    // tsato: This column only occurs in ReadGroupTable ..... so consider inheriting RecalDatum.


    /**
     * the empirical quality for datums that have been collapsed together (by read group and reported quality, for example)
     */ // tsato: defaults to the "QualityScore"
    private int empiricalQuality; // tsato: seems like it should have type integers

    /**
     * number of bases seen in total
     */
    private long numObservations;

    /**
     * Number of bases seen that didn't match the reference
     * (actually sum of the error weights - so not necessarily a whole number)
     * Stored with an internal multiplier to keep it closer to the floating-point sweet spot and avoid numerical error
     * (see https://github.com/broadinstitute/gatk/wiki/Numerical-errors ).
     * However, the value of the multiplier influences the results.
     * For example, you get different results for 1000.0 and 10000.0
     * See MathUtilsUnitTest.testAddDoubles for a demonstration.
     * The value of the MULTIPLIER that we found to give consistent results insensitive to sorting is 10000.0;
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
     * @param reportedQuality     Qreported ???
     */
    public RecalDatum(final long _numObservations, final double _numMismatches, final byte reportedQuality) {
        if ( _numObservations < 0 ) throw new IllegalArgumentException("numObservations < 0");
        if ( _numMismatches < 0.0 ) throw new IllegalArgumentException("numMismatches < 0");
        if ( reportedQuality < 0 ) throw new IllegalArgumentException("reportedQuality < 0");

        numObservations = _numObservations;
        numMismatches = (_numMismatches*MULTIPLIER);
        estimatedQReported = reportedQuality; // tsato: estimated is very confusing here
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
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
    public void combine(final RecalDatum other) { // tsato: when is data combined?
        final double sumErrors = this.calcExpectedErrors() + other.calcExpectedErrors(); // tsato: ok I think I get it now but why would we combine different quals?
        increment(other.getNumObservations(), other.getNumMismatches()); // tsato: very questionable here for combining, but it seems that the only time we combine is GatherBQSRTables (confirm)
        estimatedQReported = -10 * Math.log10(sumErrors / getNumObservations()); // tsato: aren't we mixing up estimatedQReported with EmpiricalQuality here?
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    public void setEstimatedQReported(final double estimatedQReported) {
        if ( estimatedQReported < 0 ) throw new IllegalArgumentException("estimatedQReported < 0");
        if ( Double.isInfinite(estimatedQReported) ) throw new IllegalArgumentException("estimatedQReported is infinite");
        if ( Double.isNaN(estimatedQReported) ) throw new IllegalArgumentException("estimatedQReported is NaN");

        this.estimatedQReported = estimatedQReported;
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
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

    public void setEmpiricalQuality(final int empiricalQuality) {
        if ( empiricalQuality < 0 ) throw new IllegalArgumentException("empiricalQuality < 0");
        if ( Double.isInfinite(empiricalQuality) ) throw new IllegalArgumentException("empiricalQuality is infinite");
        if ( Double.isNaN(empiricalQuality) ) throw new IllegalArgumentException("empiricalQuality is NaN");

        this.empiricalQuality = empiricalQuality;
    }

    public double getEmpiricalQuality() {
        return getEmpiricalQuality(getEstimatedQReported());
    }

    /**
     * Each datum (either for read group covariate, or special covariate) has a way to compute
     * the empirical quality.
     *
     * @param conditionalPrior
     * @return
     */
    public double getEmpiricalQuality(final double conditionalPrior) {
        if (empiricalQuality == UNINITIALIZED_EMPIRICAL_QUALITY) {
            calcEmpiricalQuality(conditionalPrior);
        }
        return empiricalQuality;
    }

    public byte getEmpiricalQualityAsByte() {
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
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    public final double getNumMismatches() {
        return numMismatches/MULTIPLIER;
    }

    public final void setNumMismatches(final double numMismatches) {
        if ( numMismatches < 0 ) throw new IllegalArgumentException("numMismatches < 0");
        this.numMismatches = (numMismatches*MULTIPLIER);
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    public final void incrementNumObservations(final long by) {
        numObservations += by;
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    public final void incrementNumMismatches(final double by) {
        numMismatches += (by*MULTIPLIER);
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY; // tsato: ditto below
    }

    public final void increment(final long incObservations, final double incMismatches) {
        numObservations += incObservations;
        numMismatches += (incMismatches*MULTIPLIER); // tsato: Multiplier used to avoid underflow, or something like that.
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY; // tsato: why not just delete this?
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
        final long mismatches = (long)(getNumMismatches() + 0.5) + SMOOTHING_CONSTANT; // tsato: why add 0.5? can I get rid of it?
        final long observations = getNumObservations() + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT; // tsato: why twice...a typo? No it's ok.
        // tsato: crux ---- done.
        final int empiricalQual = RecalDatum.bayesianEstimateOfEmpiricalQuality(observations, mismatches, conditionalPrior);

        empiricalQuality = Math.min(empiricalQual, MAX_RECALIBRATED_Q_SCORE);
    }

    /**
     * tsato: to fill in
     * @param nObservations
     * @param nErrors
     * @param QReported (is "conditionalPrior" ....
     * @return
     */
    public static int bayesianEstimateOfEmpiricalQuality(final long nObservations, final long nErrors, final double QReported) {
        // tsato: maxReasonable q score = 60. resolution = 1...
        final int numBins = (QualityUtils.MAX_REASONABLE_Q_SCORE + 1);
        // tsato: so num bins is 61
        final double[] logPosteriors = new IndexRange(0, numBins).mapToDouble(bin -> {
            final double QEmpOfBin = bin; // tsato: but this is bin is ... empirical? // tsato: this just casts int to double; necessary?
            return logQempPrior(QEmpOfBin, QReported) + logQempLikelihood(QEmpOfBin, nObservations, nErrors);
        });
        final int MAPbin = MathUtils.maxElementIndex(logPosteriors);
        return MAPbin;
    }

    /**
     * Quals above this value should be capped down to this value (because they are too high)
     * in the base quality score recalibrator
     */
    public static final byte MAX_GATK_USABLE_Q_SCORE = 40;
    private static final double[] logQempPriorCache = new double[MAX_GATK_USABLE_Q_SCORE + 1];
    static { // tsato: huh...difference in log space is modeled as a Gaussian distribution...
        // normal distribution describing P(Q empirical - Q reported).  Its mean is zero because a priori we expect // tsato: this is David's comment
        // no systematic bias in the reported quality score
        final double mean = 0.0;
        final double sigma = 0.5;   // with these parameters, deltas can shift at most ~20 Q points
        final NormalDistribution gaussian = new NormalDistribution(null, mean, sigma);

        for ( int i = 0; i <= MAX_GATK_USABLE_Q_SCORE; i++ ) {
            logQempPriorCache[i] = gaussian.logDensity(i);
        }
    }

    /**
     *
     * @param Qempirical ??
     * @param Qreported ???
     * @return
     */
    @VisibleForTesting
    protected static double logQempPrior(final double Qempirical, final double Qreported) { // tsato: David's change is to go from log10 to log
        final int difference = Math.min(Math.abs((int) (Qempirical - Qreported)), MAX_GATK_USABLE_Q_SCORE); // tsato: why cast to integer here...
        return logQempPriorCache[difference]; // tsato: revisit
    }

    private static final long MAX_NUMBER_OF_OBSERVATIONS = Integer.MAX_VALUE - 1;

    @VisibleForTesting
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
        final double logLikelihood = MathUtils.logBinomialProbability((int) nObservations, (int) nErrors, QualityUtils.qualToErrorProb(Qempirical));
        return ( Double.isInfinite(logLikelihood) || Double.isNaN(logLikelihood) ) ? -Double.MAX_VALUE : logLikelihood;
    }
}
