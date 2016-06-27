package org.broadinstitute.hellbender.tools.exome.coveragestats;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

/**
 * Per sample coverage summary statistics.
 * <p>
 *     This tuple holds a immutable reference to the sample name and the mean and sample variance across targets.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class SampleCoverageStats {

    public static final String SAMPLE_COLUMN_NAME = "SAMPLE";

    public static final String MEAN_COLUMN_NAME = "MEAN";

    public static final String VARIANCE_COLUMN_NAME = "VARIANCE";

    public static final TableColumnCollection COLUMNS = new TableColumnCollection(SAMPLE_COLUMN_NAME, MEAN_COLUMN_NAME, VARIANCE_COLUMN_NAME);

    /**
     * Name of the sample the statistics refer to.
     * <p>This member is never {@code null}</p>
     */
    public final String sample;

    /**
     * The mean coverage across targets.
     */
    public final double mean;

    /**
     * The sample variance of the coverage across targets.
     */
    public final double variance;

    /**
     * Create a instance given values for the sample name and its statistics.
     *
     * <p>
     *     If the variance pass is negative, the output object will have the corresponding positive value
     *     as the variance.
     * </p>
     *
     * @param sample the sample name.
     * @param mean the mean across targets.
     * @param variance the sample variance across targets.
     *
     * @throws IllegalArgumentException if {@code sample} is {@code null}.
     */
    public SampleCoverageStats(final String sample, final double mean, final double variance) {
        this.sample = Utils.nonNull(sample, "the sample cannot be null");
        this.mean = mean;
        this.variance = variance < 0 ? -variance : variance;
    }

    /**
     * Creates a sample coverage stats from the sum of all the element coverages and their squares.
     * <p>
     *     This method, to some extends, checks that the values for the sum and squared sum are
     *     consistent between it each other and with the size provided.
     * </p>
     * <p>
     *     Therefore, for some trivial inputs, this method might throw an IllegalArgumentException
     *     if indeed the sum, squaresSum and size provided do not correspond to a real sequence of
     *     values of that size that would generate such a
     *     sum and squaresSum.
     * </p>
     * <p>
     *     However this method won't go as far as to establish whether such a sequence really exist for
     *     every possible input.
     * </p>
     * @param sample sample name.
     * @param size number of elements.
     * @param sum the sum of coverages.
     * @param squaresSum the sum of squares.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code sample} is {@code null}.
     */
    public static SampleCoverageStats fromSums(final String sample, final long size, final double sum,
                                               final double squaresSum) {
        if (size < 0) {
            throw new IllegalArgumentException("the size cannot be less than 0");
        } else if (size == 0) {
            if (MathUtils.compareDoubles(sum, 0) != 0 && MathUtils.compareDoubles(squaresSum, 0) != 0) {
                throw new IllegalArgumentException("when the size is 0, then the sums must be zero");
            }
        } else if (size == 1) {
            if (MathUtils.compareDoubles(sum * sum, squaresSum) != 0) {
                throw new IllegalArgumentException("when the size is 1, then the sum squared must be the same as the squaresSum provided");
            }
        } else if (squaresSum < 0) {
            throw new IllegalArgumentException("the sum of squares cannot be negative");
        }
        final double mean = sum / size;
        final double squaresMean = squaresSum / size;
        final double populationVariance = squaresMean - (mean * mean);
        final double sampleVariance = size == 1 ? 0 : populationVariance * size / (size - 1);
        return new SampleCoverageStats(sample, mean, sampleVariance);
    }
}
