package org.broadinstitute.hellbender.tools.exome.coveragestats;

import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableColumn;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Target coverage statistics across samples.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetCoverageStats {

    /**
     * Name of the column for the means.
     */
    public static final String MEAN_COLUMN_NAME = TargetTableColumn.MEAN_COVERAGE.toString();

    /**
     * Name of the column for the variance.
     */
    public static final String VARIANCE_COLUMN_NAME = TargetTableColumn.COVERAGE_VARIANCE.toString();

    /**
     * Name of the column for the interquartile range.
     */
    public static final String INTERQUARTILE_RANGE_COLUMN_NAME = TargetTableColumn.COVERAGE_INTERQUARTILE_RANGE.toString();

    /**
     * Target the coverage statistics makes reference to.
     * <p>
     * This is never {@code null}.
     * </p>
     */
    public final Target target;

    /**
     * The mean coverage for {@link #target} across samples.
     */
    public final double mean;

    /**
     * The coverage sample variance for {@link #target} across samples.
     */
    public final double variance;

    /**
     * The coverage interquartile range for {@link #target} across samples.
     */
    public final double interquartileRange;


    /**
     * Creates a new instance given its member values.
     * @param target the target this statistics make reference to.
     * @param mean the mean across samples.
     * @param variance the sample variance across samples.
     * @param interquartileRange the interquartile range across samples
     */
    TargetCoverageStats(final Target target, final double mean, final double variance, final double interquartileRange) {
        this.target = Utils.nonNull(target, "the target cannot be null");
        this.mean = mean;
        this.variance = variance;
        this.interquartileRange = interquartileRange;
    }

    /**
     * Creates a new instance given the target and the coverage for that target across count columns.
     * @param target the target.
     * @param coverage the coverage for the target across count columns.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code target} is {@code null} or {@code coverage} is {@code null}.
     */
    public static TargetCoverageStats fromCoverage(final Target target, final double ... coverage)  {
        Utils.nonNull(coverage, "the coverage cannot be null");
        final int size = coverage.length;
        double sum = 0;
        double squaresSum = 0;
        for (final double value : coverage) {
            sum += value;
            squaresSum += value * value;
        }
        final double mean = sum / size;
        final double squaresMean = squaresSum / size;
        final double populationVariance = squaresMean - (mean * mean);
        final double sampleVariance = size == 1 ? 0 : populationVariance * size / (size - 1);
        final double interquartileRange = GATKProtectedMathUtils.interquartileRange(coverage);
        return new TargetCoverageStats(target, mean, sampleVariance, interquartileRange);
    }
}
