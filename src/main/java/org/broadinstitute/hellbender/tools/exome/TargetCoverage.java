package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;

/**
 * Exome analysis target with coverage
 *
 * <p>A name entity with genomic interval and coverage</p>
 *
 * @author David Benjamin
 */

public final class TargetCoverage extends Target {

    static final long serialVersionUID = 1111337337337L;

    private double coverage;

    /**
     * Construct a new TargetCoverage given all its properties.
     * @param interval the interval.
     * @param name the name of the interval.
     * @param coverage the coverage (read count or normalized coverage) of this interval
     *
     * @throws IllegalArgumentException if either {@code interval} or {@code start}
     *    and {@code end} do not represent a valid interval
     */
    public TargetCoverage(final String name, final SimpleInterval interval, final double coverage) {
        super(name, interval);
        this.coverage = coverage;
    }

    /**
     * Returns the coverage.
     */
    public double getCoverage() {
        return coverage;
    }

    /**
     * Sets the coverage.
     */
    public void setCoverage(final double coverage) {
        this.coverage = coverage;
    }
}
