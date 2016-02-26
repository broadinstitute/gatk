package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Enumeration of target annotations.
 * <p>
 *     A target annotation is a property of the target but not including the standard target name and coordinates.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public enum TargetAnnotation {

    /**
     * Fraction of bases that are G or C in the reference overlapping the target.
     *
     * <p>This value is always a number between 0 and 1 except in very rare cases with small targets with all, undefined bases (N or X),
     * in this case it value will be {@link Double#NaN NaN}</p>.
     */
    GC_CONTENT(TargetTableColumns.GC_CONTENT),

    /**
     * Fraction of bases in the target that have been marked as repeated somewhere else in the genome.
     *
     * <p>This value is always a number between 0 and 1 without exceptions.</p>
     */
    REPEAT_FRACTION(TargetTableColumns.REPEAT_FRACTION),

    /**
     * Average coverage across samples.
     */
    MEAN_COVERAGE(TargetTableColumns.MEAN_COVERAGE),

    /**
     * The coverage variance across samples.
     */
    COVERAGE_VARIANCE(TargetTableColumns.COVERAGE_VARIANCE),

    /**
     * The coverage interquartile range across samples.
     */
    COVERAGE_INTERQUARTILE_RANGE(TargetTableColumns.COVERAGE_INTERQUARTILE_RANGE);

    /**
     * Reference to the corresponding column in the targets table.
     */
    public final TargetTableColumns column;

    /**
     * Constructs a new annotation enum value.
     * @param column the corresponding column in the targets table.
     * @throws IllegalArgumentException if {@code column} is {@code null}.
     */
    TargetAnnotation(final TargetTableColumns column) {
        this.column = Utils.nonNull(column);
    }
}
