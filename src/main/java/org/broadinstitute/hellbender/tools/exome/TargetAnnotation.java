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
     * in this case it value will be {@link Double#NaN NaN}.</p>
     */
    GC_CONTENT(TargetTableColumn.GC_CONTENT),

    /**
     * Number of baits overlapping with the target. Only relevant in the context of target-enriched (e.g. WES) sequencing data.
     *
     * <p>This value is never {@link Double#NaN} and is non-negative. It can be zero or fractional (if a bait is shared by more
     * than a single target).</p>
     */
    BAIT_COUNT(TargetTableColumn.BAIT_COUNT);

    /**
     * Reference to the corresponding column in the targets table.
     */
    public final TargetTableColumn column;

    /**
     * Constructs a new annotation enum value.
     * @param column the corresponding column in the targets table.
     * @throws IllegalArgumentException if {@code column} is {@code null}.
     */
    TargetAnnotation(final TargetTableColumn column) {
        this.column = Utils.nonNull(column);
    }
}
