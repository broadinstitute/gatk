package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Exome analysis target.
 *
 * <p>A name entity with genomic interval optionally attached to it</p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class Target {

    /**
     * The target name. Always specified; never {@code null}.
     */
    private final String name;

    /**
     * Can be left unspecified ({@code null}).
     */
    private final SimpleInterval interval;

    /**
     * Construct a new interval given all its properties.
     * @param interval the interval.
     * @param name the name of the interval.
     *
     * @throws IllegalArgumentException if either {@code interval} or {@code start}
     *    and {@code end} do not represent a valid interval as described in
     */
    public Target(final String name, final SimpleInterval interval) {
        Utils.nonNull(name, "the name cannot be null");
        this.name = name;
        this.interval = interval;
    }

    /**
     * Returns the interval.
     * @return may be {@code null}.
     */
    public SimpleInterval getInterval() {
        return interval;
    }

    /**
     * Returns the name.
     * @return never{@code null}.
     */
    public String getName() {
        return name;
    }

    @Override
    public boolean equals(final Object other) {
        return other instanceof Target ? equals((Target)other) : false;
    }

    /**
     * Compare this with another target.
     *
     * <p>
     *     Only the name is compared; the interval is not considered par of the identification of the target.
     * </p>
     * @param other the other target to compare to.
     * @return {@code true} iff they have the same interval and name.
     */
    public boolean equals(final Target other) {
        if (other == null) {
            return false;
        } else {
            return name.equals(other.name);
        }
    }

    @Override
    public int hashCode() {
        return name.hashCode();
    }

}
