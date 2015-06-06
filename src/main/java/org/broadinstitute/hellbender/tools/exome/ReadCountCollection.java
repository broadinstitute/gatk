package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Represents a read-count collections.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ReadCountCollection {

    /**
     * Unmodifiable target list in the row order of their counts in {@link #counts}.
     */
    private final List<Target> targets;

    /**
     * Unmodifiable column name list in the column order of their counts in {@link #counts}.
     */
    private final List<String> columnNames;

    /**
     * Read count per target and column name.
     *
     * <p>
     *     Each row represents the counts for the ith target in {@link #targets}.
     * </p>
     *
     * <p>
     *     Each column represents the counts for the ith column name in {@link #columnNames}.
     * </p>
     */
    private final RealMatrix counts;

    /**
     * Creates a new read-counts collection.
     *
     * <p>
     *     All values in the input are copied, so these can be changed after invoking the constructor without
     *     affecting the new instance.
     * </p>
     *
     * <p>
     *     The new instance will have its own copy of the target, column-name list and counts. Therefore the input arguments
     *     can be modified after this call safely.
     * </p>
     *
     * @param targets the targets.
     * @param columnNames the count column names.
     * @param counts the read counts, with the same number for rows as {@code targets}, the same number of columns as
     *               {@code column names}.
     * @throws IllegalArgumentException if any of these is true:
     * <ul>
     *     <li>{@code targets} is {@code null},</li>
     *     <li>{@code columnNames} is {@code null},</li>
     *     <li>{@code counts} is {@code null},</li>
     *     <li>{@code targets} contains any {@code null}},</li>
     *     <li>{@code columnNames} contains any {@code null}},</li>
     *     <li>{@code counts} contains any {@code null}},</li>
     *     <li>{@code targets} length does not match {@code counts} length or</li>
     *     <li>{@code counts} elements length does not match {@code columnNames} length</li>
     * </ul>
     */
    public ReadCountCollection(final SetUniqueList<Target> targets, final SetUniqueList<String> columnNames, final RealMatrix counts) {
        Utils.nonNull(targets,"the input targets cannot be null");
        Utils.nonNull(columnNames,"the column names cannot be null");
        Utils.nonNull(counts,"the counts cannot be null");
        if (columnNames.contains(null)) {
            throw new IllegalArgumentException("column names contains nulls");
        } else if (targets.contains(null)) {
            throw new IllegalArgumentException("there is some null targets");
        } else if (counts.getRowDimension() != targets.size()) {
            throw new IllegalArgumentException("number of count rows does not match the number of targets");
        } else if (counts.getColumnDimension() != columnNames.size()) {
            throw new IllegalArgumentException("number of count columns does not match the number of column names");
        }
        this.targets = Collections.unmodifiableList(new ArrayList<>(targets));
        this.columnNames = Collections.unmodifiableList(new ArrayList<>(columnNames));
        this.counts = counts.copy();
    }

    /**
     * Returns the targets in the order they are found in this collection.
     * @return never {@code null}, and unmodifiable and immutable list of non-null targets.
     */
    public List<Target> targets() {
        return targets;
    }

    /**
     * Returns the list of count column names.
     * @return never {@code null}.
     */
    public List<String> columnNames() {
        return columnNames;
    }

    /**
     * Returns a live {@link RealMatrix} representation of the counts in this collection.
     *
     * <p>
     *     The ith row corresponds to the counts for the ith target in the list returned
     *     by {@link #targets()}.
     * </p>
     *
     * <p>
     *     The jth column corresponds to the counts for the jth target in list returned
     *     by {@link #columnNames()}.
     * </p>
     *
     * <p>
     *     This a live object and modifications of its contents modifies the count on this collection.
     * </p>
     *
     * @return the result matrix is a mutable live copy of the counts in this collection.
     */
    public RealMatrix counts() {
        return counts;
    }
}
