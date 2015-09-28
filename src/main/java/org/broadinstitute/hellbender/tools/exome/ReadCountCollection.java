package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

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
     * Creates a new collection without verifying field values and without copying inputs.
     *
     * <p>
     * The field values are supposed to be compatible with a consistent state.
     * </p>
     * @param targets target list, not a {@code null}, does not contain any {@code null}, does not contain repeats.
     * @param columnNames column name list, not a {@code null}, does not contain any {@code null}, does not contain repeats.
     * @param counts count matrix, not a {@code null}, has as many rows as {@code targets} elements and as many columns as {@code columnNames} elements.
     */
    private ReadCountCollection(final List<Target> targets, final List<String> columnNames, final RealMatrix counts) {
        this.targets = targets;
        this.columnNames = columnNames;
        this.counts = counts;
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

    /**
     * Subsets the targets in the read-count collection.
     * <p>
     *     Creates  brand-new read-count collection. Changes in the new read-count collection
     *     counts won't affect the this read-count collection and vice-versa.
     * </p>
     *
     * @param targetsToKeep the new target subset.
     * @return never {@code null}. The order of targets in the result is guaranteed to
     *  follow the traversal order of {@code targetsToKeep}. The order of count columns is
     *  guaranteed to follow the original order of count columns.
     * @throws IllegalArgumentException if {@code targetsToKeep}:
     * <ul>
     *     <li>is {@code null},</li>
     *     <li>contains {@code null}s</li>
     *     <li>or contains targets that are not part of the read-count collection</li>
     * </ul>
     */
    public ReadCountCollection subsetTargets(final Set<Target> targetsToKeep) {
        Utils.nonNull(targetsToKeep, "the input target set cannot be null");
        if (targetsToKeep.isEmpty()) {
            throw new IllegalArgumentException("the input target subset size must be greater than 0");
        }
        if (targetsToKeep.size() > targets.size()) {
            throw unknownTargetsToKeep(targetsToKeep);
        }
        if (targetsToKeep.size() == targets.size())  {
            if (targets.stream().anyMatch(name -> !targetsToKeep.contains(name))) {
                throw unknownTargetsToKeep(targetsToKeep);
            } else {
                return new ReadCountCollection(targets, columnNames, counts.copy());
            }
        }
        final int[] targetsToKeepIndices = new int[targetsToKeep.size()];
        int nextIndex = 0;
        final List<Target> resultTargets = new ArrayList<>(targetsToKeep.size());
        for (int i = 0; i < targets.size(); i++) {
            final Target target = targets.get(i);
            if (!targetsToKeep.contains(target)) {
                continue;
            }

            if (nextIndex >= targetsToKeepIndices.length) {
                throw unknownTargetsToKeep(targetsToKeep);
            } else {
                targetsToKeepIndices[nextIndex++] = i;
                resultTargets.add(target);
            }
        }
        // check that all targets to be kept were found in this collection:
        if (nextIndex < targetsToKeep.size()) {
            throw unknownTargetsToKeep(targetsToKeep);
        }
        // compose the new counts:
        final double[][] resultCounts = new double[targetsToKeepIndices.length][columnNames.size()];
        for (int i = 0; i < resultCounts.length; i++) {
            resultCounts[i] = counts.getRow(targetsToKeepIndices[i]);
        }
        return new ReadCountCollection(Collections.unmodifiableList(resultTargets), columnNames, new Array2DRowRealMatrix(resultCounts));
    }

    /**
     * Subsets the count columns in the read-count collection.
     *
     * <p>
     *     Creates a brand-new read-count collection. Changes in the new instance won't affect this one and vice-versa.
     * </p>
     *
     * @param columnsToKeep column names to keep in the result read-count collection.
     * @return never {@code null}.
     */
    public ReadCountCollection subsetColumns(final Set<String> columnsToKeep) {
        Utils.nonNull(columnsToKeep, "the set of input columns to keep cannot be null.");
        if (columnsToKeep.isEmpty()) {
            throw new IllegalArgumentException("the number of columns to keep must be greater than 0");
        }
        if (columnsToKeep.size() > columnNames.size()) {
            throw unknownColumnToKeepNames(columnsToKeep);
        }
        if (columnsToKeep.size() == columnNames.size())  {
            if (columnNames.stream().anyMatch(name -> !columnsToKeep.contains(name))) {
                throw unknownColumnToKeepNames(columnsToKeep);
            } else {
                return new ReadCountCollection(targets, columnNames, counts.copy());
            }
        }
        final int[] columnsToKeepIndices = new int[columnsToKeep.size()];
        int nextIndex = 0;
        final List<String> resultColumnNames = new ArrayList<>(columnsToKeep.size());
        for (int i = 0; i < columnNames.size(); i++) {
            final String column = columnNames.get(i);
            if (!columnsToKeep.contains(column)) {
                continue;
            }
            if (nextIndex >= columnsToKeepIndices.length) {
                throw unknownColumnToKeepNames(columnsToKeep);
            } else {
                columnsToKeepIndices[nextIndex++] = i;
                resultColumnNames.add(column);
            }
        }
        // check that all the columns to kept were found in the this collection:
        if (nextIndex < columnsToKeep.size()) {
            throw unknownColumnToKeepNames(columnsToKeep);
        }
        // compose the new counts:
        final double[][] resultCounts = new double[counts.getRowDimension()][columnsToKeepIndices.length];
        for (int i = 0; i < resultCounts.length; i++) {
            final double[] rowCounts = counts.getRow(i);
            for (int j = 0; j < columnsToKeepIndices.length; j++) {
                resultCounts[i][j] = rowCounts[columnsToKeepIndices[j]];
            }
        }
        return new ReadCountCollection(targets, Collections.unmodifiableList(resultColumnNames), new Array2DRowRealMatrix(resultCounts, false));
    }

    /**
     * Constructs the appropriate exception to report the presence of column names in the columns-to-keep set that are not present in
     * this read-count collection.
     * @param columnsToKeep the columns to keep set.
     * @return never {@code null}.
     */
    private IllegalArgumentException unknownColumnToKeepNames(final Set<String> columnsToKeep) {
        return new IllegalArgumentException("some column names in the column keep set that are not part of this read count collection: e.g. "
                + columnsToKeep.stream().filter(name -> !columnNames.contains(name)).limit(5).collect(Collectors.joining(", ")));
    }

    /**
     * Constructs the appropriate exception to report the presence of targets in the targets-to-keep set that are not present in
     * this read-count collection.
     * @param targetsToKeep the columns to keep set.
     * @return never {@code null}.
     */
    private IllegalArgumentException unknownTargetsToKeep(final Set<Target> targetsToKeep) {
        return new IllegalArgumentException("some column names in the column keep set that are not part of this read count collection: e.g. "
                + targetsToKeep.stream().filter(name -> !targets.contains(name)).map(Target::getName).limit(5).collect(Collectors.joining(", ")));
    }
}
