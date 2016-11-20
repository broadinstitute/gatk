package org.broadinstitute.hellbender.tools.exome;

import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.DefaultRealMatrixPreservingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents a read-count collections.
 *
 * Developer note: any public constructor of this class must verify that targets and column names do not contain duplicates.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ReadCountCollection implements Serializable {

    static final long serialVersionUID = 337337337L;

    /**
     * Unmodifiable target list in the row order of their counts in {@link #counts}.
     */
    private final List<Target> targets;

    /**
     * A map from targets to their indices in the list
     */
    private final Map<Target, Integer> targetIndexMap;

    /**
     * Unmodifiable column name list in the column order of their counts in {@link #counts}.
     */
    private final List<String> columnNames;

    /**
     * Read counts matrix with one row per target in {@link #targets} and columns corresponding to {@link #columnNames}.
     */
    private final RealMatrix counts;

    /**
     * Creates a new read-counts collection.
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
     *     <li>{@code targets} contains any duplicates,</li>
     *     <li>{@code columnNames} contains any duplicates,</li>
     *     <li>{@code targets} length does not match {@code counts} length or</li>
     *     <li>{@code counts} elements length does not match {@code columnNames} length</li>
     * </ul>
     */
    public ReadCountCollection(final List<Target> targets, final List<String> columnNames, final RealMatrix counts) {
        this(targets, columnNames, counts, true);
    }

    /**
     * Creates a new collection with or without verifying field values and copying inputs.
     *
     * <p>
     * The field values are supposed to be compatible with a consistent state.
     * </p>
     * @param targets target list, not a {@code null}, does not contain any {@code null}, does not contain repeats.
     * @param columnNames column name list, not a {@code null}, does not contain any {@code null}, does not contain repeats.
     * @param counts count matrix, not a {@code null}, has as many rows as {@code targets} elements and as many columns as {@code columnNames} elements.
     * @param verifyInput whether to check input for nulls and duplicates and make defensive copies
     */
    private ReadCountCollection(final List<Target> targets, final List<String> columnNames, final RealMatrix counts, final boolean verifyInput) {
        if (verifyInput) {
            Utils.nonNull(targets,"the input targets cannot be null");
            Utils.nonNull(columnNames,"the column names cannot be null");
            Utils.nonNull(counts,"the counts cannot be null");
            Utils.containsNoNull(columnNames, "column names contain nulls");
            Utils.containsNoNull(targets, "there are some null targets");
            Utils.validateArg(counts.getRowDimension() == targets.size(), "number of count rows does not match the number of targets");
            Utils.validateArg(counts.getColumnDimension() == columnNames.size(), "number of count columns does not match the number of column names");
            Utils.validateArg(new HashSet<>(targets).size() == targets.size(), "targets contain duplicates");
            Utils.validateArg(new HashSet<>(columnNames).size() == columnNames.size(), "column names contain duplicates");
            this.targets = Collections.unmodifiableList(new ArrayList<>(targets));
            this.columnNames = Collections.unmodifiableList(new ArrayList<>(columnNames));
            this.counts = counts.copy();
        } else {
            this.targets = targets;
            this.columnNames = columnNames;
            this.counts = counts;
        }
        this.targetIndexMap = createTargetIndexMap(this.targets);
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
     * Convert to a List of ReadCountRecords.  Note that this does not contain the information in the header i.e. column names.
     * @return never {@code null}
     */
    public List<ReadCountRecord> records() {
        return IntStream.range(0, targets.size())
                .mapToObj(t -> new ReadCountRecord(targets.get(t), counts.getRow(t)))
                .collect(Collectors.toList());
    }

    /**
     * Returns a column from the collection
     * @param columnIndex index of the column to return
     * @return a double array
     */
    public double[] getColumn(final int columnIndex) {
        ParamUtils.inRange(columnIndex, 0, counts.getColumnDimension() - 1, "Column index is out of range");
        return counts.getColumn(columnIndex);
    }

    /**
     * TODO unit test
     *
     * Returns a column from the collection on a given list of targets with the same ordering as {@code targetsToKeep}
     *
     * @param columnIndex index of the column to return
     * @param targetsToKeep list of targets
     * @return a double array
     */
    public double[] getColumnOnSpecifiedTargets(final int columnIndex, @Nonnull final List<Target> targetsToKeep,
                                                final boolean verifyInput) {
        if (verifyInput) {
            ParamUtils.inRange(columnIndex, 0, counts.getColumnDimension() - 1, "Column index out of range");
            Utils.nonEmpty(targetsToKeep, "The input target list can not be empty");
            if (!new HashSet<>(targets).containsAll(targetsToKeep)) {
                throw unknownTargetsToKeep(new HashSet<>(targetsToKeep));
            }
        }

        final double[] fullColumn = counts.getColumn(columnIndex);
        return targetsToKeep.stream().mapToDouble(t -> fullColumn[targetIndexMap.get(t)]).toArray();
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
    public RealMatrix counts() { return counts; }

    /**
     * Returns a map from targets to their indices
     *
     * @param targets target list
     * @return map
     */
    private Map<Target, Integer> createTargetIndexMap(@Nonnull final List<Target> targets) {
        final Map<Target, Integer> targetIndexMap = new HashMap<>();
        IntStream.range(0, targets.size()).forEach(i -> targetIndexMap.put(targets.get(i), i));
        return targetIndexMap;
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
        Utils.nonEmpty(targetsToKeep, "the input target subset size must be greater than 0");
        if (!new HashSet<>(targets).containsAll(targetsToKeep)) {
            throw unknownTargetsToKeep(targetsToKeep);
        }

        if (targetsToKeep.size() == targets.size())  {
            return new ReadCountCollection(targets, columnNames, counts.copy(), false);
        }
        final int[] targetsToKeepIndices = IntStream.range(0, targets.size())
                .filter(i -> targetsToKeep.contains(targets.get(i))).toArray();
        final List<Target> resultTargets = Arrays.stream(targetsToKeepIndices).mapToObj(targets::get).collect(Collectors.toList());

        // compose the new counts:
        final double[][] resultCounts = new double[targetsToKeepIndices.length][columnNames.size()];
        for (int i = 0; i < resultCounts.length; i++) {
            resultCounts[i] = counts.getRow(targetsToKeepIndices[i]);
        }
        return new ReadCountCollection(Collections.unmodifiableList(resultTargets), columnNames, new Array2DRowRealMatrix(resultCounts), false);
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
        Utils.nonEmpty(columnsToKeep, "the number of columns to keep must be greater than 0");
        if (!new HashSet<>(columnNames).containsAll(columnsToKeep)) {
            throw unknownColumnToKeepNames(columnsToKeep);
        }

        if (columnsToKeep.size() == columnNames.size())  {
            return new ReadCountCollection(targets, columnNames, counts.copy(), false);
        }

        final int[] columnsToKeepIndices = IntStream.range(0, columnNames.size())
                .filter(i -> columnsToKeep.contains(columnNames.get(i))).toArray();
        final List<String> resultColumnNames = Arrays.stream(columnsToKeepIndices).mapToObj(columnNames::get).collect(Collectors.toList());

        final RealMatrix resultCountsM = new Array2DRowRealMatrix(counts.getRowDimension(), columnsToKeepIndices.length);
        for (int i = 0; i < columnsToKeepIndices.length; i++) {
            resultCountsM.setColumn(i, counts.getColumn(columnsToKeepIndices[i]));
        }

        return new ReadCountCollection(targets, Collections.unmodifiableList(resultColumnNames), resultCountsM, false);
    }

    /**
     * Rearrange the targets so that they are in a particular order.
     * @return a new collection.
     * @throws IllegalArgumentException if any of the following is true:
     * <ul>
     *     <li>{@code targetsInOrder} is {@code null},</li>
     *     <li>is empty,</li>
     *     <li>it contains {@code null},</li>
     *     <li>contains any target not present in this collection.</li>
     * </ul>
     */
    public ReadCountCollection arrangeTargets(final List<Target> targetsInOrder) {
        Utils.nonNull(targetsInOrder);
        Utils.nonEmpty(targetsInOrder, "the input targets list cannot be empty");
        final RealMatrix counts = new Array2DRowRealMatrix(targetsInOrder.size(), columnNames.size());
        final Object2IntMap<Target> targetToIndex = new Object2IntOpenHashMap<>(targets.size());
        for (int i = 0; i < targets.size(); i++) {
            targetToIndex.put(targets.get(i), i);
        }
        for (int i = 0; i < targetsInOrder.size(); i++) {
            final Target target = targetsInOrder.get(i);
            Utils.validateArg(targetToIndex.containsKey(target), () -> String.format("target '%s' is not present in the collection", target.getName()));
            counts.setRow(i, this.counts.getRow(targetToIndex.getInt(target)));
        }
        return new ReadCountCollection(new ArrayList<>(targetsInOrder), columnNames, counts, false);
    }

    /**
     * Divide coverage at each target and each column by the average of that column.
     * @param weightByTargetSize whether to use a weighted average with weights given by target sizes
     * @return a new collection.
     * @throws IllegalArgumentException if {@code weightByTargetSize} is {@code true} but any target
     * is missing an interval.
     */
    public ReadCountCollection normalizeByColumnAverages(final boolean weightByTargetSize) {
        final RealMatrix normalizedCounts = counts().copy();
        final double[] columnWeightedMeans = calculateColumnWeightedMeans(weightByTargetSize);
        normalizedCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(final int target, final int column, final double coverage) {
                    return coverage / columnWeightedMeans[column];
                }
        });

        return new ReadCountCollection(targets, columnNames, normalizedCounts, false);
    }

    private double[] calculateColumnWeightedMeans(final boolean weightByTargetSize) {
        if (!weightByTargetSize) {
            return GATKProtectedMathUtils.columnMeans(counts);
        } else {
            final long[] weights;
            try {
                weights = targets.stream().mapToLong(Target::length).toArray();
            } catch (final IllegalStateException e) {
                throw new IllegalArgumentException("Weighting by target size requested but at least one target lacks an interval");
            }
            final long totalWeight = Arrays.stream(weights).sum();

            final double[] result = new double[counts.getColumnDimension()]; // elements initialized to 0.0

            counts.walkInOptimizedOrder(new DefaultRealMatrixPreservingVisitor() {
                @Override
                public void visit(final int target, final int column, final double coverage) {
                    result[column] += weights[target] * coverage / totalWeight;
                }
            });
            return result;
        }
    }

    /**
     * Express coverage in terms of Z scores with respect to the coverage distribution of the corresponding target.
     * @return a new collection.
     */
    public ReadCountCollection zScoreCounts() {
        final RealMatrix zScoreCounts = counts().copy();    //we will edit the matrix in-place

        final double[] targetMeans = GATKProtectedMathUtils.rowMeans(zScoreCounts);
        final double[] targetVariances = GATKProtectedMathUtils.rowVariances(zScoreCounts);
        final double[] targetStandardDeviations = Arrays.stream(targetVariances).map(Math::sqrt).toArray();

        zScoreCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int target, final int column, final double coverage) {
                return (coverage - targetMeans[target]) / targetStandardDeviations[target];
            }
        });

        return new ReadCountCollection(targets, columnNames, zScoreCounts, false);
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
