package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;
import java.util.stream.IntStream;

/**
 * Target collection.
 * <p>
 *      This interface includes operations to query targets based on genomic coordinates and their names (if provided).
 * </p>
 *
 * <p>
 *      Additional per target meta-data can be added by implementing class through the customizable type parameter
 *      &lt;E&gt;.
 * </p>
 *
 * @param <E> target meta-data type.
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */

public interface TargetCollection<E> {

    /**
     * Returns the number of targets in this collection.
     *
     * @return 0 or greater.
     */
    int targetCount();

    /**
     * Returns the number of targets that overlap a genomic location.
     * <p>
     * The number of targets that overlap the input {@code location} by at least one base are returned</p>
     *
     * @param overlapRegion the query location.
     * @return the number of targets that overlap the input {@code overlapRegion}
     * @throws IllegalArgumentException if {@code overlapRegion} is {@code null}.
     */
    default int targetCount(final Locatable overlapRegion) {
        return indexRange(overlapRegion).size();
    }

    /**
     * Returns a target given its index.
     * <p>
     * Indices are 0-based, thus the first target has index 0, the second 1 and so forth.
     * </p>
     *
     * @param index the query index.
     * @throws IndexOutOfBoundsException if {@code index} is not within
     *                                   valid bounds <code>[0 .. ({@link #targetCount()})</code>.
     */
    E target(final int index);

    /**
     * Returns the name of a target.
     *
     * <p>
     *     This method is guaranteed to return the name of the target in this collection as long as it is indeed
     *     included in it.
     * </p>
     *
     * <p>
     *     Otherwise it should return a reasonable name for it or {@code null}.
     * </p>
     *
     * @param target a target.
     *
     * @throws IllegalArgumentException if {@code target} is {@code null}.
     */
    String name(final E target);

    /**
     * Returns the target given its name.
     *
     * @param name the query target name.
     * @return {@code null} if there is not such a target.
     */
    default E target(final String name) {
        final int index = index(name);
        if (index < 0) {
            return null;
        } else {
            return target(index);
        }
    }

    /**
     * Returns the target index given its name.
     *
     * @param name the query target name.
     * @throws IllegalArgumentException if {@code name} is {@code null}.
     * @return {@code -1} if there is not such a target.
     */
    int index(final String name);

    /**
     * Returns the target that overlap a particular location.
     * <p>
     *     The returned target coordinates does not need to match exactly the
     * the query {@code location} just overlap it by at least one base</p>.
     * <p>
     *     This method will return {@code null} if there is no such a target.
     * </p>
     * <p>
     *     If there is more than on target that overlap with the given region it will
     *     throw a {@link AmbiguousTargetException} instead.
     * </p>
     *
     * @param overlapRegion the query location.
     * @return {@code null} if there are no overlapping targets at {@code overlapRegion}.
     * @throws IllegalArgumentException if {@code overlapRegion} is {@code null}.
     * @throws AmbiguousTargetException if the query location overlaps more than one target.
     */
    default E target(final Locatable overlapRegion) {
        final int index = index(overlapRegion);
        return index < 0 ? null : target(index);
    }

    /**
     * Returns the index of the target that overlap a particular location.
     * <p>
     *    The returned index's target coordinates does not need to match
     *    exactly the query {@code location} just overlap it by at
     *    least one base.
     * </p>
     * <p>
     *    This method will return a negative value if there is no such a target.
     *    This value can be transformed to recover the position where the query location would be inserted
     *    if it was present as <code>-(index(loc) + 1)</code>.
     * </p>
     *
     * <p>
     *    This method will result in a {@link AmbiguousTargetException} if there is more than one
     *    target that overlaps that query {@code location}.
     * </p>
     *
     * @param location the query location.
     * @return less than 0 if there are no overlapping targets at {@code location}.
     * @throws IllegalArgumentException   if {@code location} is {@code null}.
     * @throws AmbiguousTargetException if the query location overlaps more than one target.
     */
    default int index(final Locatable location) {
        final IndexRange range = indexRange(location);
        switch (range.size()) {
            case 1:
                return range.from;
            case 0:
                return - (range.from + 1);
            default:
                throw new AmbiguousTargetException(
                        String.format("location '%s' overlaps with %d targets: from '%s' to '%s'.",
                                location,range.size(), target(range.from), target(range.to - 1)));
        }
    }

    /**
     * Run a task (lambda function) on the set of targets that overlap a location.
     *
     * <p>
     *     This method may result in an arbitrary unchecked exception when one is
     *     thrown by one of the task executions.
     * </p>
     *
     * @param overlapRegion the query region that covers all targets to be processed.
     * @param consumer the task to be run on each overlapped target.
     * @throws IllegalArgumentException if any, {@code overlapRegion} or {@code targetTask}, is {@code null}.
     */
    default void forEachTarget(final Locatable overlapRegion, final IndexedTargetConsumer<E> consumer) {
        if (consumer == null) {
            throw new IllegalArgumentException("the indexed target consumer cannot be null");
        }
        final IndexRange indexRange = indexRange(overlapRegion);
        for (int i = indexRange.from; i < indexRange.to; i++) {
            consumer.accept(i, target(i));
        }
    }

    /**
     * Returns the total size of this collection in base-pairs (computed as the sum of all its targets sizes).
     *
     * @return 0 or greater.
     */
    default long totalSize() {
        return IntStream.range(0, targetCount())
                .map(i -> location(i).size()).sum();
    }

    /**
     * Returns the genome location of the target given its numeric index in the collection.
     *
     * @param index the target index.
     * @return never {@code null}.
     * @throws java.lang.IndexOutOfBoundsException if {@code index} is not valid.
     */
    SimpleInterval location(final int index);

    /**
     * Returns the genome location of a target.
     * <p>
     * It guarantees to return the location associated with the target in this collection if indeed it is included in
     * it.
     * </p>
     * <p>
     * Otherwise a plausible location might be returned based on the {@code target} value alone or {@code null}.
     * </p>
     *
     * @param target a target.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code target} is {@code null}.
     */
    SimpleInterval location(final E target);

    /**
     * Returns targets that overlap a genomic locations.
     * <p>
     * All targets that overlap the input {@code location} by at least one base are returned</p>
     *
     * <p>The returned collection
     * might not be modifiable and might change if this {@link TargetCollection} is mutable.
     *
     * @param overlapRegion the query location.
     * @return never {@code null} but an empty list if there are no such targets. The returned list is
     *  unmodifiable.
     * @throws IllegalArgumentException if {@code overlapRegion} is {@code null}.
     */
    default List<E> targets(final Locatable overlapRegion) {
        final IndexRange range = indexRange(overlapRegion);
        return targets().subList(range.from,range.to);
    }

    /**
     * Returns targets within an index interval.
     *
     * @param range the target index range.
     * @return never {@code null} but perhaps an empty list. The returned list is umodifiable.
     * @throws IllegalArgumentException if {@code range} is {@code null} or
     *   {@code range} contains invalid indices in this target collection.
     */
    default List<E> targets(final IndexRange range) {
        if (range == null) {
            throw new IllegalArgumentException("range cannot be null");
        } else if (!range.isValidLength(targetCount())) {
            throw new IllegalArgumentException(
                    String.format("index range '%s' is invalid for a targetCount of %d.", range, targetCount()));
        }
        return targets().subList(range.from, range.to);
    }

    /**
     * Returns targets within an index interval.
     *
     * @param from the index of the first requested target.
     * @param to   the index of the target right after the last one requested.
     * @return never {@code null} but perhaps an empty list.
     * @throws IllegalArgumentException if {@code from} and {@code to} is not a valid target index range.
     */
    default List<E> targets(final int from, final int to) {
        if (from < 0) {
            throw new IllegalArgumentException(
                    String.format("from index (%d) cannot be negative",from));
        } else if (to > targetCount()) {
            throw new IllegalArgumentException(
                    String.format("to index (%d) cannot larger than the target-count (%d)",to, targetCount()));
        }
        return targets().subList(from,to);
    }

    /**
     * Returns index range of targets that overlap a location.
     * <p>
     * If no target overlap the location, the returned range would have size 0 and would indicate
     * the insert index where such targets would be found if they were part of the collection.
     * </p>
     *
     * @param location the target location.
     * @throws IllegalArgumentException if {@code location} is {@code null}.
     * @return never {@code null} but perhaps an empty range.
     */
    IndexRange indexRange(final Locatable location);

    /**
     * Returns a list of all targets sorted by genomic coordinates.
     *
     * @return never {@code null} but perhaps empty if there is no target in this collection.
     *   The returned list is unmodifiable.
     */
    List<E> targets();

    /**
     * Common interface for lambda functions that perform
     * operations per target.
     *
     * @param <E> the target type.
     */
    @FunctionalInterface
    interface IndexedTargetConsumer<E> {

        /**
         * Performs the corresponding operation on a target
         *
         * @param index the index of the target.
         * @param target  the target value.
         */
        void accept(final int index, final E target);
    }

    /**
     * Indicates that more than one target was found when at most one was expected.
     */
    final class AmbiguousTargetException extends GATKException {

        private static final long serialVersionUID = -11123131L;

        /**
         * Creates an exception object given an explanatory message.
         *
         * @param message explanatory message.
         */
        public AmbiguousTargetException(final String message) {
            super(message);
        }
    }
}
