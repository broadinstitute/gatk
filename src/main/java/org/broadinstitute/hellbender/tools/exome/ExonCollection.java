package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;
import java.util.stream.IntStream;

/**
 * Exon collection.
 * <p>
 *      This interface includes operations to query exons based on genomic coordinates and their names (if provided).
 * </p>
 *
 * <p>
 *      Additional per exon meta-data can be added by implementing class through the customizable type parameter
 *      &lt;E&gt;.
 * </p>
 *
 * @param <E> exon meta-data type.
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public interface ExonCollection<E> {

    /**
     * Returns the number of exons in this Exome.
     *
     * @return 0 or greater.
     */
    int exonCount();

    /**
     * Returns a exon given is index in the exon.
     * <p>
     * Indices are 0-based, thus the first exon has index 0, the second 1 and so forth.
     * </p>
     *
     * @param index the query index.
     * @throws IndexOutOfBoundsException if {@code index} is not within
     *                                   valid bounds <code>[0 .. ({@link #exonCount()})</code>.
     */
    E exon(final int index);

    /**
     * Returns the name of an exon.
     *
     * <p>
     *     This method is guaranteed to return the name of the exon in this collection as long as it is indeed
     *     included in it.
     * </p>
     *
     * <p>
     *     Otherwise it should return a reasonable name for it or {@code null}.
     * </p>
     *
     * @param exon target exon.
     *
     * @throws IllegalArgumentException if {@code exon} is {@code null}.
     */
    String name(final E exon);

    /**
     * Returns the exon given its name.
     *
     * @param name the query exon name.
     * @return {@code null} if there is not such an exon.
     */
    default E exon(final String name) {
        final int index = index(name);
        if (index < 0) {
            return null;
        } else {
            return exon(index);
        }
    }

    /**
     * Returns the exon index given its name.
     *
     * @param name the query exon name.
     * @throws IllegalArgumentException if {@code name} is {@code null}.
     * @return {@code -1} if there is not such an exon.
     */
    int index(final String name);

    /**
     * Returns the exon that overlap a particular location.
     * <p>
     *     The returned exon coordinates does not need to match exactly the
     * the query {@code location} just overlap it by at least one base</p>.
     * <p>
     *     This method will return {@code null} if there is no such a exon.
     * </p>
     * <p>
     *     If there is more than on exon that overlap with the given region it will
     *     throw a {@link AmbiguousExonException} instead.
     * </p>
     *
     * @param overlapRegion the query location.
     * @return {@code null} if there is no overlapping exons at {@code overlapRegion}.
     * @throws IllegalArgumentException if {@code overlapRegion} is {@code null}.
     * @throws AmbiguousExonException if the query location overlaps more than one exon.
     */
    default E exon(final SimpleInterval overlapRegion) {
        final int index = index(overlapRegion);
        return index < 0 ? null : exon(index);
    }

    /**
     * Returns the index of the exon that overlap a particular location.
     * <p>
     *    The returned index's exon coordinates does not need to match
     *    exactly the query {@code location} just overlap it by at
     *    least one base.
     * </p>
     * <p>
     *    This method will return a negative value if there is no such a exon.
     *    This value can be transformed to recover the position where the query location would be inserted
     *    if it was present as <code>-(index(loc) + 1)</code>.
     * </p>
     *
     * <p>
     *    This method will result in a {@link AmbiguousExonException} if there is more than one
     *    exon that overlaps that query {@code location}.
     * </p>
     *
     * @param location the query location.
     * @return less than 0 if there is no overlapping exons at {@code location}.
     * @throws IllegalArgumentException   if {@code location} is {@code null}.
     * @throws AmbiguousExonException if the query location overlaps more than one exon.
     */
    default int index(final SimpleInterval location) {
        final IndexRange range = indexRange(location);
        switch (range.size()) {
            case 1:
                return range.from;
            case 0:
                return - (range.from + 1);
            default:
                throw new AmbiguousExonException(
                        String.format("location '%s' overlaps with %d exons: from '%s' to '%s'.",
                                location,range.size(),exon(range.from),exon(range.to - 1)));
        }
    }

    /**
     * Run a task (lambda function) on the set of exons that overlap a location.
     *
     * <p>
     *     This method may result in an arbitrary unchecked exception when one is
     *     thrown by one of the task executions.
     * </p>
     *
     * @param overlapRegion the query region that covers all exons to be processed.
     * @param consumer the task to be run on each overlapped exon.
     * @throws IllegalArgumentException if any, {@code overlapRegion} or {@code exonTask}, is {@code null}.
     */
    default void forEachExon(final SimpleInterval overlapRegion, final IndexedExonConsumer<E> consumer) {
        if (consumer == null) {
            throw new IllegalArgumentException("the indexed exon consumer cannot be null");
        }
        final IndexRange indexRange = indexRange(overlapRegion);
        for (int i = indexRange.from; i < indexRange.to; i++) {
            consumer.accept(i, exon(i));
        }
    }

    /**
     * Returns the size of the exome in base-pairs (computed as the sum of all its exons sizes).
     *
     * @return 0 or greater.
     */
    default long exomeSize() {
        return IntStream.range(0,exonCount())
                .map(i -> location(i).size()).sum();
    }

    /**
     * Returns the genome location of the exon given its numeric index in the collection.
     *
     * @param index the target exon index.
     * @return never {@code null}.
     * @throws java.lang.IndexOutOfBoundsException if {@code index} is not valid.
     */
    SimpleInterval location(final int index);

    /**
     * Returns the genome location of an exon.
     * <p>
     * It guarantees to return the location associated with the exon in this collection if indeed it is included in
     * it.
     * </p>
     * <p>
     * Otherwise a plausible location might be returned based on the {@code exon} value alone or {@code null}.
     * </p>
     *
     * @param exon the target exon.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code exon} is {@code null}.
     */
    SimpleInterval location(final E exon);

    /**
     * Returns exons that overlap a genomic locations.
     * <p>
     * All exons that overlap the input {@code location} by at least one base are returned</p>
     *
     * <p>The returned collection
     * might not be modifiable and might change if this {@link ExonCollection} is mutable.
     *
     * @param overlapRegion the query location.
     * @return never {@code null} but an empty list if there is no such exons. The returned list is
     *  unmodifiable.
     * @throws IllegalArgumentException if {@code overlapRegion} is {@code null}.
     */
    default List<E> exons(final SimpleInterval overlapRegion) {
        final IndexRange range = indexRange(overlapRegion);
        return exons().subList(range.from,range.to);
    }

    /**
     * Returns exons within an index interval.
     *
     * @param range the target index range.
     * @return never {@code null} but perhaps an empty list. The returned list is umodifiable.
     * @throws IllegalArgumentException if {@code range} is {@code null} or
     *   {@code range} contains invalid indices in this exon collection.
     */
    default List<E> exons(final IndexRange range) {
        if (range == null) {
            throw new IllegalArgumentException("range cannot be null");
        } else if (!range.isValidLength(exonCount())) {
            throw new IllegalArgumentException(
                    String.format("index range '%s' is invalid for an exonCount of %d.", range, exonCount()));
        }
        return exons().subList(range.from,range.to);
    }


    /**
     * Returns exons within an index interval.
     *
     * @param from the index of the first requested exon.
     * @param to   the index of the exon right after the last one requested.
     * @return never {@code null} but perhaps an empty list.
     * @throws IllegalArgumentException if {@code from} and {@code to} is not a valid exon index range.
     */
    default List<E> exons(final int from, final int to) {
        if (from < 0) {
            throw new IllegalArgumentException(
                    String.format("from index (%d) cannot be negative",from));
        } else if (to > exonCount()) {
            throw new IllegalArgumentException(
                    String.format("to index (%d) cannot larger than the exon-count (%d)",to,exonCount()));
        }
        return exons().subList(from,to);
    }

    /**
     * Returns index range of exons that overlap a location.
     * <p>
     * If no exon overlap the location, the returned range would have size 0 and would indicate
     * the insert index where such exons would be found if they were part of the collection.
     * </p>
     *
     * @param location the target location.
     * @throws IllegalArgumentException if {@code location} is {@code null}.
     * @return never {@code null} but perhaps an empty range.
     */
    IndexRange indexRange(final SimpleInterval location);

    /**
     * Returns a list of all exons sorted by genomic coordinates.
     *
     * @return never {@code null} but perhaps empty if there is no exon in this collection.
     *   The returned list is unmodifiable.
     */
    List<E> exons();

    /**
     * Common interface for lambda functions that perform
     * operations per exon.
     *
     * @param <E> the exon type.
     */
    @FunctionalInterface
    interface IndexedExonConsumer<E> {

        /**
         * Performs the corresponding operation on a exon
         *
         * @param index the index of the target exon.
         * @param exon  the target exon value.
         */
        void accept(final int index, final E exon);
    }

    /**
     * Indicates that more than one exon was found when at most one was expected.
     */
    final class AmbiguousExonException extends GATKException {

        private static final long serialVersionUID = -11123131L;

        /**
         * Creates an exception object given an explanatory message.
         *
         * @param message explanatory message.
         */
        public AmbiguousExonException(final String message) {
            super(message);
        }
    }
}
