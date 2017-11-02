package org.broadinstitute.hellbender.tools.copynumber.utils;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Collection of {@link Locatable}s that supports cached binary search of {@link Locatable} objects. In particular,
 * this class provides an interface to get a range of {@link Locatable}s that intersects with a given {@link Locatable}
 * query location.
 *
 * The search is optimized based on the assumption that the incoming queries of {@link Locatable}s are in sorted order.
 *
 * @param <E> a locatable object
 */
public class CachedBinarySearchIntervalList<E extends Locatable> {

    private List<E> sortedIntervals;

    private final Comparator<Locatable> intervalComparator = IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR;

    private int lastBinarySearchPosition = -1;

    public CachedBinarySearchIntervalList(final List<E> unsortedIntervalList) {
        Utils.nonEmpty(unsortedIntervalList, "Interval list cannot be empty");
        Utils.containsNoNull(unsortedIntervalList, "intervals may not be null");
        sortedIntervals = unsortedIntervalList.stream().sorted(intervalComparator).collect(Collectors.toList());
        checkForIntervalOverlaps(sortedIntervals);
    }

    private void checkForIntervalOverlaps(final List<E> sortedIntervals) {
        final OptionalInt failureIndex = IntStream.range(1, sortedIntervals.size())
                .filter(i -> IntervalUtils.overlaps(sortedIntervals.get(i-1), sortedIntervals.get(i)))
                .findFirst();

        if (failureIndex.isPresent()) {
            final int index = failureIndex.getAsInt();
            throw new IllegalArgumentException(
                    String.format("input intervals contain at least two overlapping intervals: %s and %s",
                            sortedIntervals.get(index-1), sortedIntervals.get(index)));
        }
    }

    /**
     * Returns index range of {@link Locatable}s that overlap a query location.
     * <p>
     * If no {@link Locatable} in the collection overlaps the location, the returned range would have size 0 and would indicate
     * the insert index where such {@link Locatable} would be found if they were part of the collection.
     * </p>
     *
     * @param location query location
     * @throws IllegalArgumentException if {@code location} is {@code null}.
     * @return never {@code null} but perhaps an empty range.
     */
    public IndexRange findIntersectionRange(final Locatable location) {
        Utils.nonNull(location, "the input location cannot be null");
        final int searchIndex = cachedBinarySearch(location);
        if (searchIndex < 0) {
            return new IndexRange(-searchIndex - 1, -searchIndex - 1);
        } else {
            final int firstOverlappingIndex = extendSearchIndexBackwards(location, searchIndex);
            final int lastOverlappingIndex = extendSearchIndexForward(location, searchIndex);
            return new IndexRange(firstOverlappingIndex, lastOverlappingIndex + 1);
        }
    }

    /**
     * Get list of sorted intervals
     *
     * @return sorted intervals
     */
    public List<E> getSortedIntervals() {
        return sortedIntervals;
    }

    /**
     * Check the position from the last search for overlap, if it's unsuccessful revert to binary search
     */
    private int cachedBinarySearch(final Locatable location) {
        if (lastBinarySearchPosition < 0) {
            final int candidate = -(lastBinarySearchPosition + 1);
            if (candidate >= sortedIntervals.size()) {
                return lastBinarySearchPosition = uncachedBinarySearch(location);
            } else if (IntervalUtils.overlaps(sortedIntervals.get(candidate),location)) {
                return lastBinarySearchPosition = candidate;
            } else {
                return lastBinarySearchPosition = uncachedBinarySearch(location);
            }
        } else {
            if (IntervalUtils.overlaps(sortedIntervals.get(lastBinarySearchPosition), location)) {
                return lastBinarySearchPosition;
            } else {
                final int candidate = lastBinarySearchPosition + 1;
                if (candidate == sortedIntervals.size()) {
                    return lastBinarySearchPosition = uncachedBinarySearch(location);
                } else if (IntervalUtils.overlaps(sortedIntervals.get(candidate), location)) {
                    return lastBinarySearchPosition = candidate;
                } else {
                    return lastBinarySearchPosition = uncachedBinarySearch(location);
                }
            }
        }
    }

    /**
     * Perform a binary search of the interval collection. Positive return value indicates a search hit, whereas negative
     * return value indicates where would query location be inserted were it to be added to the sorted list of intervals
     */
    private int uncachedBinarySearch(final Locatable location) {
        final int searchResult = Collections.binarySearch(sortedIntervals, location, intervalComparator);
        if (searchResult >= 0) {
            return searchResult;
        } else {
            final int insertIndex = - (searchResult + 1);
            if (insertIndex < sortedIntervals.size() && IntervalUtils.overlaps(sortedIntervals.get(insertIndex), location)) {
                return insertIndex;
            } if (insertIndex > 0 && IntervalUtils.overlaps(sortedIntervals.get(insertIndex - 1), location)) {
                return insertIndex - 1;
            } else {
                return searchResult;
            }
        }
    }

    /**
     * Looks for the last index in {@link #sortedIntervals} that has an overlap with the input {@code location}.
     * starting at {@code startIndex} and assuming that the element at that index has an
     * overlap with {@code location}.
     */
    private int extendSearchIndexForward(final Locatable location, final int startIndex) {
        final ListIterator<E> it = sortedIntervals.listIterator(startIndex + 1);
        while (it.hasNext()) {
            final E next = it.next();
            if (!IntervalUtils.overlaps(location, next)) {
                return it.previousIndex() - 1;
            }
        }
        return it.previousIndex();
    }

    /**
     * Looks for the first index in {@link #sortedIntervals} that has an overlap with the input {@code location}
     * starting at {@code startIndex} and assuming that the element at that index has an overlap with {@code location}.
     */
    private int extendSearchIndexBackwards(final Locatable location, final int startIndex) {
        final ListIterator<E> it = sortedIntervals.listIterator(startIndex);
        while (it.hasPrevious()) {
            final E previous = it.previous();
            if (!IntervalUtils.overlaps(location, previous)) {
                return it.nextIndex() + 1;
            }
        }
        return it.nextIndex();
    }
}
