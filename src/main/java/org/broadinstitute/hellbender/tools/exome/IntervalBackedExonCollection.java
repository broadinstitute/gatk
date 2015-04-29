package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Exon data-base based on a list of intervals.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class IntervalBackedExonCollection implements ExonCollection<SimpleInterval> {

    /**
     * Map from interval name to interval object.
     */
    private final Map<String,SimpleInterval> intervalsByName;

    /**
     * Sorted list of intervals.
     */
    private final List<SimpleInterval> sortedIntervals;

    /**
     * Cached index of the last overlapping interval found by
     * {@link #cachedBinarySearch(SimpleInterval)}.
     *
     * <p>
     *     This results in a noticeable performance gain when processing data in
     *     genomic order.
     * </p>
     */
    private int lastBinarySearchResult = -1;

/**
     * Creates a exon data-base give a sorted list of intervals.
     *
     * <p>
     *     Intervals will be sorted by their contig id lexicographical order.
     * </p>
     *
     * @param intervals the input interval list. Is assumed to be sorted
     * @throws IllegalArgumentException if {@code intervals} is {@code null}.
     */
    public IntervalBackedExonCollection(final List<SimpleInterval> intervals) {
        if (intervals == null) {
            throw new IllegalArgumentException("the input intervals cannot be null");
        }
        if (intervals.contains(null)) {
            throw new IllegalArgumentException("the input cannot contain null");
        }
        sortedIntervals = intervals.stream().sorted(SimpleInterval.LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList());
        checkForOverlaps(sortedIntervals);
        this.intervalsByName = composeIntervalsByName(sortedIntervals);
    }

    /**
     * Fails with an exception if the intervals collection has overlapping intervals.
     * @param sortedIntervals the intervals sorted.
     */
    private static void checkForOverlaps(final List<SimpleInterval> sortedIntervals) {
        final OptionalInt failureIndex = IntStream.range(1, sortedIntervals.size())
                .filter(i -> sortedIntervals.get(i-1).overlaps(sortedIntervals.get(i)))
                .findFirst();

        if (failureIndex.isPresent()) {
            final int index = failureIndex.getAsInt();
            throw new IllegalArgumentException(
                    String.format("input intervals contain at least two overlapping intervals: %s and %s",
                            sortedIntervals.get(index-1),sortedIntervals.get(index)));
        }
    }

    /**
     * Composes a map from name to interval.
     *
     * <p>
     *     Also check whether there are more than one interval with the same name and in that case
     *     it fails with an exception.
     * </p>
     *
     * @param intervals the input intervals.
     * @throws IllegalArgumentException if {@code intervals} is {@code null} or it contains {@code nulls},
     *                  or is not a coherent interval list as per the criteria above.
     */
    private Map<String,SimpleInterval> composeIntervalsByName(final List<SimpleInterval> intervals) {
        final Map<String,SimpleInterval> result = new HashMap<>(intervals.size());
        for (final SimpleInterval location : intervals) {
            final String name = name(location);
            if (name != null) {
                final SimpleInterval previous = result.put(name,location);
                if (previous != null) {
                    throw new IllegalStateException(
                            String.format("more than one interval in the input list results in the same name (%s); " +
                                            "perhaps repeated: '%s' and '%s'.",
                                    name,previous,location));
                }
            }
        }
        return result;
    }

    /**
     * Produces the unique name for a genomic interval.
     *
     * <p>This name can depend in either the index of the interval, its location or both</p>
     *
     * <p>
     *     The implementation of this method is such that the name will be unique, i.e. that two different locations
     *     won't ever result in the same name.
     * </p>
     *
     * @param location the genome location of the interval. Is assumed not to be null.
     *
     * @return can return {@code null} indicating that this interval should remain anonymous.
     */
    @Override
    public String name(final SimpleInterval location) {
        return location.toString();
    }

    @Override
    public List<SimpleInterval> exons() {
        return sortedIntervals;
    }

    @Override
    public int exonCount() {
        return sortedIntervals.size();
    }

    @Override
    public SimpleInterval exon(final int index) {
        if (index < 0) {
            throw new IllegalArgumentException("index cannot be negative");
        } else if (index >= sortedIntervals.size()) {
            throw new IllegalArgumentException(
                    String.format("index (%d) cannot be equal or greater than the exon-count (%d)",index,exonCount()));
        }
        return sortedIntervals.get(index);
    }

    @Override
    public SimpleInterval exon(final String name) {
        if (name == null) {
            throw new IllegalArgumentException("the input name cannot be null");
        }
        return intervalsByName.get(name);
    }

    @Override
    public int index(final String name) {
        final SimpleInterval exon = intervalsByName.get(name);
        if (exon == null) {
            return -1;
        } else {
            final int searchIndex = uncachedBinarySearch(exon);
            if (searchIndex < 0) { // checking just in case.
                throw new IllegalStateException("could not found named interval amongst sorted intervals, impossible");
            }
            return searchIndex;
        }
    }

    @Override
    public IndexRange indexRange(final SimpleInterval location) {
        if (location == null) {
            throw new IllegalArgumentException("the input location cannot be null");
        } else {
            final int searchIndex = cachedBinarySearch(location);
            if (searchIndex < 0) {
                return new IndexRange(-searchIndex - 1, -searchIndex - 1);
            } else {
                final int firstOverlappingIndex = extendSearchIndexBackwards(location, searchIndex);
                final int lastOverlappingIndex = extendSearchIndexForward(location, searchIndex);
                return new IndexRange(firstOverlappingIndex, lastOverlappingIndex + 1);
            }
        }
    }

    @Override
    public SimpleInterval location(final SimpleInterval exon) {
        if (exon == null) {
            throw new IllegalArgumentException("the exon cannot be null");
        }
        return exon;
    }

    @Override
    public SimpleInterval location(final int index) {
        if (index < 0 || index >= sortedIntervals.size()) {
            throw new IllegalArgumentException(
                    String.format("the index provided, %d, is not within the valid range [%d,%d).",index,0,sortedIntervals.size()));
        }
        return sortedIntervals.get(index);
    }

    @Override
    public SimpleInterval exon(final SimpleInterval overlapRegion) {
        final int searchIndex = index(overlapRegion);
        return searchIndex < 0 ? null : sortedIntervals.get(searchIndex);
    }

    @Override
    public int index(final SimpleInterval location) {
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

    @Override
    public List<SimpleInterval> exons(final SimpleInterval overlapRegion) {
        final IndexRange range = indexRange(overlapRegion);
        return sortedIntervals.subList(range.from, range.to);
    }

    /**
     * Implements a cached binary search of the overlapping intervals.
     *
     * <p>
     *     This was found to improve performance significantly when analyzing empirical
     *     exome data in sequential order.
     * </p>
     *
     * <p>
     *     First it checks
     *     whether the query interval overlaps with the result of the
     *     last search (whether a hit or a miss (taking the insertion position of the miss)).
     * </p>
     * <p>
     *     If the last search was a hit but the new query location does not overlap,
     *     it checks on the next interval and if not fails over to the regular binary search.
     * </p>
     *
     * <p>
     *     A positive (0 or greater) returned value indicates a hit where the corresponding interval
     *     is guaranteed to overlap the query {@code location}. There might be more intervals that
     *     overlap this location and they must all be contiguous to the returned index.
     * <p>
     *     One can use {@link #extendSearchIndexBackwards(SimpleInterval, int)}
     *     and {@link #extendSearchIndexForward(SimpleInterval, int)}
     *     to find the actual overlapping index range.
     * </p>
     *
     * <p>
     *     In contrast, a negative result indicates that there is no overlapping interval.
     *     This value encodes the insertion position for the query location; that is, where would
     *     the query location be inserted were it to be added to the sorted list of intervals:
     * </p>
     * <p>
     *     <code>insertion index == -(result + 1)</code>
     * </p>
     *
     * @param location the query location.
     * @return any integer between <code>-{@link #exonCount()}-1</code> and <code>{@link #exonCount()} - 1</code>.
     */
    private int cachedBinarySearch(final SimpleInterval location) {

        if (lastBinarySearchResult < 0) {
            final int candidate = -(lastBinarySearchResult + 1);
            if (candidate >= sortedIntervals.size()) {
                return lastBinarySearchResult = uncachedBinarySearch(location);
            } else if (sortedIntervals.get(candidate).overlaps(location)) {
                return lastBinarySearchResult = candidate;
            } else {
                return lastBinarySearchResult = uncachedBinarySearch(location);
            }
        } else {
            if (sortedIntervals.get(lastBinarySearchResult).overlaps(location)) {
                return lastBinarySearchResult;
            } else {
                final int candidate = lastBinarySearchResult + 1;
                if (candidate == sortedIntervals.size()) {
                    return lastBinarySearchResult = uncachedBinarySearch(location);
                } else if (sortedIntervals.get(candidate).overlaps(location)) {
                    return lastBinarySearchResult = candidate;
                } else {
                    return lastBinarySearchResult = uncachedBinarySearch(location);
                }
            }
        }
    }

    /**
     * Implements a binary search of the overlapping intervals.
     *
     * <p>
     *     A positive (0 or greater) returned value indicates a hit where the corresponding interval
     *     is guaranteed to overlap the query {@code location}. There might be more intervals that
     *     overlap this location and they must all be contiguous to the returned index.
     * <p>
     *     One can use {@link #extendSearchIndexBackwards(SimpleInterval, int)}
     *     and {@link #extendSearchIndexForward(SimpleInterval, int)}
     *     to find the actual overlapping index range.
     * </p>
     *
     * <p>
     *     In contrast, a negative result indicates that there is no overlapping interval.
     *     This value encodes the insertion position for the query location; that is, where would
     *     the query location be inserted were it to be added to the sorted list of intervals:
     * </p>
     * <p>
     *     <code>insertion index == -(result + 1)</code>
     * </p>
     *
     * @param location the query location.
     * @return any integer between <code>-{@link #exonCount()}-1</code> and <code>{@link #exonCount()} - 1</code>.
     */
    private int uncachedBinarySearch(final SimpleInterval location) {
        if (sortedIntervals.size() == 0) {
            return -1;
        }

        final int searchResult = Collections.binarySearch(sortedIntervals, location, SimpleInterval.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        if (searchResult >= 0) {
            return searchResult;
        } else {
            final int insertIndex = - (searchResult + 1);
            if (insertIndex < sortedIntervals.size() && sortedIntervals.get(insertIndex).overlaps(location)) {
                return insertIndex;
            } if (insertIndex > 0 && sortedIntervals.get(insertIndex - 1).overlaps(location)) {
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
    private int extendSearchIndexForward(final SimpleInterval location, final int startIndex) {
        final ListIterator<SimpleInterval> it = sortedIntervals.listIterator(startIndex + 1);
        while (it.hasNext()) {
            final SimpleInterval next = it.next();
            if (!location.overlaps(next)) {
                return it.previousIndex() - 1;
            }
        }
        return it.previousIndex();
    }

    /**
     * Looks for the first index in {@link #sortedIntervals} that has an overlap with the input {@code location}
     * starting at {@code startIndex} and assuming that the element at that index has an overlap with {@code location}.
     */
    private int extendSearchIndexBackwards(final SimpleInterval location, final int startIndex) {
        final ListIterator<SimpleInterval> it = sortedIntervals.listIterator(startIndex);
        while (it.hasPrevious()) {
            final SimpleInterval previous = it.previous();
            if (!location.overlaps(previous)) {
                return it.nextIndex() + 1;
            }
        }
        return it.nextIndex();
    }

}
