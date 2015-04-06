package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.IndexRange;

import java.util.*;

/**
 * Exon data-base based on a list of intervals.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class IntervalBackedExonCollection implements ExonCollection<GenomeLoc> {

    /**
     * Map from interval name to interval object.
     */
    private final Map<String,GenomeLoc> intervalsByName;

    /**
     * Sorted list of intervals.
     */
    private final List<GenomeLoc> sortedIntervals;

    /**
     * Cached index of the last overlapping interval found by
     * {@link #cachedBinarySearch(GenomeLoc)}.
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
     * @param intervals the input interval list. Is assumed to be sorted
     * @throws IllegalArgumentException if {@code intervals} is {@code null}.
     */
    public IntervalBackedExonCollection(final List<GenomeLoc> intervals) {
        if (intervals == null) {
            throw new IllegalArgumentException("the input intervals cannot be null");
        }
        checkInputIntervalCoherence(intervals);
        this.sortedIntervals = Collections.unmodifiableList(new ArrayList<>(intervals));
        this.intervalsByName = composeIntervalsByName(intervals);
    }

    /**
     * Throws a exception if there some interval location incoherence.
     *
     * <p>
     *     These include:
     *     <ul>
     *         <li>{@code} null intervals,</li>
     *         <li>intervals are out of order,</li>
     *         <li>intervals don't overlap each other,</li>
     *         <li>some interval is unmapped or</li>
     *         <li>is the special value {@link GenomeLoc#WHOLE_GENOME}</li>
     *     </ul>
     * </p>
     *
     * @param intervals the input intervals.
     * @throws IllegalArgumentException if there is any inconsistency in {@code intervals}.
     */
    private void checkInputIntervalCoherence(final List<GenomeLoc> intervals) {
        GenomeLoc lastLocation = null;
        for (final GenomeLoc location : intervals) {
            checkInputIntervalCoherence(lastLocation, location);
            lastLocation = location;
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
    private Map<String,GenomeLoc> composeIntervalsByName(final List<GenomeLoc> intervals) {
        final Map<String,GenomeLoc> result = new HashMap<>(intervals.size());
        for (final GenomeLoc location : intervals) {
            final String name = intervalName(location);
            if (name != null) {
                final GenomeLoc previous = result.put(name,location);
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
     * Throws a exception if there some interval location incoherence.
     *
     * See {@link #checkInputIntervalCoherence(List)} documenation for details.
     *
     * @param lastLocation previous location.
     * @param nextLocation next location to add.
     *
     * @throws IllegalArgumentException if there is some problem with a location given the previous one.
     */
    private static void checkInputIntervalCoherence(final GenomeLoc lastLocation, final GenomeLoc nextLocation) {
        if (nextLocation == null) {
            throw new IllegalArgumentException("the input intervals contain a null location");
        } else if (nextLocation.isUnmapped()) {
            throw new IllegalArgumentException("the input intervals contain unmapped locations");
        } else if (nextLocation == GenomeLoc.WHOLE_GENOME) {
            throw new IllegalArgumentException("the input intervals contain a whole-genome location");
        } else if (lastLocation != null && !lastLocation.isBefore(nextLocation)) {
            throw new IllegalArgumentException(
                    String.format("the input intervals list is not in order or " +
                            "contains overlapped intervals: '%s' is not strictly before '%s'",lastLocation,nextLocation));
        }
    }

    /**
     * Produces the unique name for a genomic interval.
     *
     * <p>This name can depend in either the index of the interval, its location or both</p>
     *
     * <p>The implementation can assume that input interval is a proper one with start and stop positions and a contig,
     * thus the result is unpredictable with special intervals values such as {@link GenomeLoc#UNMAPPED}
     *   or {@link GenomeLoc#WHOLE_GENOME}.
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
    protected String intervalName(final GenomeLoc location) {
        return location.toString();
    }

    @Override
    public List<GenomeLoc> exons() {
        return sortedIntervals;
    }

    @Override
    public int exonCount() {
        return sortedIntervals.size();
    }

    @Override
    public GenomeLoc exon(final int index) {
        if (index < 0) {
            throw new IllegalArgumentException("index cannot be negative");
        } else if (index >= sortedIntervals.size()) {
            throw new IllegalArgumentException(
                    String.format("index (%d) cannot be equal or greater than the exon-count (%d)",index,exonCount()));
        }
        return sortedIntervals.get(index);
    }

    @Override
    public GenomeLoc exon(final String name) {
        if (name == null) {
            throw new IllegalArgumentException("the input name cannot be null");
        }
        return intervalsByName.get(name);
    }

    @Override
    public int index(final String name) {
        final GenomeLoc exon = intervalsByName.get(name);
        if (exon == null) {
            return -1;
        } else {
            final int searchIndex = Collections.binarySearch(sortedIntervals, exon);
            if (searchIndex < 0) { // checking just in case.
                throw new IllegalStateException("could not found named interval amongst sorted intervals, impossible");
            }
            return searchIndex;
        }
    }

    @Override
    public IndexRange indexRange(final GenomeLoc location) {
        if (location == null) {
            throw new IllegalArgumentException("the input location cannot be null");
        } else if (location == GenomeLoc.WHOLE_GENOME) {
            return new IndexRange(0, sortedIntervals.size());
        } else if (location.isUnmapped()) {
            return new IndexRange(sortedIntervals.size(),sortedIntervals.size());
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
    public GenomeLoc location(final GenomeLoc exon) {
        if (exon == null) {
            throw new IllegalArgumentException("the exon cannot be null");
        }
        return exon;
    }

    @Override
    public GenomeLoc location(final int index) {
        if (index < 0 || index >= sortedIntervals.size()) {
            throw new IllegalArgumentException(
                    String.format("the index provided, %d, is not within the valid range [%d,%d).",index,0,sortedIntervals.size()));
        }
        return sortedIntervals.get(index);
    }

    /**
     * Comparator used to search overlapping genomic locations.
     *
     * <p>
     *     In contrast to the default sorting of {@link GenomeLoc} instances, this comparator
     *     will return 0 of there is any base overlap.
     * </p>
     */
    private static final Comparator<GenomeLoc> OVERLAP_SEARCH_COMPARATOR = (left, right) -> {
        if (left.isUnmapped()) {
            return right.isUnmapped() ? 0 : 1;
        } else if (right.isUnmapped()) {
            return -1;
        } else if (left.isBefore(right)) {
            return -1;
        } else if (left.isPast(right)) {
            return 1;
        } else {
            return 0;
        }
    };

    @Override
    public GenomeLoc exon(final GenomeLoc overlapRegion) {
        final int searchIndex = index(overlapRegion);
        return searchIndex < 0 ? null : sortedIntervals.get(searchIndex);
    }

    @Override
    public int index(final GenomeLoc location) {
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
    public List<GenomeLoc> exons(final GenomeLoc overlapRegion) {
        final IndexRange range = indexRange(overlapRegion);
        return sortedIntervals.subList(range.from, range.to);
    }

    /**
     * Implement a cached binary search of the overlapping intervals.
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
     *     One can use {@link #extendSearchIndexBackwards(GenomeLoc, int)}
     *     and {@link #extendSearchIndexForward(GenomeLoc, int)}
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
     * <p>
     *     The input location is assumed to be proper location thus fully mapped and never the special
     *     value {@link GenomeLoc#WHOLE_GENOME}. Otherwise the outcome is unspecified.
     * </p>
     *
     * @param location the query location.
     * @return any integer between <code>-{@link #exonCount()}-1</code> and <code>{@link #exonCount()} - 1</code>.
     */
    private int cachedBinarySearch(final GenomeLoc location) {

        if (lastBinarySearchResult < 0) {
            final int candidate = -(lastBinarySearchResult + 1);
            if (candidate >= sortedIntervals.size()) {
                return lastBinarySearchResult = Collections.binarySearch(sortedIntervals, location, OVERLAP_SEARCH_COMPARATOR);
            } else if (sortedIntervals.get(candidate).overlapsP(location)) {
                return lastBinarySearchResult = candidate;
            } else {
                return lastBinarySearchResult = Collections.binarySearch(sortedIntervals, location, OVERLAP_SEARCH_COMPARATOR);
            }
        } else {
            if (sortedIntervals.get(lastBinarySearchResult).overlapsP(location)) {
                return lastBinarySearchResult;
            } else {
                final int candidate = lastBinarySearchResult + 1;
                if (candidate == sortedIntervals.size()) {
                    return lastBinarySearchResult = Collections.binarySearch(sortedIntervals, location, OVERLAP_SEARCH_COMPARATOR);
                } else if (sortedIntervals.get(candidate).overlapsP(location)) {
                    return lastBinarySearchResult = candidate;
                } else {
                    return lastBinarySearchResult = Collections.binarySearch(sortedIntervals, location, OVERLAP_SEARCH_COMPARATOR);
                }
            }
        }
    }

    /**
     * Looks for the last index in {@link #sortedIntervals} that has an overlap with the input {@code location}.
     * starting at {@code startIndex} and assuming that the element at that index has an
     * overlap with {@code location}.
     */
    private int extendSearchIndexForward(final GenomeLoc location, final int startIndex) {
        final ListIterator<GenomeLoc> it = sortedIntervals.listIterator(startIndex + 1);
        while (it.hasNext()) {
            final GenomeLoc next = it.next();
            if (next.isPast(location)) {
                return it.previousIndex() - 1;
            }
        }
        return it.previousIndex();
    }

    /**
     * Looks for the first index in {@link #sortedIntervals} that has an overlap with the input {@code location}
     * starting at {@code startIndex} and assuming that the element at that index has an overlap with {@code location}.
     */
    private int extendSearchIndexBackwards(final GenomeLoc location, final int startIndex) {
        final ListIterator<GenomeLoc> it = sortedIntervals.listIterator(startIndex);
        while (it.hasPrevious()) {
            final GenomeLoc previous = it.previous();
            if (previous.isBefore(location)) {
                return it.nextIndex() + 1;
            }
        }
        return it.nextIndex();
    }

}
