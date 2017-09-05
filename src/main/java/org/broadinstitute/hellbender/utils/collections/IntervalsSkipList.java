package org.broadinstitute.hellbender.utils.collections;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Iterators;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Holds many intervals in memory, with an efficient operation to get
 * intervals that overlap a given query interval.
 *
 * This version allows intervals to lie on different contigs.
 */
public final class IntervalsSkipList<T extends Locatable> extends AbstractCollection<T> implements Serializable {
    private static final long serialVersionUID = 1L;

    private final int size;
    private final Map<String, IntervalsSkipListOneContig<T>> intervals;

    /**
     * Creates an IntervalsSkipList that holds a copy of the given intervals, sorted
     * and indexed.
     *
     * @param loc Locatables, not necessarily sorted. Will be iterated over exactly once.
     */
    public IntervalsSkipList(final Iterable<T> loc) {
        final Map<String, List<T>> variantsPerContig = new LinkedHashMap<>();
        int size = 0;
        for (final T v : loc) {
            final String k = v.getContig();
            variantsPerContig.putIfAbsent(k, new ArrayList<>());
            variantsPerContig.get(k).add(v);
            size++;
        }
        intervals = new LinkedHashMap<>();
        for (String k : variantsPerContig.keySet()) {
            intervals.put(k, new IntervalsSkipListOneContig<>(variantsPerContig.get(k)));
        }
        this.size = size;
    }

    /**
     * Returns all the intervals that overlap with the query.
     * The query doesn't *have* to be in the same contig as any interval we
     * hold, but of course if it isn't you'll get an empty result.
     * You may modify the returned list.
     */
    public List<T> getOverlapping(final SimpleInterval query) {
        final String k = query.getContig();
        final IntervalsSkipListOneContig<T> result = intervals.get(k);
        if (result == null){
            return new ArrayList<>();
        }
        return result.getOverlapping(query);
    }

    @Override
    @SuppressWarnings("unchecked")
    public Iterator<T> iterator() {
        return Iterators.concat(intervals.values().stream().map(Iterable::iterator).toArray(i -> (Iterator<T>[]) new Iterator[i]));
    }

    @Override
    public int size() {
        return size;
    }

    /**
     * Returns the list of contigs for which this list has shards.
     * <p>
     *     The order of contigs in the returned list reflect the order these would be traverse by
     *     the iterator returned by {@link #iterator()}.
     * </p>
     * <p>
     *     The returned list cannot be modified.
     * </p>
     * @return never {@code null}.
     */
    public List<String> contigs() {
        return Collections.unmodifiableList(intervals.keySet().stream().collect(Collectors.toList()));
    }

    /**
     * Returns an iterator over the shards on a particular contig.
     * <p>
     *     Request on a (unknown) contig without shards will return an iterator without elements to
     *     traverse ({@code it.hasNext() == false}).
     * </p>
     * @param contig the target contig.
     * @throws IllegalArgumentException if the input contig is {@code null}.
     * @return never {@code null}.
     */
    public Iterator<T> contigIterator(final String contig) {
        Utils.nonNull(contig);
        final IntervalsSkipListOneContig<T> contigIntervals = intervals.get(contig);
        if (contigIntervals == null) {
            return ImmutableSet.<T>of().iterator();
        } else {
            return contigIntervals.iterator();
        }
    }

}
