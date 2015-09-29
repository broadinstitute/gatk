package org.broadinstitute.hellbender.utils.collections;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Holds many intervals in memory, with an efficient operation to get
 * intervals that overlap a given query interval.
 *
 * This version allows intervals to lie on different contigs.
 */
public final class IntervalsSkipList<T extends Locatable> implements Serializable {
    private static final long serialVersionUID = 1L;

    private final Map<String, IntervalsSkipListOneContig<T>> intervals;

    /**
     * Creates an IntervalsSkipList that holds a copy of the given intervals, sorted
     * and indexed.
     *
     * @param loc Locatables, not necessarily sorted. Will be iterated over exactly once.
     */
    public IntervalsSkipList(final Iterable<T> loc) {
        final Map<String,ArrayList<T>> variantsPerContig = new HashMap<>();
        for (T v : loc) {
            String k = v.getContig();
            ArrayList<T> l;
            if (variantsPerContig.containsKey(k)) {
                l = variantsPerContig.get(k);
            } else {
                l = new ArrayList<>();
                variantsPerContig.put(k,l);
            }
            l.add(v);
        }
        intervals = new HashMap<>();
        for (final String k : variantsPerContig.keySet()) {
            intervals.put(k, new IntervalsSkipListOneContig<>(variantsPerContig.get(k)));
        }
    }

    /**
     * Returns all the intervals that overlap with the query.
     * The query doesn't *have* to be in the same contig as any interval we
     * hold, but of course if it isn't you'll get an empty result.
     * You may modify the returned list.
     */
    public List<T> getOverlapping(final SimpleInterval query) {
        final String k = query.getContig();
        if (!intervals.containsKey(k)) return new ArrayList<>();
        return intervals.get(k).getOverlapping(query);
    }

}
