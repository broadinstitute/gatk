package org.broadinstitute.hellbender.utils.collections;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Holds many intervals in memory, with an efficient operation to get
 * intervals that overlap a given query interval.
 *
 * This version assumes that all the intervals lie on the same contig.
 */
public final class IntervalsSkipListOneContig<T extends Locatable> implements Serializable {
    private static final long serialVersionUID = 1L;

    // approx number of buckets we're aiming for.
    private static final int NUMBUCKETS = 1000;
    // each bucket contains 2**shift entries.
    private final int shift;

    // input intervals, sorted by start location
    private final List<T> vs;
    // the contig all the intervals are in.
    private final String contig;

    // reach: bucket# -> how far that bucket reaches.
    // e.g. bucket 0 contains the first 2**shift locatables. reach[0] is the max over their .getEnd()
    //      reach[x] is the max over the .getEnd for that bucket and all the ones before it.
    private final int[] reach;
    private final int reachLength;

    /**
     * Creates an IntervalsSkipList that holds a copy of the given intervals, sorted
     * and indexed.
     *
     * @param loc Locatables, not necessarily sorted. Will be iterated over exactly once.
     */
    public IntervalsSkipListOneContig(final Iterable<T> loc) {
        Utils.nonNull(loc);
        vs = Lists.newArrayList(loc);

        final Set<String> contigs = vs.stream().map(l -> l.getContig()).collect(Collectors.toSet());
        if (contigs.size() > 1){
            throw new IllegalArgumentException("Only one contig expected but got " + contigs);
        }

        if (vs.isEmpty()) {
            contig="";
        } else {
            contig=vs.get(0).getContig();
        }
        int bSize = vs.size() / NUMBUCKETS;
        // heuristic: if we have too many small buckets then we're better off instead
        // taking fewer but bigger steps, and then iterating through a few values.
        // Thus, put a lower bound on bucket size.
        if (bSize < 32) {
            bSize = 32;
        }
        shift = floorLog2(bSize);

        vs.sort(Comparator.comparing(Locatable::getContig).thenComparingInt(Locatable::getStart).thenComparing(Locatable::getEnd));

        reach = buildIndexAndCheck();
        reachLength = reach.length;
    }

    /**
     * Returns all the intervals that overlap with the query.
     * The query doesn't *have* to be in the same contig as the intervals we
     * hold, but of course if it isn't you'll get an empty result.
     * You may modify the returned list.
     */
    public List<T> getOverlapping(final SimpleInterval query) {
        Utils.nonNull(query);
        if (!contig.equals(query.getContig())) {
            // different contig, so we know no one'll overlap.
            return new ArrayList<>();
        }
        final List<T> ret = new ArrayList<>();
        // use index to skip early non-overlapping entries.
        int idx = firstPotentiallyReaching(query.getStart());
        if (idx<0) {
            idx=0;
        }
        for (;idx<vs.size();idx++) {
            final T v = vs.get(idx);
            // they are sorted by start location, so if this one starts too late
            // then all of the others will, too.
            if (v.getStart() > query.getEnd()) {
                break;
            }
            if (query.overlaps(v)) {
                ret.add(v);
            }
        }
        return ret;
    }

    // returns all the intervals that overlap with the query.
    // (use the optimized version instead, unless you're testing it and need something to compare against)
    @VisibleForTesting
    List<T> getOverlappingIgnoringIndex(final SimpleInterval query) {
        if (!contig.equals(query.getContig())) {
          // different contig, so we know no one'll overlap.
          return new ArrayList<>();
        }
        final List<T> ret = new ArrayList<>();
        for (final T v : vs) {
            // they are sorted by start location, so if this one starts too late
            // then all of the others will, too.
            if (v.getStart() > query.getEnd()) {
                break;
            }
            if (query.overlaps(v)) {
                ret.add(v);
            }
        }
        return ret;
    }

    // returns an index into the vs array s.t. no entry before that index
    // reaches (or extends beyond) the given position.
    private int firstPotentiallyReaching(final int position) {
        for (int i=0; i<reachLength; i++) {
            if (reach[i]>=position) {
                return i<<shift;
            }
        }
        // no one reaches to the given position.
        return vs.size()-1;
    }

    // build index and check everyone's in the same contig
    private int[] buildIndexAndCheck() {
        int max = 0;
        int key = 0;
        int idx = 0;
        // result: bucket# -> how far that bucket reaches
        final int[] result = new int[(vs.size()>>shift)+1];
        for (Locatable v : vs) {
            int k = idx>>shift;
            if (k>key) {
                result[key]=max;
                key=k;
            }
            if (v.getEnd()>max) {
                max=v.getEnd();
            }
            idx++;
        }
        result[key]=max;
        return result;
    }

    private static int floorLog2(final int n){
        if (n <= 0) {
            throw new IllegalArgumentException();
        }
        // size of int is 32 bits
        return 31 - Integer.numberOfLeadingZeros(n);
    }

}
