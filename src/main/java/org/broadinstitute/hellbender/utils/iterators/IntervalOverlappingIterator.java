package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Wraps an iterator of {@link htsjdk.samtools.util.Locatable} with a list of sorted intervals
 * to return only the objects which overlaps with them
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class IntervalOverlappingIterator<T extends Locatable> implements Iterator<T> {

    // underlying iterator
    private final Iterator<T> iterator;

    // sorted intervals
    private final Iterator<SimpleInterval> intervals;

    // the current interval to check
    private SimpleInterval currentInterval;

    // the sequence dictionary to get contig ordering
    private final SAMSequenceDictionary dictionary;

    // the next object to return
    private T next;

    /**
     * Wraps an iterator to be filtered by a sorted list of intervals
     *
     * @param iterator   underlying iterator
     * @param intervals  sorted list of intervals to traverse
     * @param dictionary sequence dictionary indicating the ordering of contigs
     */
    public IntervalOverlappingIterator(Iterator<T> iterator, List<SimpleInterval> intervals,
        SAMSequenceDictionary dictionary) {
        Utils.nonNull(iterator);
        Utils.nonEmpty(intervals);
        Utils.nonNull(dictionary);
        this.iterator = iterator;
        this.intervals = intervals.iterator();
        this.dictionary = dictionary;
        currentInterval = this.intervals.next();
        advance();
    }

    @Override
    public boolean hasNext() {
        return next != null;
    }

    @Override
    public T next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        T toReturn = next;
        advance();
        return toReturn;
    }

    /**
     * Advance to the next location, setting next to null if there are no more records
     */
    private void advance() {
        // get the next record to check
        next = (iterator.hasNext()) ? iterator.next() : null;
        // iterate till next or current interval is null
        while(next != null && currentInterval != null) {
            if(currentInterval.overlaps(next)) {
                // keep the next because it overlaps
                return;
            } else {
                int comparison = IntervalUtils.compareLocatables(currentInterval, next, dictionary);
                // only advance if the current interval is before the next record
                if(comparison < 0) {
                    // advance the interval and try with th next one
                    currentInterval = (intervals.hasNext()) ? intervals.next() : null;
                } else if(comparison > 0)  {
                    // advance the location
                    next = (iterator.hasNext()) ? iterator.next() : null;
                }
            }
        }
        // if the value of next overlaps some interval, the method should return before this point
        next = null;
    }
}
