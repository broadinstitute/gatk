package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;


import java.util.*;

/**
 * Returns a SimpleInterval for each locus in a set of intervals.  I.e. returns each genomic point location in an interval list.
 */
public class IntervalLocusIterator implements Iterator<SimpleInterval> {

    private Iterator<SimpleInterval> intervalIterator;

    public SimpleInterval currentInterval = null;

    /** Encoded as an interval that i of size 1 in the current interval */
    private Iterator<SimpleInterval> baseLocationIterator;

    public IntervalLocusIterator(final Iterator<SimpleInterval> intervalIterator) {
        Utils.nonNull(intervalIterator, "Input iterator cannot be null");
        this.intervalIterator = intervalIterator;
        advanceCurrentInterval();
    }

    @Override
    public boolean hasNext() {
        if (currentInterval == null) {
            return false;
        }
        if (!baseLocationIterator.hasNext()) {
            if (!intervalIterator.hasNext()) {
                return false;
            }
        }
        return true;
    }

    @Override
    public SimpleInterval next() {
        if (baseLocationIterator.hasNext()) {
            return baseLocationIterator.next();
        }
        if (intervalIterator.hasNext()) {
            advanceCurrentInterval();
            return baseLocationIterator.next();
        }
        throw new NoSuchElementException();
    }

    private void advanceCurrentInterval() {
        if (intervalIterator.hasNext()) {
            currentInterval = intervalIterator.next();
        } else {
            // Typically, this code block should only get hit when this class is
            //  initalized with an empty interval list.
            currentInterval = null;
        }
        baseLocationIterator = createBaseLocationIterator(currentInterval);
    }


    private Iterator<SimpleInterval> createBaseLocationIterator(final SimpleInterval fullInterval) {
        if (fullInterval == null) {
            return Collections.<SimpleInterval>emptyList().iterator();
        }

        List<SimpleInterval> dummyList = Collections.singletonList(fullInterval);
        return new ShardedIntervalIterator(dummyList.iterator(), 1);
    }
}
