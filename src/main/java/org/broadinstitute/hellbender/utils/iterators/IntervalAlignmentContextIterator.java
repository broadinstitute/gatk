package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * For special cases where we want to emit AlignmentContexts regardless of whether we have an overlap with a given interval.
 *
 * In other words, this will iterate over each base of each interval and emit an alignment context.  If no reads overlap,
 *  it will emit an empty alignment context.  Each empty alignment context will be a new instance.
 */
public class IntervalAlignmentContextIterator implements Iterable<AlignmentContext>, Iterator<AlignmentContext> {
    private Iterator<AlignmentContext> alignmentContextIterator;
    private IntervalLocusIterator intervalLocusIterator;
    private SimpleInterval currentInterval;
    private AlignmentContext currentAlignmentContext;
    private SAMSequenceDictionary dictionary;

    public IntervalAlignmentContextIterator(Iterator<AlignmentContext> alignmentContextIterator, IntervalLocusIterator intervalLocusIterator, SAMSequenceDictionary dictionary) {
        this.alignmentContextIterator = alignmentContextIterator;
        this.intervalLocusIterator = intervalLocusIterator;
        this.dictionary = dictionary;

        // Iterate once.
        advanceIntervalLocus();
        advanceAlignmentContext();
        advanceAlignmentContextToCurrentInterval(this.currentInterval);
    }

    @Override
    public Iterator<AlignmentContext> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        return (intervalLocusIterator.hasNext()) || (currentInterval != null);
    }

    @Override
    public AlignmentContext next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }

        final boolean isOverlaps = currentInterval.overlaps(currentAlignmentContext);
        AlignmentContext result;

        // Get ready to return result and determine what should go the next time next() is called.
        if (isOverlaps) {
            result = currentAlignmentContext;
            advanceIntervalLocus();
            advanceAlignmentContextToCurrentInterval(currentInterval);
        } else {
            result = createEmptyAlignmentContext(currentInterval);
            advanceIntervalLocus();

            if (currentInterval != null) {
                final int comparison = IntervalUtils.compareLocatables(currentInterval, currentAlignmentContext, dictionary);
                // Comparison should not be zero and if less than zero, we'd just advance the intervalLocus iterator,
                //  which we are already doing.
                //  Interval is after the next alignment context when comparison is greater than zero.
                if (comparison > 0) {
                    advanceAlignmentContextToCurrentInterval(currentInterval);
                }
            }
        }

        return result;

    }

    private AlignmentContext createEmptyAlignmentContext(final SimpleInterval interval) {
        return new AlignmentContext(interval, new ReadPileup(interval, new ArrayList<>()));
    }

    private void advanceAlignmentContextToCurrentInterval(final SimpleInterval interval) {
        if (interval == null) {
            currentAlignmentContext = null;
            return;
        }
        while (IntervalUtils.compareLocatables(currentInterval, currentAlignmentContext, dictionary) > 0) {
            advanceAlignmentContext();
        }
    }

    private void advanceAlignmentContext() {
        if (alignmentContextIterator.hasNext()) {
            currentAlignmentContext = this.alignmentContextIterator.next();
        } else {
            currentAlignmentContext = createEmptyAlignmentContext(currentInterval);
        }
    }

    private void advanceIntervalLocus() {
        if (intervalLocusIterator.hasNext()) {
            currentInterval = intervalLocusIterator.next();
        } else {
            currentInterval = null;
        }
    }
}
