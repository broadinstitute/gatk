package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentContextIteratorBuilder;
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
public class IntervalAlignmentContextIterator implements Iterator<AlignmentContext> {
    private Iterator<AlignmentContext> alignmentContextIterator;
    private IntervalLocusIterator intervalLocusIterator;
    private SimpleInterval currentInterval;
    private AlignmentContext currentAlignmentContext;
    private SAMSequenceDictionary dictionary;


    /**
     *  Note:  Typically, if you are calling this from a walker tool, you want to use {@link AlignmentContextIteratorBuilder}
     *
     * @param alignmentContextIterator iterator for alignment context
     * @param intervalLocusIterator an iterator that will traverse all relevant loci
     * @param dictionary reference or best available dictionary.
     */
    public IntervalAlignmentContextIterator(Iterator<AlignmentContext> alignmentContextIterator, IntervalLocusIterator intervalLocusIterator, SAMSequenceDictionary dictionary) {
        this.alignmentContextIterator = alignmentContextIterator;
        this.intervalLocusIterator = intervalLocusIterator;
        this.dictionary = dictionary;

        // Iterate once.
        advanceIntervalLocus();

        // Just to make sure that current alignment context is not null.  We only want to do this during the initialization.
        advanceAlignmentContext();

        // Advance the alignment context to overlap or be in front of the current interval.
        advanceAlignmentContextToCurrentInterval();
    }

    @Override
    public boolean hasNext() {
        return currentInterval != null;
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
            advanceAlignmentContextToCurrentInterval();
        } else {
            result = createEmptyAlignmentContext(currentInterval);
            advanceIntervalLocus();

            if (currentInterval != null) {
                final int comparison = IntervalUtils.compareLocatables(currentInterval, currentAlignmentContext, dictionary);
                // Comparison should not be zero and if less than zero, we'd just advance the intervalLocus iterator,
                //  which we are already doing.
                //  Interval is after the next alignment context when comparison is greater than zero.
                if (comparison > 0) {
                    advanceAlignmentContextToCurrentInterval();
                }
            }
        }

        return result;

    }

    private AlignmentContext createEmptyAlignmentContext(final SimpleInterval interval) {
        return new AlignmentContext(interval, new ReadPileup(interval, new ArrayList<>()));
    }

    private void advanceAlignmentContextToCurrentInterval() {
        if (currentInterval == null) {
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
