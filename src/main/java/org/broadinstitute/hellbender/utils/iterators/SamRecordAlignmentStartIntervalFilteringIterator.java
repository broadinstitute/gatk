package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Iterator that will filter the {@link SAMRecord}s from the given {@link CloseableIterator<SAMRecord>} by whether they appear
 * in the given {@link Iterator<SimpleInterval>}.
 * If a {@link SAMRecord} starts in any of the intervals given, then the record passes filtering and is emitted.
 *
 * NOTE:
 *      Assumes that both the given {@link CloseableIterator<SAMRecord>} and {@link Iterator<SimpleInterval>} contain
 *      sorted data.
 * Created by jonn on 7/23/19.
 */
public class SamRecordAlignmentStartIntervalFilteringIterator implements CloseableIterator<SAMRecord> {

    final private CloseableIterator<SAMRecord> samRecordIterator;
    final private Iterator<SimpleInterval> intervalIterator;

    private SAMSequenceDictionary sequenceDictionary;

    private SimpleInterval currentInterval = null;

    private SAMRecord next = null;

    public SamRecordAlignmentStartIntervalFilteringIterator(
            final SAMSequenceDictionary sequenceDictionary,
            final Iterator<SimpleInterval> intervalIterator,
            final CloseableIterator<SAMRecord> samRecordIterator) {

        Utils.nonNull(sequenceDictionary, "Sequence dictionary cannot be null");
        Utils.nonNull(intervalIterator, "Input interval iterator cannot be null");
        Utils.nonNull(samRecordIterator, "Input sam record iterator cannot be null");
        Utils.validate(intervalIterator.hasNext(), "Given intervals iterator should never be empty.");

        this.sequenceDictionary = sequenceDictionary;
        this.intervalIterator = intervalIterator;
        this.samRecordIterator = samRecordIterator;

        if (samRecordIterator.hasNext()) {
            advanceCurrentInterval();
            advanceIterators();
        }
    }

    @Override
    public void close() {
        if ( samRecordIterator != null ) {
            samRecordIterator.close();
        }
    }

    @Override
    public boolean hasNext() {
        return next != null;
    }

    @Override
    public SAMRecord next() {

        if ( !hasNext() ) {
            throw new NoSuchElementException();
        }
        final SAMRecord toReturn = next;
        advanceIterators();
        return toReturn;
    }

    private void advanceIterators() {
        // Only continue if we have more data:
        if ( !samRecordIterator.hasNext() ) {
            next = null;
            return;
        }

        // Get our record and the start position.
        // These will be updated in the while loop below if they do not overlap our intervals:
        SAMRecord samRecord = samRecordIterator.next();
        SimpleInterval recordStart =
                new SimpleInterval( samRecord.getContig(), samRecord.getStart(), samRecord.getStart() );

        // Since we are sorted we can do the following to keep it fast:

        // While we are not in an interval that contains the current samRecord, advance the sam records.
        // If the sam records jump contigs, we jump contigs and try again.
        while ( (currentInterval != null) && (!currentInterval.overlaps( recordStart )) ) {

            // If our interval no longer can align with our record, we increment the interval:
            if ( mustAdvanceIntervals(samRecord) ) {
                advanceCurrentInterval();
            }
            else {
                if ( !samRecordIterator.hasNext() ) {
                    next = null;
                    return;
                }
                samRecord = samRecordIterator.next();
                recordStart = new SimpleInterval( samRecord.getContig(), samRecord.getStart(), samRecord.getStart() );
            }
        }

        // We either have a containing interval or we're done:
        if ( currentInterval == null ) {
            next = null;
        }
        else {
            next = samRecord;
        }
    }

    /**
     * Determines whether to advance the intervals iterator or the reads iterator.
     * @param samRecord The current {@link SAMRecord} in the reads iterator.
     * @return {@code true} if the interval iterator should be updated.  {@code false} if the reads iterator should be updated.
     */
    private boolean mustAdvanceIntervals(final SAMRecord samRecord) {

        if ( !currentInterval.getContig().equals( samRecord.getContig() )) {
            // Different contigs.
            // We must figure out whether the reads are behind or the intervals are behind:

            final int intervalContigIndex = sequenceDictionary.getSequenceIndex( currentInterval.getContig() );
            final int readContigIndex = sequenceDictionary.getSequenceIndex( samRecord.getContig() );

            if ( intervalContigIndex < readContigIndex ) {
                return true;
            }
            else {
                return false;
            }
        }

        if ( !currentInterval.overlaps(samRecord) ) {
            // The interval is on the same contig but does not overlap at all.
            // Because the intervals and reads are in order, this means we need to move to the next interval.
            return true;
        }

        return false;
    }

    private void advanceCurrentInterval() {
        if (intervalIterator.hasNext()) {
            currentInterval = intervalIterator.next();
        } else {
            // This code block should only get hit when there are no more intervals in the intervalIterator.
            currentInterval = null;
        }
    }
}
