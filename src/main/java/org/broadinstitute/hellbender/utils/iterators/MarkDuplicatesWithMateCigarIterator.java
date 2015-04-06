package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.markduplicates.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

/**
 * This will iterate through a coordinate sorted SAM file (iterator) and either mark or
 * remove duplicates as appropriate.  This class relies on the coordinate sort order as
 * well as the mate cigar (MC) optional SAM tag.
 */
public class MarkDuplicatesWithMateCigarIterator implements SAMRecordIterator {

    private SAMFileHeader header = null;

    /** The iterator from which records are read. */
    private PeekableIterator<SAMRecord> backingIterator = null;

    /** The ordinal of the next record to be read from the backing iterator */
    private int backingIteratorRecordIndex = 0;

    private boolean removeDuplicates = false;

    /** Should we skip pairs with no mate cigars or should be throw an error? */
    private boolean skipPairsWithNoMateCigar = true;
    private int numRecordsWithNoMateCigar = 0;

    /** When we hit unmapped reads that are just before the EOF, we can greedily process them as they will not have coordinates */
    private boolean foundUnmappedEOFReads = false;

    /** We can flush our queues and buffers if we move to a different reference index */
    private int referenceIndex = 0;

    /**
     * This buffer contains all the records read from input in the same order.  Nonetheless, each record
     * must be examined for duplicate marking, and so we may need to wait for this process to occur.  This
     * buffer stores the records in coordinate order, whether or not they can be emitted, and their associated
     * duplicate marking flag.  By definition, any record in the toMarkQueue will also be in the outputBuffer,
     * so we can omit checking the size of the toMarkQueue in some cases.
     */
    private SamRecordTrackingBuffer<SamRecordWithOrdinalAndSetDuplicateReadFlag> outputBuffer = null;

    /**
     * The queue that stores the records that currently are not marked as duplicates.  These need to be kept until
     * they cannot proven not to be duplicates, with the latter records having greater coordinate.  The queue is stored in 5' unclipped
     * ordering, along with keeping the record with the best score, defined by the scoring strategies.  If any record
     * is added to this queue and can be identified as a duplicate, the outputBuffer is notified of its
     * status and it can be emitted.  Therefore, we limit the amount of records in this queue to only those that will NOT
     * be duplicates.
     */
    private final MarkQueue toMarkQueue;

    /** The next record to be returned by next * */
    private SAMRecord nextRecord = null;

    /** This gets various information about the library id for a given record */
    private final LibraryIdGenerator libraryIdGenerator;

    /** This is used to identify optical duplicates among sets of duplicates */
    private OpticalDuplicateFinder opticalDuplicateFinder = null;

    /** We use this to check that the input data was in coordinate sort order */
    private final SAMRecordCoordinateComparator sortComparator = new SAMRecordCoordinateComparator();

    boolean isClosed = false;

    /**
     * Initializes the mark duplicates iterator.
     *
     * @param header                     the SAM header
     * @param iterator                   an iterator over the SAM records to consider
     * @param opticalDuplicateFinder     the algorithm for optical duplicate detection
     * @param duplicateScoringStrategy   the scoring strategy for choosing duplicates.  This cannot be SUM_OF_BASE_QUALITIES.
     * @param toMarkQueueMinimumDistance minimum distance for which to buffer
     * @param removeDuplicates           true to remove duplicates, false to mark duplicates
     * @param skipPairsWithNoMateCigar   true to not return mapped pairs with no mate cigar, false otherwise
     * @param blockSize                  the size of the blocks in the underlying buffer/queue
     * @param tmpDirs                    the temporary directories to use if we spill records to disk
     * @throws UserException if the inputs are not in coordinate sort order
     */
    public MarkDuplicatesWithMateCigarIterator(final SAMFileHeader header,
                                               final CloseableIterator<SAMRecord> iterator,
                                               final OpticalDuplicateFinder opticalDuplicateFinder,
                                               final DuplicateScoringStrategy.ScoringStrategy duplicateScoringStrategy,
                                               final int toMarkQueueMinimumDistance,
                                               final boolean removeDuplicates,
                                               final boolean skipPairsWithNoMateCigar,
                                               final int maxRecordsInRam,
                                               final int blockSize,
                                               final List<File> tmpDirs) throws UserException {
        if (header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new UserException(getClass().getName() + " expects the input to be in coordinate sort order.");
        }

        this.header = header;
        backingIterator = new PeekableIterator<SAMRecord>(iterator);
        outputBuffer = new SamRecordTrackingBuffer<SamRecordWithOrdinalAndSetDuplicateReadFlag>(maxRecordsInRam, blockSize, tmpDirs, header, SamRecordWithOrdinalAndSetDuplicateReadFlag.class);

        this.removeDuplicates = removeDuplicates;
        this.skipPairsWithNoMateCigar = skipPairsWithNoMateCigar;
        this.opticalDuplicateFinder = opticalDuplicateFinder;
        toMarkQueue = new MarkQueue(duplicateScoringStrategy);
        libraryIdGenerator = new LibraryIdGenerator(header);

        // Check for supported scoring strategies
        if (duplicateScoringStrategy == DuplicateScoringStrategy.ScoringStrategy.SUM_OF_BASE_QUALITIES)
            throw new UserException("SUM_OF_BASE_QUALITIES not supported as this may cause inconsistencies across ends in a pair.  Please use a different scoring strategy.");

        // set up metrics
        for (final SAMReadGroupRecord readGroup : header.getReadGroups()) {
            final String library = readGroup.getLibrary();
            DuplicationMetrics metrics = libraryIdGenerator.getMetricsByLibrary(library);
            if (metrics == null) {
                metrics = new DuplicationMetrics();
                metrics.LIBRARY = library;
                libraryIdGenerator.addMetricsByLibrary(library, metrics);
            }
        }

        // This sets the window size we need to keep to guarantee we can mark duplicates correctly
        toMarkQueue.setToMarkQueueMinimumDistance(toMarkQueueMinimumDistance);

        // get the first samRecordWithOrdinal
        nextRecord = markDuplicatesAndGetTheNextAvailable(); // get one directly, or null
    }

    public void logMemoryStats(final Log log) {
        System.gc();
        final Runtime runtime = Runtime.getRuntime();
        log.info("freeMemory: " + runtime.freeMemory() +
                        "; totalMemory: " + runtime.totalMemory() +
                        "; maxMemory: " + runtime.maxMemory() +
                        "; output buffer size: " + outputBuffer.size() +
                        "; duplicate queue size: " + toMarkQueue.size()
        );
    }

    /**
     * Establishes that records returned by this iterator are expected to
     * be in the specified sort order.  If this method has been called,
     * then implementers must throw an IllegalStateException from tmpReadEnds()
     * when a samRecordWithOrdinal is read that violates the sort order.  This method
     * may be called multiple times over the course of an iteration,
     * changing the expected sort, if desired -- from the time it is called,
     * it validates whatever sort is set, or stops validating if it
     * is set to null or SAMFileHeader.SortOrder.unsorted.  If this method
     * is not called, then no validation of the iterated records is done.
     *
     * @param sortOrder The order in which records are expected to be returned
     * @return This SAMRecordIterator
     */
    @Override
    public SAMRecordIterator assertSorted(final SAMFileHeader.SortOrder sortOrder) {
        if (sortOrder != SAMFileHeader.SortOrder.coordinate) {
            throw new IllegalStateException("Cannot assort " + sortOrder + " when expecting coordinate sorted input");
        }
        return this;
    }

    @Override
    public boolean hasNext() {
        // fast succeed
        if (null != nextRecord) return true;

        // We would need to get another record, so check if we can either a record read from the input to the mark queue, or we have more that we should return.
        // There should be at no time records in the mark queue that are not tracked in the output buffer.
        return (backingIterator.hasNext() || !outputBuffer.isEmpty());
    }

    @Override
    public SAMRecord next() throws UserException {
        final SAMRecord toReturn = nextRecord; // save for return


        // This should always return an element
        // NB: it should be the case that nextRecord != null
        if (null == toReturn) {
            throw new NoSuchElementException();
        }

        // Get the next record, if possible
        // NB: it should be the case that (nextRecord != null), due to the (null == toReturn) above
        if (hasNext()) {
            nextRecord = markDuplicatesAndGetTheNextAvailable(); // get one more, if possible
        } else {
            nextRecord = null;
        }

        // Check for sorted order
        if (null != nextRecord && 0 < sortComparator.fileOrderCompare(toReturn, nextRecord)) {
            System.err.print("Previous record: " + toReturn.getSAMString());
            System.err.print("Current record:" + nextRecord.getSAMString());
            throw new UserException("Records were not found coordinate sort order");
        }

        return toReturn;
    }

    /**
     * Handles records that are paired with both ends mapped, but lacking a mate cigar.  This returns true if we
     * can ignore this record after calling this method (when reading input), false otherwise.
     */
    private boolean ignoreDueToMissingMateCigar(final SamRecordWithOrdinal samRecordWithOrdinal) {
        final SAMRecord record = samRecordWithOrdinal.getRecord();
        // ignore/except-on paired records with mapped mate and no mate cigar
        if (record.getReadPairedFlag() &&
                !record.getMateUnmappedFlag() && null == SAMUtils.getMateCigar(record)) { // paired with one end unmapped and no mate cigar

            // NB: we are not truly examining these records. Do we want to count them?
            if (!record.isSecondaryOrSupplementary()) {
                // update metrics
                final DuplicationMetrics metrics = getMetrics(record);
                if (record.getReadUnmappedFlag()) {
                    ++metrics.UNMAPPED_READS;
                } else if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                    ++metrics.UNPAIRED_READS_EXAMINED;
                } else {
                    ++metrics.READ_PAIRS_EXAMINED;
                }
            }

            if (skipPairsWithNoMateCigar) { // pseudo-silently ignores them
                // NB: need to addRecordToTheOutputBuffer/set-flag as chunking/flushing of the toMarkQueue may need to occur
                addRecordToTheOutputBuffer(samRecordWithOrdinal); // now samRecordWithOrdinal will be stored in outputBuffer for return
                backingIteratorRecordIndex++;
                outputBuffer.setResultState(samRecordWithOrdinal, false); // indicate the present wrapped samRecordWithOrdinal is available for return
                numRecordsWithNoMateCigar++;
                backingIterator.next(); // remove it, since we called backingIterator.peek()
                return true;
            } else {
                throw new UserException("Read " + record.getReadName() + " was mapped and had a mapped mate, but no mate cigar (\"MC\") tag.");
            }
        }
        return false;
    }

    /**
     * This handles unmapped records at the end of the file.  If this is the first time we have found them, then we
     * can empty the toMarkQueue and call markDuplicatesAndGetTheNextAvailable, otherwise we can just emit them.  The
     * duplication metrics will be updated.
     */
    private SAMRecord nextIfRecordIsUnmappedAtEOF(final SAMRecord record) {
        // when we find unmapped reads with -1 as their reference index, there should be no mapped reads later in the file.
        if (foundUnmappedEOFReads) { // previously found unmapped reads at the end of the file
            final SAMRecord unmappedRecord = backingIterator.next(); // since we called backingIterator.peek()

            if (!record.isSecondaryOrSupplementary()) {
                // update metrics
                final DuplicationMetrics metrics = getMetrics(record);
                ++metrics.UNMAPPED_READS;
            }

            // We should have no more in the queue
            if (!outputBuffer.isEmpty()) {
                throw new UserException("Encountered unmapped reads at the end of the file, but the alignment start buffer was not empty.");
            }
            return unmappedRecord; // unmapped end of file records can simply be emitted - no need to duplicate mark them
        } else {
            foundUnmappedEOFReads = true;
            // move past all mapped reads
            referenceIndex = header.getSequenceDictionary().getSequences().size();

            // do the final round of duplicate marking
            tryPollingTheToMarkQueue(true, null);

            // NB: we do not call next here since we will recurse and perhaps hit the flush, or re-enter the if with unmapped EOF reads
            return markDuplicatesAndGetTheNextAvailable(); // this should flush the buffer
        }
    }

    /**
     * Check that we are not incorrectly performing any duplicate marking, by having too few of the records.  This
     * can happen if the alignment start is increasing but 5' soft-clipping is increasing such that we miss reads with
     * the same 5' unclipped alignment start.  This is especially true for RNAseq.
     */
    private void checkForMinimumDistanceFailure(final ReadEndsForMateCigar current) {
        if (!toMarkQueue.isEmpty()) {
            final ReadEndsForMateCigar other = toMarkQueue.peek();
            if (other.read1ReferenceIndex == current.read1ReferenceIndex && toMarkQueue.getToMarkQueueMinimumDistance() <= other.read1Coordinate - current.read1Coordinate) {
                if (checkCigarForSkips(other.getRecord().getCigar())) {
                    throw new GATKException("Found a samRecordWithOrdinal with sufficiently large code length that we may have\n"
                            + " missed including it in an early duplicate marking iteration.  Alignment contains skipped"
                            + " reference bases (N's). If this is an\n RNAseq aligned bam, please use MarkDuplicates instead,"
                            + " as this tool does not work well with spliced reads.\n Minimum distance set to "
                            + toMarkQueue.getToMarkQueueMinimumDistance() + " but " + (other.read1Coordinate - current.read1Coordinate - 1)
                            + " would be required.\n" + "Record was: " + other.getRecord().getSAMString());
                } else {
                    System.err.print("record #1: " + other.getRecord().getSAMString());
                    System.err.print("record #2: " + current.getRecord().getSAMString());
                    throw new GATKException("Found a samRecordWithOrdinal with sufficiently large clipping that we may have\n"
                            + " missed including it in an early duplicate marking iteration.  Please increase the"
                            + " minimum distance to at least " + (other.read1Coordinate - current.read1Coordinate - 1)
                            + "bp\nto ensure it is considered (was " + toMarkQueue.getToMarkQueueMinimumDistance() + ").\n"
                            + "Record was: " + other.getRecord().getSAMString());
                }
            }
        }
    }

    /**
     * This tries to get a record that has been evaluated for duplicate marking.  It does this by first seeing if there
     * are any records that have been through duplicate marking.  If none are available, it will try to get more records
     * from the input iterator until there are reads available that have been duplicate marked.  If there are no more
     * records available from the input iterator, it will duplicate mark the final chunk of records.  Finally, if there
     * are no more records, it will return null;
     */
    private SAMRecord markDuplicatesAndGetTheNextAvailable() {

        // Check if there are any we can flush output buffer
        { // NB: braces to limit the scope of 'record'
            final SAMRecord record = flush();
            if (null != record) return record;
        }

        // Check if there are any more records to read in
        if (!backingIterator.hasNext()) { // no more records to read in

            // Check if there are any more to mark
            if (toMarkQueue.isEmpty()) {
                // check if there are any that can be outputted
                if (outputBuffer.isEmpty()) {
                    return null;
                } // no need to flush; no records in the queue and buffer
            } else {
                // force marking duplicates on the remaining records
                tryPollingTheToMarkQueue(true, null);
            }

            /** Since we have no more records to read in, and no more records that need duplicate marking run, we
             * update our coordinate to past the end of the reference
             */
            referenceIndex = header.getSequenceDictionary().getSequences().size();

            /** Now we recurse, so that we can flush from the outputBuffer until it is empty, then return null when
             * all of the input, queue, and buffer are empty */
            return markDuplicatesAndGetTheNextAvailable();
        }

        /** We need to retrieve more records from the input iterator and duplicate mark, until we can return one that
         *  has been through duplicate marking.
         */
        while (backingIterator.hasNext()) {

            // NB: we could get rid of this if we made nextRecord into a list...
            // NB: we do not actually remove this record from the backing iterator until much later, to help with processing unmapped reads at the EOF
            SAMRecord record = backingIterator.peek(); // peek: used for unmapped reads
            final SamRecordWithOrdinal samRecordWithOrdinal = new SamRecordWithOrdinalAndSetDuplicateReadFlag(record, backingIteratorRecordIndex);

            ReadEndsForMateCigar readEnds = null;
            boolean performedChunkAndMarkTheDuplicates = false;

            // remove duplicate information
            record.setDuplicateReadFlag(false);

            /** Check for pairs that have both ends mapped and missing mate cigar. */
            if (ignoreDueToMissingMateCigar(samRecordWithOrdinal)) {
                continue;
            }

            // check for an unmapped read
            if (record.getReadUnmappedFlag()) {
                // unmapped reads at the end of the file!
                if (-1 == record.getReferenceIndex()) {
                    // NB: this may call markDuplicatesAndGetTheNextAvailable if this is the first time a EOF unmapped record has been seen
                    return nextIfRecordIsUnmappedAtEOF(record);
                } else if (!record.isSecondaryOrSupplementary()) {
                    // update metrics
                    final DuplicationMetrics metrics = getMetrics(record);
                    ++metrics.UNMAPPED_READS;
                }
                // we will check for unmapped reads later so as not to add them to mark queue
            } else {
                // If not already set, this sets the minimum distance to twice the read length, or 100, whichever is larger
                if (-1 == toMarkQueue.getToMarkQueueMinimumDistance()) {
                    // use twice the first read's length
                    toMarkQueue.setToMarkQueueMinimumDistance(Math.max(2 * record.getReadBases().length, 100));
                }

                // build a read end for use in the toMarkQueue
                readEnds = new ReadEndsForMateCigar(header, samRecordWithOrdinal, opticalDuplicateFinder, libraryIdGenerator.getLibraryId(samRecordWithOrdinal.getRecord()));

                // check that the minimumDistance was not too small
                checkForMinimumDistanceFailure(readEnds);

                /**
                 * If we can do some duplicate marking, lets do it!
                 * IMPORTANT: this does not flush the to-mark-queue, so the minimum distance needs to be set for us to infer
                 * which records will never be supplemented (i.e. are non-duplicate).
                 */
                performedChunkAndMarkTheDuplicates = tryPollingTheToMarkQueue(false, readEnds);
            }

            // We can now remove the record from the input
            backingIterator.next();

            // Add this to the outputBuffer so it can be tracked.  It will not be available to emit until it has been through duplicate marking.
            addRecordToTheOutputBuffer(samRecordWithOrdinal);
            backingIteratorRecordIndex++; // Each record is has an index and is emitted in the same order. This helps that.

            // We do not consider secondary, supplementary, or unmapped alignments for duplicate marking. We can thus mark that duplicate marking on them has been completed.
            if (record.isSecondaryOrSupplementary() || record.getReadUnmappedFlag()) {
                outputBuffer.setResultState(samRecordWithOrdinal, false);
            } else {
                // Bring the simple metrics up to date
                final DuplicationMetrics metrics = getMetrics(record);
                if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                    ++metrics.UNPAIRED_READS_EXAMINED;
                } else {
                    ++metrics.READ_PAIRS_EXAMINED; // will need to be divided by 2 at the end
                }

                // Add the record for duplicate marking, which may in fact cause it to be duplicate marked or stored for later
                toMarkQueue.add(readEnds, outputBuffer, getMetrics(readEnds.getRecord()));
            }

            // Check if there are any we can flush, which happens if we just performed duplicate marking
            if (performedChunkAndMarkTheDuplicates) {
                record = flush();
                if (null != record) return record;
            }
        }

        // try again, as we may have records we can flush, or we want to see if we are at the EOF
        return markDuplicatesAndGetTheNextAvailable();
    }


    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void close() {
        // close the input and output
        backingIterator.close();
        outputBuffer.close();
        isClosed = true;
    }

    /**
     * Checks a Cigar for the presence of N operators. Reads with skipped bases may be spliced RNAseq reads
     *
     * @param cigar
     */
    private boolean checkCigarForSkips(final Cigar cigar) {
        final List<CigarElement> elements = cigar.getCigarElements();
        for (final CigarElement el : elements) {
            if (el.getOperator() == CigarOperator.N) return true;
        }
        return false;
    }

    private void enforceClosed() {
        if (!isClosed) throw new GATKException("Calling a method that assumes the iterator is closed");
    }

    /** Useful for statistics after the iterator has been exhausted and closed. */
    public int getNumRecordsWithNoMateCigar() {
        enforceClosed();
        return numRecordsWithNoMateCigar;
    }

    public int getNumDuplicates() {
        enforceClosed();
        return toMarkQueue.getNumDuplicates();
    }

    public LibraryIdGenerator getLibraryIdGenerator() {
        enforceClosed();
        return libraryIdGenerator;
    }

    public Histogram<Short> getOpticalDupesByLibraryId() {
        enforceClosed();
        return libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap();
    }

    /**
     * Gets a SAMRecord if one is available after marking.  This enforces that we return records in the original
     * coordinate sort order in a stable fashion.
     *
     * @return record representing the head of the alignment-start sorted buffer, or null if the head record has not yet been duplicate marked
     */
    private SAMRecord flush() {
        // Check that there is at least one record in the coordinate-sorted buffer, and that the head record has been through duplicate-marking
        while (!outputBuffer.isEmpty() && outputBuffer.canEmit()) {
            // the buffer contains wrapped SAMRecords, which we want to unwrap
            final SAMRecord record = outputBuffer.next().getRecord();

            // If this read is a duplicate, do we want to remove it (continue the loop) or return it for emission?
            if (!removeDuplicates || !record.getDuplicateReadFlag()) {
                return record;
            }
        }
        return null;
    }

    /**
     * Adds a samRecordWithOrdinal to the output buffer.  This does not mean that it is ready to be emitted, since it may need to be
     * duplicate marked.
     *
     * @param samRecordWithOrdinal the index of the record of which to track.
     * @throws GATKException if the records are added out of order
     */
    private void addRecordToTheOutputBuffer(final SamRecordWithOrdinal samRecordWithOrdinal) throws GATKException {
        final int recordReferenceIndex = samRecordWithOrdinal.getRecord().getReferenceIndex();
        if (recordReferenceIndex < referenceIndex) {
            throw new GATKException("Records out of order: " + recordReferenceIndex + " < " + referenceIndex);
        } else if (referenceIndex < recordReferenceIndex) {
            // new reference, so we need to mark duplicates on the current ones
            // NB: we will not miss inter-chromosomal alignments since presumably one end will have been mapped to this chromosome and processed, and we do not need the other end to do so.
            tryPollingTheToMarkQueue(true, null);
            // update genomic coordinate to the next reference index
            referenceIndex = recordReferenceIndex;
        }

        // add the samRecordWithOrdinal to the output buffer so that it can be tracked
        outputBuffer.add(samRecordWithOrdinal);
    }

    /**
     * Tries to get a record from the toMarkQueue that has been successfully through duplicate marking.  Note, either flush is true or
     * current must be non-null.
     *
     * @param flush   true if we should empty the toMarkQueue fully.
     * @param current the current end to ensure we consider all possible ends for a duplicate
     * @return true if we did get at least one record, false otherwise
     */
    private boolean tryPollingTheToMarkQueue(final boolean flush, final ReadEndsForMateCigar current) {
        boolean performedChunkAndMarkTheDuplicates = false;

        if (!flush && null == current) throw new GATKException("Flush cannot be false and current be null");

        if (toMarkQueue.isEmpty()) return false;

        if (!toMarkQueue.isEmpty() && outputBuffer.isEmpty()) {
            throw new GATKException("0 < toMarkQueue && outputBuffer.isEmpty()");
        }

        /**
         * Try to poll the toMarkQueue.  If we are flushing all the records from it, just do so until empty.  Otherwise, we need to
         * make sure we only poll those a certain distance away from current.
         */
        while (!toMarkQueue.isEmpty() &&
                (flush || referenceIndex != current.read1ReferenceIndex ||
                        toMarkQueue.getToMarkQueueMinimumDistance() < current.read1Coordinate - toMarkQueue.peek().read1Coordinate)) {

            // Poll will track that this samRecordWithOrdinal has been through duplicate marking. It is not marked as a duplicate :)
            final ReadEndsForMateCigar next = toMarkQueue.poll(outputBuffer, header, opticalDuplicateFinder, libraryIdGenerator); // get the first one!
            performedChunkAndMarkTheDuplicates = true;

            // track optical duplicates using only those reads that are the first end...
            if (toMarkQueue.shouldBeInLocations(next) && next.getRecord().getFirstOfPairFlag()) {
                final Set<ReadEnds> locations = toMarkQueue.getLocations(next);

                if (!locations.isEmpty()) {
                    AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(new ArrayList<ReadEnds>(locations),
                            opticalDuplicateFinder, libraryIdGenerator);
                }
            }
            // NB: we could try to greedily return a record if one is available here.  Instead we continue processing the mark queue */
        }
        return performedChunkAndMarkTheDuplicates;
    }

    /** Get the duplication metrics for the library associated with end. */
    private DuplicationMetrics getMetrics(final SAMRecord record) {
        final String library = LibraryIdGenerator.getLibraryName(header, record);
        DuplicationMetrics metrics = libraryIdGenerator.getMetricsByLibrary(library);
        if (metrics == null) {
            metrics = new DuplicationMetrics();
            metrics.LIBRARY = library;
            libraryIdGenerator.addMetricsByLibrary(library, metrics);
        }
        return metrics;
    }
}
