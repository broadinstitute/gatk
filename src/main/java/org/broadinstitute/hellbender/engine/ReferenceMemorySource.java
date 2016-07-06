package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.iterators.ByteArrayIterator;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;

/**
 * Manages traversals and queries over in-memory reference data.
 *
 * Supports targeted queries over the reference by interval, but does not
 * yet support complete iteration over the entire reference.
 */
public final class ReferenceMemorySource implements ReferenceDataSource {

    private final ReferenceBases bases;
    private final SAMSequenceDictionary sequenceDictionary;

    /**
     * Initialize this data source using ReferenceBases and corresponding sequence dictionary.
     */
    public ReferenceMemorySource(final ReferenceBases bases, final SAMSequenceDictionary referenceSequenceDictionary) {
        this.bases = Utils.nonNull(bases);
        this.sequenceDictionary = referenceSequenceDictionary;
    }

    /**
     * Start an iteration over the entire reference. Not yet supported!
     *
     * See the BaseUtils class for guidance on how to work with bases in this format.
     *
     * @return iterator over all bases in this reference
     */
    @Override
    public Iterator<Byte> iterator() {
        throw new UnsupportedOperationException("Iteration over entire reference not yet implemented");
    }

    /**
     * Query a specific interval on this reference, and get back an iterator over the bases spanning that interval.
     *
     * See the BaseUtils class for guidance on how to work with bases in this format.
     *
     * @param interval query interval
     * @return iterator over the bases spanning the query interval
     */
    @Override
    public Iterator<Byte> query( final SimpleInterval interval ) {

        int startIndex = (interval.getStart() - bases.getInterval().getStart());
        int stopIndex = startIndex + interval.size();

        return new ByteArrayIterator(bases.getBases(), startIndex, stopIndex);
    }

    /**
     * Query a specific interval on this reference, and get back all bases spanning that interval at once.
     * Call getBases() on the returned ReferenceSequence to get the actual reference bases. See the BaseUtils
     * class for guidance on how to work with bases in this format.
     *
     * @param interval query interval
     * @return a ReferenceSequence containing all bases spanning the query interval, prefetched
     */
    @Override
    public ReferenceSequence queryAndPrefetch( final SimpleInterval interval ) {
        return queryAndPrefetch(interval.getContig(), interval.getStart(), interval.getEnd());
    }

    /**
     * Query a specific interval on this reference, and get back all bases spanning that interval at once.
     * Call getBases() on the returned ReferenceSequence to get the actual reference bases. See the BaseUtils
     * class for guidance on how to work with bases in this format.
     *
     * @param contig query interval contig
     * @param start query interval start
     * @param stop query interval stop (included)
     * @return a ReferenceSequence containing all bases spanning the query interval, prefetched
     */
    @Override
    public ReferenceSequence queryAndPrefetch( final String contig, final long start , final long stop) {
        final int contigIndex = sequenceDictionary.getSequenceIndex(contig);
        int startIndex = (int)(start - bases.getInterval().getStart());
        int length = (int)(stop - start + 1);
        byte[] basesBytes = bases.getBases();
        if (startIndex==0 && length==basesBytes.length) {
            // special case: no need to make a copy
            return new ReferenceSequence(contig, contigIndex, basesBytes);
        }
        Utils.validIndex(startIndex, basesBytes.length);
        Utils.validateArg(startIndex+length <= basesBytes.length, () -> String.format("Asking for stop %d on contig %s but the ReferenceData only has data until %d.", stop, contig, bases.getInterval().getEnd()));
        Utils.validateArg(length >= 0, () -> String.format("Asking for stop<start (%d < %d)", stop, start));
        return new ReferenceSequence(contig, contigIndex, Arrays.copyOfRange(basesBytes, startIndex, startIndex+length));
    }


    /**
     * Get the sequence dictionary for this reference
     *
     * @return SAMSequenceDictionary for this reference
     */
    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return sequenceDictionary;
    }

    /**
     * no-op (nothing's open)
     */
    @Override
    public void close() {}
}
