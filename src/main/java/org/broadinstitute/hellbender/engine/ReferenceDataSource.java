package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.iterators.ByteArrayIterator;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * Manages traversals and queries over reference data (for now, fasta files only)
 *
 * Supports targeted queries over the reference by interval, but does not
 * yet support complete iteration over the entire reference.
 */
public class ReferenceDataSource implements GATKDataSource<Byte>, AutoCloseable {

    /**
     * Our reference file. Uses the caching version of IndexedFastaSequenceFile
     * so that repeated queries over nearby locations will be efficient (this
     * is the primary reference access pattern in most traversals).
     */
    private CachingIndexedFastaSequenceFile reference;

    /**
     * Initialize this data source using a fasta file.
     *
     * The provided fasta file must have companion .fai and .dict files.
     *
     * @param fastaFile reference fasta file
     */
    public ReferenceDataSource( final File fastaFile ) {
        if ( fastaFile == null )
            throw new IllegalArgumentException("fastaFile must be non-null");

        // Will throw a UserException if the .fai and/or .dict are missing
        reference = CachingIndexedFastaSequenceFile.checkAndCreate(fastaFile);
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
    public Iterator<Byte> query( final GenomeLoc interval ) {
        // TODO: need a way to iterate lazily over reference bases without necessarily loading them all into memory at once
        return new ByteArrayIterator(queryAndPrefetch(interval).getBases());
    }

    /**
     * Query a specific interval on this reference, and get back all bases spanning that interval at once.
     * Call getBases() on the returned ReferenceSequence to get the actual reference bases. See the BaseUtils
     * class for guidance on how to work with bases in this format.
     *
     * @param interval query interval
     * @return a ReferenceSequence containing all bases spanning the query interval, prefetched
     */
    public ReferenceSequence queryAndPrefetch( final GenomeLoc interval ) {
        return queryAndPrefetch(interval.getContig(), interval.getStart(), interval.getStop());
    }

    /**
     * Query a specific interval on this reference, and get back all bases spanning that interval at once.
     * Call getBases() on the returned ReferenceSequence to get the actual reference bases. See the BaseUtils
     * class for guidance on how to work with bases in this format.
     *
     * @param contig query interval contig
     * @param start query interval start
     * @param stop query interval stop
     * @return a ReferenceSequence containing all bases spanning the query interval, prefetched
     */
    public ReferenceSequence queryAndPrefetch( final String contig, final long start , final long stop) {
        return reference.getSubsequenceAt(contig, start, stop);
    }


    /**
     * Get the sequence dictionary for this reference
     *
     * @return SAMSequenceDictionary for this reference
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        return reference.getSequenceDictionary();
    }

    /**
     * Permanently close this data source
     */
    public void close() {
        try {
            reference.close();
        }
        catch ( IOException e ) {
            throw new GATKException("Error closing reference file", e);
        }
    }
}
