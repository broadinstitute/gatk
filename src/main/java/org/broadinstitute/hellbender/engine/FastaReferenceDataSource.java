package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.iterators.ByteArrayIterator;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * Fasta-backed implementation of the reference data source.
 *
 * NOTE: does not yet support complete iteration over the entire reference.
 */
public final class FastaReferenceDataSource implements ReferenceDataSource {

    /**
     * Our reference file. Uses the caching version of IndexedFastaSequenceFile
     * so that repeated queries over nearby locations will be efficient (this
     * is the primary reference access pattern in most traversals).
     */
    private final CachingIndexedFastaSequenceFile reference;

    /**
     * Initialize this data source using a fasta file.
     *
     * The provided fasta file must have companion .fai and .dict files.
     *
     * @param fastaFile reference fasta file
     */
    public FastaReferenceDataSource( final File fastaFile ) {
        Utils.nonNull(fastaFile, "fastaFile must be non-null");

        // Will throw a UserException if the .fai and/or .dict are missing
        reference = CachingIndexedFastaSequenceFile.checkAndCreate(fastaFile);
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
    @Override
    public ReferenceSequence queryAndPrefetch( final String contig, final long start , final long stop) {
        return reference.getSubsequenceAt(contig, start, stop);
    }

    /**
     * Get the sequence dictionary for this reference
     *
     * @return SAMSequenceDictionary for this reference
     */
    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return reference.getSequenceDictionary();
    }

    /**
     * Permanently close this data source
     */
    @Override
    public void close() {
        try {
            reference.close();
        }
        catch ( final IOException e ) {
            throw new GATKException("Error closing reference file", e);
        }
    }
}
