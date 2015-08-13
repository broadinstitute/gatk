package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.iterators.ByteArrayIterator;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * Manages traversals and queries over reference data (for now, fasta files only)
 *
 * Supports targeted queries over the reference by interval, but does not
 * yet support complete iteration over the entire reference.
 */
public interface ReferenceDataSource extends GATKDataSource<Byte>, AutoCloseable {

    /**
     * Initialize this data source using a fasta file.
     *
     * The provided fasta file must have companion .fai and .dict files.
     *
     * @param fastaFile reference fasta file
     */
    public static ReferenceDataSourceFromFile of(final File fastaFile) {
        return new ReferenceDataSourceFromFile(fastaFile);
    }

    /**
     * Initialize this data source using ReferenceBases and corresponding sequence dictionary.
     */
    public static ReferenceDataSourceFromReferenceBases of(final ReferenceBases bases, final SAMSequenceDictionary referenceSequenceDictionary) {
        return new ReferenceDataSourceFromReferenceBases(bases, referenceSequenceDictionary);
    }


    /**
     * Query a specific interval on this reference, and get back all bases spanning that interval at once.
     * Call getBases() on the returned ReferenceSequence to get the actual reference bases. See the BaseUtils
     * class for guidance on how to work with bases in this format.
     *
     * @param interval query interval
     * @return a ReferenceSequence containing all bases spanning the query interval, prefetched
     */
    public ReferenceSequence queryAndPrefetch( final SimpleInterval interval );

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
    public ReferenceSequence queryAndPrefetch( final String contig, final long start , final long stop);


    /**
     * Get the sequence dictionary for this reference
     *
     * @return SAMSequenceDictionary for this reference
     */
    public SAMSequenceDictionary getSequenceDictionary();

    @Override
    public void close();

}
