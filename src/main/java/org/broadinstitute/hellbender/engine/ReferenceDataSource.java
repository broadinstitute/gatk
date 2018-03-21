package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.iterators.ByteArrayIterator;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.nio.file.Path;
import java.util.Iterator;

/**
 * Manages traversals and queries over reference data.
 *
 * Supports targeted queries over the reference by interval and over the entire reference.
 */
public interface ReferenceDataSource extends GATKDataSource<Byte>, AutoCloseable {

    /**
     * Initialize this data source using a fasta file.
     *
     * The provided fasta file must have companion .fai and .dict files.
     *
     * @param fastaPath reference fasta Path
     */
    public static ReferenceDataSource of(final Path fastaPath) {
        return new ReferenceFileSource(fastaPath);
    }

    /**
     * Initialize this data source using a fasta file.
     *
     * The provided fasta file must have companion .fai and .dict files.
     *
     * If {@code preserveFileBases} is {@code true}, will NOT convert IUPAC bases in the file to `N` and will NOT capitalize lower-case bases.
     *
     * NOTE: Most GATK tools do not support data created by setting {@code preserveFileBases} to {@code true}.
     *
     * @param fastaPath reference fasta Path
     * @param preserveAmbiguityCodesAndCapitalization Whether to preserve the original bases in the given reference file path.
     */
    public static ReferenceDataSource of(final Path fastaPath, final boolean preserveAmbiguityCodesAndCapitalization) {
        return new ReferenceFileSource(fastaPath, preserveAmbiguityCodesAndCapitalization);
    }

    /**
     * Initialize this data source using ReferenceBases and corresponding sequence dictionary.
     */
    public static ReferenceDataSource of(final ReferenceBases bases, final SAMSequenceDictionary referenceSequenceDictionary) {
        return new ReferenceMemorySource(bases, referenceSequenceDictionary);
    }

    /**
     * Query a specific interval on this reference, and get back all bases spanning that interval at once.
     * Call getBases() on the returned ReferenceSequence to get the actual reference bases. See the BaseUtils
     * class for guidance on how to work with bases in this format.
     *
     * The default implementation calls #queryAndPrefetch(contig, start, stop).
     *
     * @param interval query interval
     * @return a ReferenceSequence containing all bases spanning the query interval, prefetched
     */
    default public ReferenceSequence queryAndPrefetch( final SimpleInterval interval ) {
        return queryAndPrefetch(interval.getContig(), interval.getStart(), interval.getEnd());
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
    public ReferenceSequence queryAndPrefetch(final String contig, final long start , final long stop);

    /**
      * Query a specific interval on this reference, and get back an iterator over the bases spanning that interval.
      *
      * See the BaseUtils class for guidance on how to work with bases in this format.
      *
      * @param interval query interval
      * @return iterator over the bases spanning the query interval
      */
    @Override
    default public Iterator<Byte> query(final SimpleInterval interval) {
        // TODO: need a way to iterate lazily over reference bases without necessarily loading them all into memory at once
        return new ByteArrayIterator(queryAndPrefetch(interval).getBases());
    }

    /**
     * Get the sequence dictionary for this reference
     *
     * @return SAMSequenceDictionary for this reference
     */
    public SAMSequenceDictionary getSequenceDictionary();

    /**
     * Permanently close this data source. The default implementation does nothing.
     */
    @Override
    default public void close(){
        //do nothing
    }
}
