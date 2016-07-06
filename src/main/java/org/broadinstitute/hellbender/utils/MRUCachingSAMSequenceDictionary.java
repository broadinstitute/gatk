package org.broadinstitute.hellbender.utils;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * A wrapper class that provides efficient most recently used caching for the global
 * SAMSequenceDictionary underlying all of the GATK engine capabilities.  It is essential
 * that these class be as efficient as possible.  It doesn't need to be thread-safe, as
 * GenomeLocParser uses a thread-local variable to ensure that each thread gets its own MRU
 * cache.
 *
 * The MRU elements are the SAMSequenceRecord, the lastContig, and the lastIndex.  The
 * cached value is the actual SAMSequenceRecord of the most recently accessed value from
 * getSequence, along with local variables for the contig index and contig string.
 */
final class MRUCachingSAMSequenceDictionary {
    /**
     * Our sequence dictionary
     */
    private final SAMSequenceDictionary dict;

    SAMSequenceRecord lastSSR = null;
    String lastContig = "";
    int lastIndex = -1;

    /**
     * Create a new MRUCachingSAMSequenceDictionary that provides information about sequences in dict
     * @param dict a non-null, non-empty sequencing dictionary
     */
    public MRUCachingSAMSequenceDictionary(final SAMSequenceDictionary dict) {
        Utils.nonNull( dict == null, "Dictionary cannot be null");
        Utils.validateArg( !dict.isEmpty(), "Dictionary cannot have size zero");
        this.dict = dict;
    }

    /**
     * Get our sequence dictionary
     * @return a non-null SAMSequenceDictionary
     */
    public SAMSequenceDictionary getDictionary() {
        return dict;
    }

    /**
     * Is contig present in the dictionary?  Efficiently caching.
     * @param contig a non-null contig we want to test
     * @return true if contig is in dictionary, false otherwise
     */
    public final boolean hasContig(final String contig) {
        return contig.equals(lastContig) || dict.getSequence(contig) != null;
    }

    /**
     * Same as SAMSequenceDictionary.getSequence but uses a MRU cache for efficiency
     *
     * @param contig the contig name we want to get the sequence record of
     * @throws GATKException if contig isn't present in the dictionary
     * @return the sequence record for contig
     */
    public final SAMSequenceRecord getSequence(final String contig) {
        if ( isCached(contig) )
            return lastSSR;
        else
            return updateCache(contig, -1);
    }

    /**
     * Same as SAMSequenceDictionary.getSequence but uses a MRU cache for efficiency
     *
     * @param index the contig index we want to get the sequence record of
     * @throws GATKException if contig isn't present in the dictionary
     * @return the sequence record for contig
     */
    public final SAMSequenceRecord getSequence(final int index) {
        if ( isCached(index) )
            return lastSSR;
        else
            return updateCache(null, index);
    }

    /**
     * Same as SAMSequenceDictionary.getSequenceIndex but uses a MRU cache for efficiency
     *
     * @param contig the contig we want to get the sequence record of
     * @throws GATKException if index isn't present in the dictionary
     * @return the sequence record index for contig
     */
    public final int getSequenceIndex(final String contig) {
        if ( ! isCached(contig) ) {
            updateCache(contig, -1);
        }

        return lastIndex;
    }

    /**
     * Is contig the MRU cached contig?
     * @param contig the contig to test
     * @return true if contig is the currently cached contig, false otherwise
     */
    protected boolean isCached(final String contig) {
        return contig.equals(lastContig);
    }

    /**
     * Is the contig index index the MRU cached index?
     * @param index the contig index to test
     * @return true if contig index is the currently cached contig index, false otherwise
     */
    protected boolean isCached(final int index) {
        return lastIndex == index;
    }

    /**
     * The key algorithm.  Given a new record, update the last used record, contig
     * name, and index.
     *
     * @param contig the contig we want to look up.  If null, index is used instead
     * @param index the contig index we want to look up.  Only used if contig is null
     * @throws GATKException if index isn't present in the dictionary
     * @return the SAMSequenceRecord for contig / index
     */
    private SAMSequenceRecord updateCache(final String contig, int index ) {
        SAMSequenceRecord rec = contig == null ? dict.getSequence(index) : dict.getSequence(contig);
        if ( rec == null ) {
            throw new GATKException("BUG: requested unknown contig=" + contig + " index=" + index);
        } else {
            lastSSR = rec;
            lastContig = rec.getSequenceName();
            lastIndex = rec.getSequenceIndex();
            return rec;
        }
    }
}
