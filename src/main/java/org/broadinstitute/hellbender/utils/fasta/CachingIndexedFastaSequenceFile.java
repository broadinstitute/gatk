package org.broadinstitute.hellbender.utils.fasta;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.StringUtil;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;

/**
 * A caching version of the IndexedFastaSequenceFile that avoids going to disk as often as the raw indexer.
 *
 * Automatically upper-cases the bases coming in, unless the flag preserveCase is explicitly set.
 * Automatically converts IUPAC bases to Ns, unless the flag preserveIUPAC is explicitly set.
 */
public final class CachingIndexedFastaSequenceFile extends IndexedFastaSequenceFile {
    protected static final Logger logger = LogManager.getLogger(CachingIndexedFastaSequenceFile.class);

    /** do we want to print debugging information about cache efficiency? */
    private static final boolean PRINT_EFFICIENCY = false;

    /** If we are printing efficiency info, what frequency should we do it at? */
    private static final int PRINT_FREQUENCY = 10000;

    /** The default cache size in bp */
    public static final long DEFAULT_CACHE_SIZE = 1000000;

    /** The cache size of this CachingIndexedFastaSequenceFile */
    private final long cacheSize;

    /** When we have a cache miss at position X, we load sequence from X - cacheMissBackup */
    private final long cacheMissBackup;

    /**
     * If true, we will preserve the case of the original base in the genome
     */
    private final boolean preserveCase;

    /**
     * If true, we will preserve the IUPAC bases in the genome
     */
    private final boolean preserveIUPAC;

    // information about checking efficiency
    long cacheHits = 0;
    long cacheMisses = 0;

    /** Represents a specific cached sequence, with a specific start and stop, as well as the bases */
    private static class Cache {
        long start = -1, stop = -1;
        ReferenceSequence seq = null;
    }

    private final Cache cache = new Cache();

    /**
     * Same as general constructor but allows one to override the default cacheSize
     *
     * @param fasta the file we will read our FASTA sequence from.
     * @param index the index of the fasta file, used for efficient random access
     * @param cacheSize the size in bp of the cache we will use for this reader
     * @param preserveCase If true, we will keep the case of the underlying bases in the FASTA, otherwise everything is converted to upper case
     * @param preserveIUPAC If true, we will keep the IUPAC bases in the FASTA, otherwise they are converted to Ns
     */
    public CachingIndexedFastaSequenceFile(final File fasta, final FastaSequenceIndex index, final long cacheSize, final boolean preserveCase, final boolean preserveIUPAC) {
        super(fasta, index);
        if ( cacheSize < 0 ) throw new IllegalArgumentException("cacheSize must be > 0");
        this.cacheSize = cacheSize;
        this.cacheMissBackup = Math.max(cacheSize / 1000, 1);
        this.preserveCase = preserveCase;
        this.preserveIUPAC = preserveIUPAC;
    }

    /**
     * Open the given indexed fasta sequence file.  Throw an exception if the file cannot be opened.
     *
     * Looks for a index file for fasta on disk
     * Uses provided cacheSize instead of the default
     *
     * @param fasta The file to open.
     * @param cacheSize the size of the cache to use in this CachingIndexedFastaReader, must be >= 0
     * @param preserveCase If true, we will keep the case of the underlying bases in the FASTA, otherwise everything is converted to upper case
     */
    public CachingIndexedFastaSequenceFile(final File fasta, final long cacheSize, final boolean preserveCase, final boolean  preserveIUPAC) throws FileNotFoundException {
        super(fasta);
        if ( cacheSize < 0 ) throw new IllegalArgumentException("cacheSize must be > 0");
        this.cacheSize = cacheSize;
        this.cacheMissBackup = Math.max(cacheSize / 1000, 1);
        this.preserveCase = preserveCase;
        this.preserveIUPAC = preserveIUPAC;
    }

    /**
     * Same as general constructor but allows one to override the default cacheSize
     *
     * By default, this CachingIndexedFastaReader converts all incoming bases to upper case
     *
     * @param fasta the file we will read our FASTA sequence from.
     * @param index the index of the fasta file, used for efficient random access
     * @param cacheSize the size in bp of the cache we will use for this reader
     */
    public CachingIndexedFastaSequenceFile(final File fasta, final FastaSequenceIndex index, final long cacheSize) {
        this(fasta, index, cacheSize, false, false);
    }

    /**
     * Open the given indexed fasta sequence file.  Throw an exception if the file cannot be opened.
     *
     * Looks for a index file for fasta on disk.
     * This CachingIndexedFastaReader will convert all FASTA bases to upper cases under the hood
     *
     * @param fasta The file to open.
     */
    public CachingIndexedFastaSequenceFile(final File fasta) throws FileNotFoundException {
        this(fasta, false);
    }

    /**
     * Open the given indexed fasta sequence file.  Throw an exception if the file cannot be opened.
     *
     * Looks for a index file for fasta on disk
     *
     * @param fasta The file to open.
     * @param preserveCase If true, we will keep the case of the underlying bases in the FASTA, otherwise everything is converted to upper case
     */
    public CachingIndexedFastaSequenceFile(final File fasta, final boolean preserveCase) throws FileNotFoundException {
        this(fasta, DEFAULT_CACHE_SIZE, preserveCase, false);
    }

    /**
     * Create reference data source from fasta file, after performing several preliminary checks on the file.
     * This static utility was refactored from the constructor of ReferenceDataSource.
     * Possibly may be better as an overloaded constructor.
     * @param fastaFile Fasta file to be used as reference
     * @return A new instance of a CachingIndexedFastaSequenceFile.
     */
    public static CachingIndexedFastaSequenceFile checkAndCreate(final File fastaFile) {
        // does the fasta file exist? check that first...
        if (!fastaFile.exists()) {
            throw new UserException.MissingReference("The specified fasta file (" + fastaFile.getAbsolutePath() + ") does not exist.");
        }

        final boolean isGzipped = fastaFile.getAbsolutePath().endsWith(".gz");
        if ( isGzipped ) {
            throw new UserException.CannotHandleGzippedRef();
        }

        final File indexFile = new File(fastaFile.getAbsolutePath() + ".fai");

        // determine the name for the dict file
        final String fastaExt = fastaFile.getAbsolutePath().endsWith("fa") ? "\\.fa$" : "\\.fasta$";
        final File dictFile = new File(fastaFile.getAbsolutePath().replaceAll(fastaExt, ".dict"));

        // It's an error if either the fai or dict file does not exist. The user is now responsible
        // for creating these files.
        if (!indexFile.exists()) {
            throw new UserException.MissingReferenceFaiFile(indexFile, fastaFile);
        }
        if (!dictFile.exists()) {
            throw new UserException.MissingReferenceDictFile(dictFile, fastaFile);
        }

        // Read reference data by creating an IndexedFastaSequenceFile.
        try {
            return new CachingIndexedFastaSequenceFile(fastaFile);
        }
        catch (IllegalArgumentException e) {
            throw new UserException.CouldNotReadInputFile(fastaFile, "Could not read reference sequence.  The FASTA must have either a .fasta or .fa extension", e);
        }
        catch (Exception e) {
            throw new UserException.CouldNotReadInputFile(fastaFile, e);
        }
    }

    /**
     * Open the given indexed fasta sequence file.  Throw an exception if the file cannot be opened.
     *
     * Looks for a index file for fasta on disk
     * Uses provided cacheSize instead of the default
     *
     * @param fasta The file to open.
     * @param cacheSize the size of the cache to use in this CachingIndexedFastaReader, must be >= 0
     */
    public CachingIndexedFastaSequenceFile(final File fasta, final long cacheSize ) throws FileNotFoundException {
        this(fasta, cacheSize, false, false);
    }

    /**
     * Print the efficiency (hits / queries) to logger with priority
     */
    public void printEfficiency(final Level priority) {
        logger.log(priority, String.format("### CachingIndexedFastaReader: hits=%d misses=%d efficiency %.6f%%", cacheHits, cacheMisses, calcEfficiency()));
    }

    /**
     * Returns the efficiency (% of hits of all queries) of this object
     * @return
     */
    public double calcEfficiency() {
        return 100.0 * cacheHits / (cacheMisses + cacheHits * 1.0);
    }

    /**
     * @return the number of cache hits that have occurred
     */
    public long getCacheHits() {
        return cacheHits;
    }

    /**
     * @return the number of cache misses that have occurred
     */
    public long getCacheMisses() {
        return cacheMisses;
    }

    /**
     * @return the size of the cache we are using
     */
    public long getCacheSize() {
        return cacheSize;
    }

    /**
     * Is this CachingIndexedFastaReader keeping the original case of bases in the fasta, or is
     * everything being made upper case?
     *
     * @return true if the bases coming from this reader are in the original case in the fasta, false if they are all upper cased
     */
    public boolean isPreservingCase() {
        return preserveCase;
    }

    /**
     * Gets the subsequence of the contig in the range [start,stop]
     *
     * Uses the sequence cache if possible, or updates the cache to handle the request.  If the range
     * is larger than the cache itself, just loads the sequence directly, not changing the cache at all
     *
     * @param contig Contig whose subsequence to retrieve.
     * @param start inclusive, 1-based start of region.
     * @param stop inclusive, 1-based stop of region.
     * @return The partial reference sequence associated with this range.  If preserveCase is false, then
     *         all of the bases in the ReferenceSequence returned by this method will be upper cased.
     */
    @Override
    public ReferenceSequence getSubsequenceAt( final String contig, long start, final long stop ) {
        final ReferenceSequence result;

        if ( (stop - start) >= cacheSize ) {
            cacheMisses++;
            result = super.getSubsequenceAt(contig, start, stop);
            if ( ! preserveCase ) StringUtil.toUpperCase(result.getBases());
            if ( ! preserveIUPAC ) BaseUtils.convertIUPACtoN(result.getBases(), true, start < 1);
        } else {
            // todo -- potential optimization is to check if contig.name == contig, as this in general will be true
            SAMSequenceRecord contigInfo = super.getSequenceDictionary().getSequence(contig);
            if (contigInfo == null){
                throw new UserException.MissingContigInSequenceDictionary(contig, super.getSequenceDictionary());
            }

            if (stop > contigInfo.getSequenceLength())
                throw new SAMException("Query asks for data past end of contig. Query contig " + contig + " start:" + start + " stop:" + stop + " contigLength:" +  contigInfo.getSequenceLength());

            if ( start < cache.start || stop > cache.stop || cache.seq == null || cache.seq.getContigIndex() != contigInfo.getSequenceIndex() ) {
                cacheMisses++;
                cache.start = Math.max(start - cacheMissBackup, 0);
                cache.stop  = Math.min(start + cacheSize + cacheMissBackup, contigInfo.getSequenceLength());
                cache.seq   = super.getSubsequenceAt(contig, cache.start, cache.stop);

                // convert all of the bases in the sequence to upper case if we aren't preserving cases
                if ( ! preserveCase ) StringUtil.toUpperCase(cache.seq.getBases());
                if ( ! preserveIUPAC ) BaseUtils.convertIUPACtoN(cache.seq.getBases(), true, cache.start == 0);
            } else {
                cacheHits++;
            }

            // at this point we determine where in the cache we want to extract the requested subsequence
            final int cacheOffsetStart = (int)(start - cache.start);
            final int cacheOffsetStop = (int)(stop - start + cacheOffsetStart + 1);

            try {
                result = new ReferenceSequence(cache.seq.getName(), cache.seq.getContigIndex(), Arrays.copyOfRange(cache.seq.getBases(), cacheOffsetStart, cacheOffsetStop));
            } catch ( ArrayIndexOutOfBoundsException e ) {
                throw new GATKException(String.format("BUG: bad array indexing.  Cache start %d and end %d, request start %d end %d, offset start %d and end %d, base size %d",
                        cache.start, cache.stop, start, stop, cacheOffsetStart, cacheOffsetStop, cache.seq.getBases().length), e);
            }
        }

        // for debugging -- print out our efficiency if requested
        if ( PRINT_EFFICIENCY && (getCacheHits() + getCacheMisses()) % PRINT_FREQUENCY == 0 )
            printEfficiency(Level.INFO);

        return result;
    }
}