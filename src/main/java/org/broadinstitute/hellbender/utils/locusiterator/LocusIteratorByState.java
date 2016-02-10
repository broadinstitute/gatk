/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.locusiterator;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecordIterator;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;

import java.util.*;

/**
 * Iterator that traverses a SAM File, accumulating information on a per-locus basis
 *
 * Produces AlignmentContext objects, that contain ReadBackedPileups of PileupElements.  This
 * class has its core job of converting an iterator of ordered SAMRecords into those
 * RBPs.
 *
 * There are a few constraints on required and ensured by LIBS:
 *
 * -- Requires the Iterator<GATKSAMRecord> to returns reads in coordinate sorted order, consistent with the ordering
 * defined by the SAM file format.  That for performance reasons this constraint isn't actually enforced.
 * The behavior of LIBS is undefined in the case where the reads are badly ordered.
 * -- The reads in the ReadBackedPileup are themselves in the order of appearance of the reads from the iterator.
 * That is, the pileup is ordered in a way consistent with the SAM coordinate ordering
 * -- Only aligned reads with at least one on-genomic cigar operator are passed on in the pileups.  That is,
 * unmapped reads or reads that are all insertions (10I) or soft clipped (10S) are not passed on.
 * -- LIBS can perform per-sample downsampling of a variety of kinds.
 * -- Because of downsampling there's no guarantee that:
 *   -- A read that could be aligned to a position will actually occur in the pileup (downsampled away)
 *   -- A read that appears in a previous pileup that could align to a future position will actually occur
 *      in that pileup.  That is, a read might show up at position i but be downsampled away in the pileup at j
 * -- LIBS can optionally capture all of the reads that come off the iterator, before any leveling downsampling
 * occurs, if requested.  This allows users of LIBS to see both a ReadBackedPileup view of the data as well as
 * a stream of unique, sorted reads
 */
public final class LocusIteratorByState extends LocusIterator {
    /** Indicates that we shouldn't do any downsampling */
    public final static LIBSDownsamplingInfo NO_DOWNSAMPLING = new LIBSDownsamplingInfo(false, -1);

    /**
     * our log, which we want to capture anything from this class
     */
    private final static Logger logger = Logger.getLogger(LocusIteratorByState.class);

    // -----------------------------------------------------------------------------------------------------------------
    //
    // member fields
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Used to create new GenomeLocs as needed
     */
    private final GenomeLocParser genomeLocParser;

    /**
     * A complete list of all samples that may come out of the reads.  Must be
     * comprehensive.
     */
    private final ArrayList<String> samples;

    /**
     * The system that maps incoming reads from the iterator to their pileup states
     */
    private final ReadStateManager readStates;

    /**
     * Should we include reads in the pileup which are aligned with a deletion operator to the reference?
     */
    private final boolean includeReadsWithDeletionAtLoci;

    /**
     * The next alignment context.  A non-null value means that a
     * context is waiting from hasNext() for sending off to the next next() call.  A null
     * value means that either hasNext() has not been called at all or that
     * the underlying iterator is exhausted
     */
    private AlignmentContext nextAlignmentContext;

    // -----------------------------------------------------------------------------------------------------------------
    //
    // constructors and other basic operations
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Create a new LocusIteratorByState
     *
     * @param samIterator the iterator of reads to process into pileups.  Reads must be ordered
     *                    according to standard coordinate-sorted BAM conventions
     * @param downsamplingMethod information about how to downsample the reads
     * @param includeReadsWithDeletionAtLoci Include reads with deletion at loci
     * @param keepUniqueReadListInLIBS Keep unique read list in LIBS
     * @param genomeLocParser used to create genome locs
     * @param samples a complete list of samples present in the read groups for the reads coming from samIterator.
     *                This is generally just the set of read group sample fields in the SAMFileHeader.  This
     *                list of samples may contain a null element, and all reads without read groups will
     *                be mapped to this null sample
     */
    public LocusIteratorByState(final Iterator<GATKSAMRecord> samIterator,
                                final DownsamplingMethod downsamplingMethod,
                                final boolean includeReadsWithDeletionAtLoci,
                                final boolean keepUniqueReadListInLIBS,
                                final GenomeLocParser genomeLocParser,
                                final Collection<String> samples) {
        this(samIterator,
                toDownsamplingInfo(downsamplingMethod),
                includeReadsWithDeletionAtLoci,
                genomeLocParser,
                samples,
                keepUniqueReadListInLIBS);
    }

    /**
     * Create a new LocusIteratorByState based on a SAMFileReader using reads in an iterator it
     *
     * Simple constructor that uses the samples in the reader, doesn't do any downsampling,
     * and makes a new GenomeLocParser using the reader.  This constructor will be slow(ish)
     * if you continually invoke this constructor, but it's easy to make.
     *
     * @param reader a non-null reader
     * @param it an iterator from reader that has the reads we want to use to create ReadBackPileups
     */
    public LocusIteratorByState(final SAMFileReader reader, final CloseableIterator<SAMRecord> it) {
        this(new GATKSAMRecordIterator(it),
                new LIBSDownsamplingInfo(false, 0),
                true,
                new GenomeLocParser(reader.getFileHeader().getSequenceDictionary()),
                ReadUtils.getSAMFileSamples(reader.getFileHeader()),
                false);
    }

    /**
     * Create a new LocusIteratorByState
     *
     * @param samIterator the iterator of reads to process into pileups.  Reads must be ordered
     *                    according to standard coordinate-sorted BAM conventions
     * @param downsamplingInfo meta-information about how to downsampling the reads
     * @param genomeLocParser used to create genome locs
     * @param samples a complete list of samples present in the read groups for the reads coming from samIterator.
     *                This is generally just the set of read group sample fields in the SAMFileHeader.  This
     *                list of samples may contain a null element, and all reads without read groups will
     *                be mapped to this null sample
     * @param maintainUniqueReadsList if true, we will keep the unique reads from off the samIterator and make them
     *                                available via the transferReadsFromAllPreviousPileups interface
     */
    public LocusIteratorByState(final Iterator<GATKSAMRecord> samIterator,
                                final LIBSDownsamplingInfo downsamplingInfo,
                                final boolean includeReadsWithDeletionAtLoci,
                                final GenomeLocParser genomeLocParser,
                                final Collection<String> samples,
                                final boolean maintainUniqueReadsList) {
        if ( samIterator == null ) throw new IllegalArgumentException("samIterator cannot be null");
        if ( downsamplingInfo == null ) throw new IllegalArgumentException("downsamplingInfo cannot be null");
        if ( genomeLocParser == null ) throw new IllegalArgumentException("genomeLocParser cannot be null");
        if ( samples == null ) throw new IllegalArgumentException("Samples cannot be null");

        // currently the GATK expects this LocusIteratorByState to accept empty sample lists, when
        // there's no read data.  So we need to throw this error only when samIterator.hasNext() is true
        if (samples.isEmpty() && samIterator.hasNext()) {
            throw new IllegalArgumentException("samples list must not be empty");
        }

        this.genomeLocParser = genomeLocParser;
        this.includeReadsWithDeletionAtLoci = includeReadsWithDeletionAtLoci;
        this.samples = new ArrayList<String>(samples);
        this.readStates = new ReadStateManager(samIterator, this.samples, downsamplingInfo, maintainUniqueReadsList);
    }

    @Override
    public Iterator<AlignmentContext> iterator() {
        return this;
    }

    /**
     * Get the current location (i.e., the bp of the center of the pileup) of the pileup, or null if not anywhere yet
     *
     * Assumes that read states is updated to reflect the current pileup position, but not advanced to the
     * next location.
     *
     * @return the location of the current pileup, or null if we're after all reads
     */
    private GenomeLoc getLocation() {
        return readStates.isEmpty() ? null : readStates.getFirst().getLocation(genomeLocParser);
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // next() routine and associated collection operations
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Is there another pileup available?
     * @return
     */
    @Override
    public boolean hasNext() {
        lazyLoadNextAlignmentContext();
        return nextAlignmentContext != null;
    }

    /**
     * Get the next AlignmentContext available from the reads.
     *
     * @return a non-null AlignmentContext of the pileup after to the next genomic position covered by
     * at least one read.
     */
    @Override
    public AlignmentContext next() {
        lazyLoadNextAlignmentContext();
        if (!hasNext())
            throw new NoSuchElementException("LocusIteratorByState: out of elements.");
        AlignmentContext currentAlignmentContext = nextAlignmentContext;
        nextAlignmentContext = null;
        return currentAlignmentContext;
    }

    /**
     * Move this LIBS until we are over position
     *
     * Will return null if cannot reach position (because we run out of data in the locus)
     *
     * @param position the start position of the AlignmentContext we want back
     * @param stopAtFirstNonEmptySiteAfterPosition if true, we will stop as soon as we find a context with data with
     *                                             position >= position, otherwise we will return a null value
     *                                             and consume the data for the next position.  This means that without
     *                                             specifying this value the LIBS will be in an indeterminate state
     *                                             after calling this function, and should be reconstructed from scratch
     *                                             for subsequent use
     * @return a AlignmentContext at position, or null if this isn't possible
     */
    public AlignmentContext advanceToLocus(final int position, final boolean stopAtFirstNonEmptySiteAfterPosition) {
        while ( hasNext() ) {
            final AlignmentContext context = next();

            if ( context == null )
                // we ran out of data
                return null;

            if ( context.getPosition() == position )
                return context;

            if ( context.getPosition() > position)
                return stopAtFirstNonEmptySiteAfterPosition ? context : null;
        }

        return null;
    }

    /**
     * Creates the next alignment context from the given state.  Note that this is implemented as a
     * lazy load method. nextAlignmentContext MUST BE null in order for this method to advance to the
     * next entry.
     */
    private void lazyLoadNextAlignmentContext() {
        while (nextAlignmentContext == null && readStates.hasNext()) {
            readStates.collectPendingReads();

            final GenomeLoc location = getLocation();
            final Map<String, ReadBackedPileupImpl> fullPileup = new HashMap<String, ReadBackedPileupImpl>();

            for (final Map.Entry<String, PerSampleReadStateManager> sampleStatePair : readStates ) {
                final String sample = sampleStatePair.getKey();
                final PerSampleReadStateManager readState = sampleStatePair.getValue();
                final Iterator<AlignmentStateMachine> iterator = readState.iterator();
                final List<PileupElement> pile = new ArrayList<PileupElement>(readState.size());

                while (iterator.hasNext()) {
                    // state object with the read/offset information
                    final AlignmentStateMachine state = iterator.next();
                    final GATKSAMRecord read = state.getRead();
                    final CigarOperator op = state.getCigarOperator();

                    if (op == CigarOperator.N) // N's are never added to any pileup
                        continue;

                    if (!dontIncludeReadInPileup(read, location.getStart())) {
                        if ( ! includeReadsWithDeletionAtLoci && op == CigarOperator.D ) {
                            continue;
                        }

                        pile.add(state.makePileupElement());
                    }
                }

                if (! pile.isEmpty() ) // if this pileup added at least one base, add it to the full pileup
                    fullPileup.put(sample, new ReadBackedPileupImpl(location, pile));
            }

            readStates.updateReadStates(); // critical - must be called after we get the current state offsets and location
            if (!fullPileup.isEmpty()) // if we got reads with non-D/N over the current position, we are done
                nextAlignmentContext = new AlignmentContext(location, new ReadBackedPileupImpl(location, fullPileup), false);
        }
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // getting the list of reads
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Transfer current list of all unique reads that have ever been used in any pileup, clearing old list
     *
     * This list is guaranteed to only contain unique reads, even across calls to the this function.  It is
     * literally the unique set of reads ever seen.
     *
     * The list occurs in the same order as they are encountered in the underlying iterator.
     *
     * Takes the maintained list of submitted reads, and transfers it to the caller of this
     * function.  The old list of set to a new, cleanly allocated list so the caller officially
     * owns the list returned by this call.  This is the only way to clear the tracking
     * of submitted reads, if enabled.
     *
     * The purpose of this function is allow users of LIBS to keep track of all of the reads pulled off the
     * underlying GATKSAMRecord iterator and that appeared at any point in the list of SAMRecordAlignmentState for
     * any reads.  This function is intended to allow users to efficiently reconstruct the unique set of reads
     * used across all pileups.  This is necessary for LIBS to handle because attempting to do
     * so from the pileups coming out of LIBS is extremely expensive.
     *
     * This functionality is only available if LIBS was created with the argument to track the reads
     *
     * @throws UnsupportedOperationException if called when keepingSubmittedReads is false
     *
     * @return the current list
     */
    @Ensures("result != null")
    public List<GATKSAMRecord> transferReadsFromAllPreviousPileups() {
        return readStates.transferSubmittedReads();
    }

    /**
     * Get the underlying list of tracked reads.  For testing only
     * @return a non-null list
     */
    @Ensures("result != null")
    protected List<GATKSAMRecord> getReadsFromAllPreviousPileups() {
        return readStates.getSubmittedReads();
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // utility functions
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Should this read be excluded from the pileup?
     *
     * Generic place to put per-base filters appropriate to LocusIteratorByState
     *
     * @param rec the read to potentially exclude
     * @param pos the genomic position of the current alignment
     * @return true if the read should be excluded from the pileup, false otherwise
     */
    @Requires({"rec != null", "pos > 0"})
    private boolean dontIncludeReadInPileup(final GATKSAMRecord rec, final long pos) {
        return ReadUtils.isBaseInsideAdaptor(rec, pos);
    }

    /**
     * Create a LIBSDownsamplingInfo object from the requested info in DownsamplingMethod
     *
     * LIBS will invoke the Reservoir and Leveling downsamplers on the read stream if we're
     * downsampling to coverage by sample. SAMDataSource will have refrained from applying
     * any downsamplers to the read stream in this case, in the expectation that LIBS will
     * manage the downsampling. The reason for this is twofold: performance (don't have to
     * split/re-assemble the read stream in SAMDataSource), and to enable partial downsampling
     * of reads (eg., using half of a read, and throwing the rest away).
     *
     * @param downsamplingMethod downsampling information about what should be done to the reads
     * @return a LIBS specific info holder about downsampling only
     */
    @Requires("downsamplingMethod != null")
    @Ensures("result != null")
    private static LIBSDownsamplingInfo toDownsamplingInfo(final DownsamplingMethod downsamplingMethod) {
        final boolean performDownsampling = downsamplingMethod != null &&
                downsamplingMethod.type == DownsampleType.BY_SAMPLE &&
                downsamplingMethod.toCoverage != null;
        final int coverage = performDownsampling ? downsamplingMethod.toCoverage : 0;

        return new LIBSDownsamplingInfo(performDownsampling, coverage);
    }

    /**
     * Create a pileup element for read at offset
     *
     * offset must correspond to a valid read offset given the read's cigar, or an IllegalStateException will be throw
     *
     * @param read a read
     * @param offset the offset into the bases we'd like to use in the pileup
     * @return a valid PileupElement with read and at offset
     */
    @Ensures("result != null")
    public static PileupElement createPileupForReadAndOffset(final GATKSAMRecord read, final int offset) {
        if ( read == null ) throw new IllegalArgumentException("read cannot be null");
        if ( offset < 0 || offset >= read.getReadLength() ) throw new IllegalArgumentException("Invalid offset " + offset + " outside of bounds 0 and " + read.getReadLength());

        final AlignmentStateMachine stateMachine = new AlignmentStateMachine(read);

        while ( stateMachine.stepForwardOnGenome() != null ) {
            if ( stateMachine.getReadOffset() == offset )
                return stateMachine.makePileupElement();
        }

        throw new IllegalStateException("Tried to create a pileup for read " + read + " with offset " + offset +
                " but we never saw such an offset in the alignment state machine");
    }

    /**
     * For testing only.  Assumes that the incoming SAMRecords have no read groups, so creates a dummy sample list
     * for the system.
     */
    public static List<String> sampleListForSAMWithoutReadGroups() {
        List<String> samples = new ArrayList<String>();
        samples.add(null);
        return samples;
    }
}