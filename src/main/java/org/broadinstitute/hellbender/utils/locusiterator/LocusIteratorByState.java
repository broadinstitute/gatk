package org.broadinstitute.hellbender.utils.locusiterator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;

/**
 * Iterator that traverses a SAM File, accumulating information on a per-locus basis
 *
 * Produces AlignmentContext objects, that contain ReadPileups of PileupElements.  This
 * class has its core job of converting an iterator of ordered GATKReads into those
 * ReadPileups.
 *
 * There are a few constraints on required and ensured by LIBS:
 *
 * -- Requires the Iterator<GATKRead> to returns reads in coordinate sorted order, consistent with the ordering
 * defined by the SAM file format.  That for performance reasons this constraint isn't actually enforced.
 * The behavior of LIBS is undefined in the case where the reads are badly ordered.
 * -- The reads in the ReadPileup are themselves in the order of appearance of the reads from the iterator.
 * That is, the pileup is ordered in a way consistent with the SAM coordinate ordering
 * -- Only aligned reads with at least one on-genomic cigar operator are passed on in the pileups.  That is,
 * unmapped reads or reads that are all insertions (10I) or soft clipped (10S) are not passed on.
 * -- LIBS can perform per-sample downsampling to a coverage target if requested.
 * -- Because of downsampling there's no guarantee that:
 *   -- A read that could be aligned to a position will actually occur in the pileup (downsampled away)
 *   -- A read that appears in a previous pileup that could align to a future position will actually occur
 *      in that pileup.  That is, a read might show up at position i but be downsampled away in the pileup at j
 * -- LIBS can optionally capture all of the reads that come off the iterator, before any leveling downsampling
 * occurs, if requested.  This allows users of LIBS to see both a ReadPileup view of the data as well as
 * a stream of unique, sorted reads
 */
public final class LocusIteratorByState implements Iterator<AlignmentContext> {
    /** Indicates that we shouldn't do any downsampling */
    public static final LIBSDownsamplingInfo NO_DOWNSAMPLING = new LIBSDownsamplingInfo(false, -1);

    /**
     * our log, which we want to capture anything from this class
     */
    private static final Logger logger = LogManager.getLogger(LocusIteratorByState.class);

    // -----------------------------------------------------------------------------------------------------------------
    //
    // member fields
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * A complete list of all samples that may come out of the reads.  Must be
     * comprehensive.
     */
    private final List<String> samples;

    /**
     * The system that maps incoming reads from the iterator to their pileup states
     */
    private final ReadStateManager readStates;

    /**
     * Should we include reads in the pileup which are aligned with a deletion operator to the reference?
     */
    private final boolean includeReadsWithDeletionAtLoci;

    /**
     * Should we include reads in the pileup which are aligned with a deletion operator to the reference?
     */
    private final boolean includeReadsWithNsAtLoci;

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
     * @param samIterator                    the iterator of reads to process into pileups.  Reads must be ordered
     *                                       according to standard coordinate-sorted BAM conventions
     * @param downsamplingMethod             information about how to downsample the reads
     * @param keepUniqueReadListInLIBS       if true, we will keep the unique reads from the samIterator and make them
     *                                       available via the transferReadsFromAllPreviousPileups interface
     * @param samples                        a complete list of samples present in the read groups for the reads coming from samIterator.
     *                                       This is generally just the set of read group sample fields in the SAMFileHeader.  This
     *                                       list of samples may contain a null element, and all reads without read groups will
     *                                       be mapped to this null sample
     * @param header                         header from the reads
     * @param includeReadsWithDeletionAtLoci Include reads with deletion at loci
     */
    public LocusIteratorByState(final Iterator<GATKRead> samIterator,
                                final DownsamplingMethod downsamplingMethod,
                                final boolean keepUniqueReadListInLIBS,
                                final Collection<String> samples,
                                final SAMFileHeader header,
                                final boolean includeReadsWithDeletionAtLoci) {
        this(samIterator,
                LIBSDownsamplingInfo.toDownsamplingInfo(downsamplingMethod),
                keepUniqueReadListInLIBS,
                samples,
                header,
                includeReadsWithDeletionAtLoci
        );
    }

    /**
     * Create a new LocusIteratorByState
     *
     * @param samIterator                    the iterator of reads to process into pileups.  Reads must be ordered
     *                                       according to standard coordinate-sorted BAM conventions
     * @param downsamplingMethod             information about how to downsample the reads
     * @param keepUniqueReadListInLIBS       if true, we will keep the unique reads from the samIterator and make them
     *                                       available via the transferReadsFromAllPreviousPileups interface
     * @param samples                        a complete list of samples present in the read groups for the reads coming from samIterator.
     *                                       This is generally just the set of read group sample fields in the SAMFileHeader.  This
     *                                       list of samples may contain a null element, and all reads without read groups will
     *                                       be mapped to this null sample
     * @param header                         header from the reads
     * @param includeReadsWithDeletionAtLoci Include reads with deletion at loci
     * @param includeReadsWithNsAtLoci       Include reads with Ns at loci (usually it is not needed)
     */
    public LocusIteratorByState(final Iterator<GATKRead> samIterator,
                                final DownsamplingMethod downsamplingMethod,
                                final boolean keepUniqueReadListInLIBS,
                                final Collection<String> samples,
                                final SAMFileHeader header,
                                final boolean includeReadsWithDeletionAtLoci,
                                final boolean includeReadsWithNsAtLoci) {
        this(samIterator,
                LIBSDownsamplingInfo.toDownsamplingInfo(downsamplingMethod),
                keepUniqueReadListInLIBS,
                samples,
                header,
                includeReadsWithDeletionAtLoci,
                includeReadsWithNsAtLoci
        );
    }

    /**
     * Create a new LocusIteratorByState
     *
     * @param samIterator                    the iterator of reads to process into pileups.  Reads must be ordered
     *                                       according to standard coordinate-sorted BAM conventions
     * @param downsamplingInfo               meta-information about how to downsample the reads
     * @param keepUniqueReadListInLIBS       if true, we will keep the unique reads from the samIterator and make them
     *                                       available via the transferReadsFromAllPreviousPileups interface
     * @param samples                        a complete list of samples present in the read groups for the reads coming from samIterator.
     *                                       This is generally just the set of read group sample fields in the SAMFileHeader.  This
     *                                       list of samples may contain a null element, and all reads without read groups will
     *                                       be mapped to this null sample
     * @param header                         header from the reads
     * @param includeReadsWithDeletionAtLoci Include reads with deletion at loci
     */
    public LocusIteratorByState(final Iterator<GATKRead> samIterator,
                                final LIBSDownsamplingInfo downsamplingInfo,
                                final boolean keepUniqueReadListInLIBS,
                                final Collection<String> samples,
                                final SAMFileHeader header,
                                final boolean includeReadsWithDeletionAtLoci) {
        this(samIterator,
                downsamplingInfo,
                keepUniqueReadListInLIBS,
                samples,
                header,
                includeReadsWithDeletionAtLoci,
                false
        );
    }

    /**
     * Create a new LocusIteratorByState
     *
     * @param samIterator                    the iterator of reads to process into pileups.  Reads must be ordered
     *                                       according to standard coordinate-sorted BAM conventions
     * @param downsamplingInfo               meta-information about how to downsample the reads
     * @param keepUniqueReadListInLIBS       if true, we will keep the unique reads from the samIterator and make them
     *                                       available via the transferReadsFromAllPreviousPileups interface
     * @param samples                        a complete list of samples present in the read groups for the reads coming from samIterator.
     *                                       This is generally just the set of read group sample fields in the SAMFileHeader.  This
     *                                       list of samples may contain a null element, and all reads without read groups will
     *                                       be mapped to this null sample
     * @param header                         header from the reads
     * @param includeReadsWithDeletionAtLoci Include reads with deletion at loci
     * @param includeReadsWithNsAtLoci       Include reads with Ns at loci (usually it is not needed)
     */
    public LocusIteratorByState(final Iterator<GATKRead> samIterator,
                                final LIBSDownsamplingInfo downsamplingInfo,
                                final boolean keepUniqueReadListInLIBS,
                                final Collection<String> samples,
                                final SAMFileHeader header,
                                final boolean includeReadsWithDeletionAtLoci,
                                final boolean includeReadsWithNsAtLoci) {
        Utils.nonNull(samIterator, "samIterator cannot be null");
        Utils.nonNull(downsamplingInfo, "downsamplingInfo cannot be null");
        Utils.nonNull(samples, "Samples cannot be null");
        Utils.nonNull(header, "header cannot be null");

        // currently the GATK expects this LocusIteratorByState to accept empty sample lists, when
        // there's no read data.  So we need to throw this error only when samIterator.hasNext() is true
        if (samples.isEmpty() && samIterator.hasNext()) {
            throw new IllegalArgumentException("samples list must not be empty");
        }

        this.includeReadsWithDeletionAtLoci = includeReadsWithDeletionAtLoci;
        this.includeReadsWithNsAtLoci = includeReadsWithNsAtLoci;
        this.samples = new ArrayList<>(samples);
        this.readStates = new ReadStateManager(samIterator, this.samples, downsamplingInfo, keepUniqueReadListInLIBS, header);
    }

    /**
     * Get the current location (i.e., the bp of the center of the pileup) of the pileup, or null if not anywhere yet
     *
     * Assumes that read states is updated to reflect the current pileup position, but not advanced to the
     * next location.
     *
     * @return the location of the current pileup, or null if we're after all reads
     */
    private SimpleInterval getLocation() {
        return readStates.isEmpty() ? null : readStates.getFirst().getLocation();
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // next() routine and associated collection operations
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Is there another pileup available?
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
        if (!hasNext()) {
            throw new NoSuchElementException("LocusIteratorByState: out of elements.");
        }
        AlignmentContext currentAlignmentContext = nextAlignmentContext;
        nextAlignmentContext = null;
        return currentAlignmentContext;
    }

    /**
     * Move this LIBS until we are over position
     *
     * Will return null if cannot reach position (because we run out of data in the locus)
     *
     * @param position                             the start position of the AlignmentContext we want back
     * @param stopAtFirstNonEmptySiteAfterPosition if true, we will stop as soon as we find a context with data with
     *                                             position >= position, otherwise we will return a null value
     *                                             and consume the data for the next position.  This means that without
     *                                             specifying this value the LIBS will be in an indeterminate state
     *                                             after calling this function, and should be reconstructed from scratch
     *                                             for subsequent use
     * @return a AlignmentContext at position, or null if this isn't possible
     */
    public AlignmentContext advanceToLocus(final int position, final boolean stopAtFirstNonEmptySiteAfterPosition) {
        while (hasNext()) {
            final AlignmentContext context = next();

            if (context == null) {
                // we ran out of data
                return null;
            }

            if (context.getPosition() == position) {
                return context;
            }

            if (context.getPosition() > position) {
                return stopAtFirstNonEmptySiteAfterPosition ? context : null;
            }
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

            final Locatable location = getLocation();

            // We don't need to keep the pileup elements separated by sample within this method,
            // since they are just going to get combined into one monolithic pileup anyway
            // when we construct the final ReadPileup below. This optimization speeds up the
            // HaplotypeCaller by quite a bit!
            final List<PileupElement> allPileupElements = new ArrayList<>(100);

            for (final Map.Entry<String, PerSampleReadStateManager> sampleStatePair : readStates) {
                final PerSampleReadStateManager readState = sampleStatePair.getValue();
                final Iterator<AlignmentStateMachine> iterator = readState.iterator();

                while (iterator.hasNext()) {
                    // state object with the read/offset information
                    final AlignmentStateMachine state = iterator.next();
                    final GATKRead read = state.getRead();
                    final CigarOperator op = state.getCigarOperator();

                    if (!includeReadsWithNsAtLoci && op == CigarOperator.N) {
                        continue;
                    }

                    if (!dontIncludeReadInPileup(read, location.getStart())) {
                        if (!includeReadsWithDeletionAtLoci && op == CigarOperator.D) {
                            continue;
                        }

                        allPileupElements.add(state.makePileupElement());
                    }
                }
            }

            readStates.updateReadStates(); // critical - must be called after we get the current state offsets and location
            if (!allPileupElements.isEmpty()) { // if we got reads with non-D/N over the current position, we are done
                nextAlignmentContext = new AlignmentContext(location, new ReadPileup(location, allPileupElements));
            }
        }
    }

    /**
     * Should this read be excluded from the pileup?
     *
     * Generic place to put per-base filters appropriate to LocusIteratorByState
     *
     * @param rec the read to potentially exclude
     * @param pos the genomic position of the current alignment
     * @return true if the read should be excluded from the pileup, false otherwise
     */
    private boolean dontIncludeReadInPileup(final GATKRead rec, final long pos) {
        return ReadUtils.isBaseInsideAdaptor(rec, pos);
    }

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
     * underlying GATKRead iterator and that appeared at any point in the list of SAMRecordAlignmentState for
     * any reads.  This function is intended to allow users to efficiently reconstruct the unique set of reads
     * used across all pileups.  This is necessary for LIBS to handle because attempting to do
     * so from the pileups coming out of LIBS is extremely expensive.
     *
     * This functionality is only available if LIBS was created with the argument to track the reads
     *
     * @return the current list
     * @throws UnsupportedOperationException if called when keepingSubmittedReads is false
     */
    public List<GATKRead> transferReadsFromAllPreviousPileups() {
        return readStates.transferSubmittedReads();
    }

    /**
     * Get the underlying list of tracked reads.  For testing only
     * @return a non-null list
     */
    @VisibleForTesting
    List<GATKRead> getReadsFromAllPreviousPileups() {
        return readStates.getSubmittedReads();
    }
}
