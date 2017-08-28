package org.broadinstitute.hellbender.tools.walkers.indels;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.*;

import java.util.*;

/**
 * A locally resorting, mate fixing sam file writer that supports an idiom where reads are only moved around if
 * the ISIZE of the pair is < X and reads are not allowed to move any more than Y bp from their original positions.
 *
 * To understand this data structure, let's begin by asking -- when are we certain we know the position of read R added
 * to the writer and its mate M given that R has been added to the writer (but M may not be), their ISIZE in R, at the
 * moment that a read K is added to the writer, under the constraints X and Y?  Complex I know.  First, because
 * reads cannot move more than Y bp in either direction, we know that R originated at most R.pos + Y bp from its
 * current position.  Also, we know that K is at most K.pos + Y bp from it's original position.  If R is maximally
 * shifted to the right, and K shifted to the left, then they could at most move 2Y together.  So if the distance
 * between R and K > 2Y, we know that there are no reads left in the original stream that could be moved before R.
 *
 * Now, we also need to be certain if we have a mate pair M, that won't emit R before we can incorporate any move of
 * M into the mate pair info R.  There are two cases to consider here:
 *
 * If ISIZE > X, we know that we won't move M when we see it, so we can safely emit R knowing that
 * M is fixed in place.
 *
 * If ISIZE <= X, M might be moved, and it we have to wait until we see M in the stream to know it's position.
 * So R must be buffered until either M arrives, or we see a read K that's more than 2Y units past the original position
 * of M.
 *
 * So the worst-case memory consumption here is proportional to the number of reads
 * occurring between R and M + 2 Y, and so is proportional to the depth of the data and X and Y.
 *
 * This leads to the following simple algorithm:
 *
 * addAlignment(newRead):
 *   addReadToListOfReads(newRead)
 *   update mate pair of newRead if present in list of reads
 *
 *   for ( read in list of reads [in order of increasing read.pos] ):
 *     if read.pos < newRead.pos - 2Y && (read.isize >= X || read.matePos < newRead.pos - 2 * Y):
 *        emit read and remove from list of reads
 *     else:
 *        break
 *
 * @author depristo, ebanks
 * @version 0.2
 */
public class ConstrainedMateFixingManager {

    protected static final Logger logger = Logger.getLogger(ConstrainedMateFixingManager.class);
    private static final boolean DEBUG = false;

    /** How often do we check whether we want to emit reads? */
    protected final static int EMIT_FREQUENCY = 1000;

    /**
     * How much could a single read move in position from its original position?
     */
    final int MAX_POS_MOVE_ALLOWED;

    /**
     * How many reads should we store in memory before flushing the queue?
     */
    final int MAX_RECORDS_IN_MEMORY;

    // header for the comparator and for access information about the read
    private final SAMFileHeader header;
    /** how we order our SAM records */
    private final ReadCoordinateComparator comparer;

    /** The place where we ultimately write out our records */
    final GATKReadWriter writer;

    /**
     * what is the maximum isize of a pair of reads that can move?  Reads with isize > this value
     * are assumes to not be allowed to move in the incoming read stream.
     */
    final int maxInsertSizeForMovingReadPairs;
    final int initialCapacity = 5000;

    private SimpleInterval lastLocFlushed = null;

    int counter = 0;

    /** read.name -> records */
    HashMap<String, GATKReadHashObject> forMateMatching = new HashMap<String, GATKReadHashObject>();
    PriorityQueue<GATKRead> waitingReads;

    private GATKRead remove(PriorityQueue<GATKRead> queue) {
        GATKRead first = queue.poll();
        if (first == null)
            throw new UserException("Error caching SAM record -- priority queue is empty, and yet there was an attempt to poll it -- which is usually caused by malformed SAM/BAM files in which multiple identical copies of a read are present.");
        return first;
    }

    private static class GATKReadHashObject {
        public GATKRead record;
        public boolean wasModified;

        public GATKReadHashObject(GATKRead record, boolean wasModified) {
            this.record = record;
            this.wasModified = wasModified;
        }
    }


    //private SimpleTimer timer = new SimpleTimer("ConstrainedWriter");
    //private long PROGRESS_PRINT_FREQUENCY = 10 * 1000;             // in milliseconds
    //private long lastProgressPrintTime = -1;                       // When was the last time we printed progress log?


    /**
     *
     * @param writer                                 actual writer
     * @param header                                 the header object
     * @param maxInsertSizeForMovingReadPairs        max insert size allowed for moving pairs
     * @param maxMoveAllowed                         max positional move allowed for any read
     * @param maxRecordsInMemory                     max records to keep in memory
     */
    public ConstrainedMateFixingManager(final GATKReadWriter writer,
                                        final SAMFileHeader header,
                                        final int maxInsertSizeForMovingReadPairs,
                                        final int maxMoveAllowed,
                                        final int maxRecordsInMemory) {
        this.writer = writer;
        this.header = header;
        this.maxInsertSizeForMovingReadPairs = maxInsertSizeForMovingReadPairs;
        this.MAX_POS_MOVE_ALLOWED = maxMoveAllowed;
        this.MAX_RECORDS_IN_MEMORY = maxRecordsInMemory;

        // initialize comparer
        this.comparer = new ReadCoordinateComparator(this.header);
        this.waitingReads = new PriorityQueue<>(initialCapacity, comparer);

        //timer.start();
        //lastProgressPrintTime = timer.currentTime();
    }

    public int getNReadsInQueue() { return waitingReads.size(); }

    /**
     * For testing purposes only
     *
     * @return the list of reads currently in the queue
     */
    protected List<GATKRead> getReadsInQueueForTesting() {
        return new ArrayList<>(waitingReads);
    }

    public boolean canMoveReads(Locatable earliestPosition) {
        if ( DEBUG ) logger.info("Refusing to realign? " + earliestPosition + " vs. " + lastLocFlushed);

        return lastLocFlushed == null ||
                ! lastLocFlushed.getContig().equals(earliestPosition.getContig()) ||
                ! lastLocFlushed.overlapsWithMargin(earliestPosition, maxInsertSizeForMovingReadPairs);
    }

    private boolean noReadCanMoveBefore(int pos, GATKRead addedRead) {
        return pos + 2 * MAX_POS_MOVE_ALLOWED < addedRead.getAssignedStart();
    }

    public void addRead(GATKRead newRead, boolean readWasModified) {
        addRead(newRead, readWasModified, true);
    }

    public void addReads(List<GATKRead> newReads, Set<GATKRead> modifiedReads) {
        for ( GATKRead newRead : newReads )
            addRead(newRead, modifiedReads.contains(newRead), false);
    }

    protected void addRead(GATKRead newRead, boolean readWasModified, boolean canFlush) {
        if ( DEBUG ) logger.info("New read pos " + newRead.getAssignedStart() + " OP = " + newRead.getAttributeAsString("OP") + " " + readWasModified);

        //final long curTime = timer.currentTime();
        //if ( curTime - lastProgressPrintTime > PROGRESS_PRINT_FREQUENCY ) {
        //    lastProgressPrintTime = curTime;
        //    System.out.println("WaitingReads.size = " + waitingReads.size() + ", forMateMatching.size = " + forMateMatching.size());
        //}

        // if the new read is on a different contig or we have too many reads, then we need to flush the queue and clear the map
        boolean tooManyReads = getNReadsInQueue() >= MAX_RECORDS_IN_MEMORY;
        if ( (canFlush && tooManyReads) || (getNReadsInQueue() > 0 && ReadUtils.getAssignedReferenceIndex(waitingReads.peek(), header) != ReadUtils.getAssignedReferenceIndex(newRead, header)) ) {
            if ( DEBUG ) logger.warn("Flushing queue on " + (tooManyReads ? "too many reads" : ("move to new contig: " + newRead.getAssignedContig() + " from " + waitingReads.peek().getAssignedContig())) + " at " + newRead.getAssignedStart());

            while ( getNReadsInQueue() > 1 ) {
                // emit to disk
                writeRead(remove(waitingReads));
            }

            GATKRead lastRead = remove(waitingReads);
            lastLocFlushed = ReadUtils.getAssignedReferenceIndex(lastRead, header) == -1 ? null : new SimpleInterval(lastRead);
            writeRead(lastRead);

            if ( !tooManyReads )
                forMateMatching.clear();
            else
                purgeUnmodifiedMates();
        }

        // fix mates, as needed
        // Since setMateInfo can move reads, we potentially need to remove the mate, and requeue
        // it to ensure proper sorting
        if ( isMateFixableRead(newRead) ) {
            GATKReadHashObject mate = forMateMatching.get(newRead.getName());
            if ( mate != null ) {
                // 1. Frustratingly, Picard's setMateInfo() method unaligns (by setting the reference contig
                // to '*') read pairs when both of their flags have the unmapped bit set.  This is problematic
                // when trying to emit reads in coordinate order because all of a sudden we have reads in the
                // middle of the bam file that now belong at the end - and any mapped reads that get emitted
                // after them trigger an exception in the writer.  For our purposes, because we shouldn't be
                // moving read pairs when they are both unmapped anyways, we'll just not run fix mates on them.
                // 2. Furthermore, when reads get mapped to the junction of two chromosomes (e.g. MT since it
                // is actually circular DNA), their unmapped bit is set, but they are given legitimate coordinates.
                // The Picard code will come in and move the read all the way back to its mate (which can be
                // arbitrarily far away).  However, we do still want to move legitimately unmapped reads whose
                // mates are mapped, so the compromise will be that if the mate is still in the queue then we'll
                // move the read and otherwise we won't.
                boolean doNotFixMates = newRead.isUnmapped() && (mate.record.isUnmapped() || !waitingReads.contains(mate.record));
                if ( !doNotFixMates ) {

                    boolean reQueueMate = mate.record.isUnmapped() && ! newRead.isUnmapped();
                    if ( reQueueMate ) {
                        // the mate was unmapped, but newRead was mapped, so the mate may have been moved
                        // to be next-to newRead, so needs to be reinserted into the waitingReads queue
                        // note -- this must be called before the setMateInfo call below
                        if ( ! waitingReads.remove(mate.record) )
                            // we must have hit a region with too much depth and flushed the queue
                            reQueueMate = false;
                    }

                    // we've already seen our mate -- set the mate info and remove it from the map;
                    // add/update the mate cigar if appropriate
                    setMateInfo(mate.record, newRead, true);
                    if ( reQueueMate ) waitingReads.add(mate.record);
                }

                forMateMatching.remove(newRead.getName());
            } else if ( pairedReadIsMovable(newRead) ) {
                forMateMatching.put(newRead.getName(), new GATKReadHashObject(newRead, readWasModified));
            }
        }

        waitingReads.add(newRead);

        if ( ++counter % EMIT_FREQUENCY == 0 ) {
            while ( ! waitingReads.isEmpty() ) { // there's something in the queue
                GATKRead read = waitingReads.peek();

                if ( noReadCanMoveBefore(read.getAssignedStart(), newRead) &&
                        (!pairedReadIsMovable(read)                               // we won't try to move such a read
                                || noReadCanMoveBefore(read.getMateStart(), newRead ) ) ) { // we're already past where the mate started

                    // remove reads from the map that we have emitted -- useful for case where the mate never showed up
                    if ( !read.isSecondaryAlignment() )
                        forMateMatching.remove(read.getName());

                    if ( DEBUG )
                        logger.warn(String.format("EMIT!  At %d: read %s at %d with isize %d, mate start %d, op = %s",
                                newRead.getAssignedStart(), read.getName(), read.getAssignedStart(),
                                read.getFragmentLength(), read.getMateStart(), read.getAttributeAsString("OP")));
                    // emit to disk
                    writeRead(remove(waitingReads));
                } else {
                    if ( DEBUG )
                        logger.warn(String.format("At %d: read %s at %d with isize %d couldn't be emited, mate start %d",
                                newRead.getAssignedStart(), read.getName(), read.getAssignedStart(), read.getFragmentLength(), read.getMateStart()));
                    break;
                }
            }

            if ( DEBUG ) logger.warn(String.format("At %d: Done with emit cycle", newRead.getAssignedStart()));
        }
    }

    private void writeRead(GATKRead read) {
        try {
            if ( writer != null )
                writer.addRead(read);
        } catch (IllegalArgumentException e) {
            throw new UserException("If the maximum allowable reads in memory is too small, it may cause reads to be written out of order when trying to write the BAM; please see the --maxReadsInMemory argument for details.  " + e.getMessage(), e);
        }
    }

    /**
     * Is the given read one for which we can fix its mate?
     *
     * @param read  the read
     * @return true if we could fix its mate, false otherwise
     */
    protected boolean isMateFixableRead(final GATKRead read) {
        return read.isPaired() && !(read.isSecondaryAlignment() || read.isSupplementaryAlignment());
    }

    /**
     * @param read  the read
     * @return true if the read shouldn't be moved given the constraints of this SAMFileWriter
     */
    public boolean iSizeTooBigToMove(GATKRead read) {
        return iSizeTooBigToMove(read, maxInsertSizeForMovingReadPairs);               // we won't try to move such a read
    }

    public static boolean iSizeTooBigToMove(GATKRead read, int maxInsertSizeForMovingReadPairs) {
        return ( read.isPaired() && ! read.mateIsUnmapped() && !read.getAssignedContig().equals(read.getMateContig()) ) // maps to different chromosomes
                || Math.abs(read.getFragmentLength()) > maxInsertSizeForMovingReadPairs;     // we won't try to move such a read
    }

    private void purgeUnmodifiedMates() {
        HashMap<String, GATKReadHashObject> forMateMatchingCleaned = new HashMap<String, GATKReadHashObject>();
        for ( Map.Entry<String, GATKReadHashObject> entry : forMateMatching.entrySet() ) {
            if ( entry.getValue().wasModified )
                forMateMatchingCleaned.put(entry.getKey(), entry.getValue());
        }

        forMateMatching.clear(); // explicitly clear the memory
        forMateMatching = forMateMatchingCleaned;
    }

    private boolean pairedReadIsMovable(GATKRead read) {
        return read.isPaired()                                          // we're a paired read
                && (!read.isUnmapped() || !read.mateIsUnmapped())  // at least one read is mapped
                && !iSizeTooBigToMove(read);                                     // insert size isn't too big

    }

    public void close() {
        // write out all of the remaining reads
        while ( ! waitingReads.isEmpty() ) { // there's something in the queue
            writeRead(remove(waitingReads));
        }
    }


    // TODO: this is copied from HTSJDK: SamPairUtil.setMateInfo
    // TODO: maybe move to ReadUtils
    public static void setMateInfo(final GATKRead rec1, final GATKRead rec2, final boolean setMateCigar) {
        // If neither read is unmapped just set their mate info
        if (!rec1.isUnmapped() && !rec2.isUnmapped()) {
            rec1.setMatePosition(rec2);
            rec1.setMateIsReverseStrand(rec2.isReverseStrand());
            rec1.setAttribute(SAMTag.MQ.name(), rec2.getMappingQuality());

            rec2.setMatePosition(rec1);
            rec2.setMateIsReverseStrand(rec1.isReverseStrand());
            rec2.setAttribute(SAMTag.MQ.name(), rec1.getMappingQuality());

            if (setMateCigar) {
                rec1.setAttribute(SAMTag.MC.name(), TextCigarCodec.encode(rec2.getCigar()));
                rec2.setAttribute(SAMTag.MC.name(), TextCigarCodec.encode(rec1.getCigar()));
            }
            else {
                rec1.clearAttribute(SAMTag.MC.name());
                rec2.clearAttribute(SAMTag.MC.name());
            }
        }
        // Else if they're both unmapped set that straight
        else if (rec1.isUnmapped() && rec2.isUnmapped()) {
            rec1.setPosition(ReadConstants.UNSET_CONTIG, ReadConstants.UNSET_POSITION);
            rec1.setMatePosition(ReadConstants.UNSET_CONTIG, ReadConstants.UNSET_POSITION);
            rec1.setMateIsReverseStrand(rec2.isReverseStrand());
            rec1.setIsUnmapped();
            rec1.setMateIsUnmapped();
            rec1.clearAttribute(SAMTag.MQ.name());
            rec1.clearAttribute(SAMTag.MC.name());
            rec1.setFragmentLength(0);

            rec2.setPosition(ReadConstants.UNSET_CONTIG, ReadConstants.UNSET_POSITION);
            rec2.setMatePosition(ReadConstants.UNSET_CONTIG, ReadConstants.UNSET_POSITION);
            rec2.setMateIsReverseStrand(rec1.isReverseStrand());
            rec2.setIsUnmapped();
            rec2.setMateIsUnmapped();
            rec2.clearAttribute(SAMTag.MQ.name());
            rec2.clearAttribute(SAMTag.MC.name());
            rec2.setFragmentLength(0);
        }
        // And if only one is mapped copy it's coordinate information to the mate
        else {
            final GATKRead mapped   = rec1.isUnmapped() ? rec2 : rec1;
            final GATKRead unmapped = rec1.isUnmapped() ? rec1 : rec2;
            unmapped.setPosition(mapped);

            mapped.setMatePosition(unmapped);
            mapped.setMateIsReverseStrand(unmapped.isReverseStrand());
            mapped.setMateIsUnmapped();
            mapped.clearAttribute(SAMTag.MQ.name());
            mapped.clearAttribute(SAMTag.MC.name());
            mapped.setFragmentLength(0);

            unmapped.setMatePosition(mapped);
            unmapped.setMateIsReverseStrand(mapped.isReverseStrand());
            unmapped.setAttribute(SAMTag.MQ.name(), mapped.getMappingQuality());
            // For the unmapped read, set mateCigar to the mate's Cigar, since the mate must be mapped
            if (setMateCigar) unmapped.setAttribute(SAMTag.MC.name(), TextCigarCodec.encode(mapped.getCigar()));
            else unmapped.clearAttribute(SAMTag.MC.name());
            unmapped.setFragmentLength(0);
        }

        final int insertSize = computeInsertSize(rec1, rec2);
        rec1.setFragmentLength(insertSize);
        rec2.setFragmentLength(-insertSize);
    }

    // TODO: this is copied from HTSJDK: SamPairUtil.computeInsertSize
    // TODO: maybe move to ReadUtils
    public static int computeInsertSize(final GATKRead firstEnd, final GATKRead secondEnd) {
        if (firstEnd.isUnmapped() || secondEnd.isUnmapped()) {
            return 0;
        }
        if (!firstEnd.getAssignedContig().equals(secondEnd.getAssignedContig())) {
            return 0;
        }

        final int firstEnd5PrimePosition = firstEnd.isReverseStrand()? firstEnd.getEnd(): firstEnd.getAssignedStart();
        final int secondEnd5PrimePosition = secondEnd.isReverseStrand()? secondEnd.getEnd(): secondEnd.getAssignedStart();

        final int adjustment = (secondEnd5PrimePosition >= firstEnd5PrimePosition) ? +1 : -1;
        return secondEnd5PrimePosition - firstEnd5PrimePosition + adjustment;
    }

}
