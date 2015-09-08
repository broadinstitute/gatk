package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.*;

import java.util.*;

/**
 * The class manages reads and splices and tries to apply overhang clipping when appropriate.
 * Important note: although for efficiency the manager does try to send reads to the underlying writer in coordinate
 * sorted order, it does NOT guarantee that it will do so in every case!  So unless there's a good reason not to,
 * methods that instantiate this manager should pass in a writer that does not assume the reads are pre-sorted.
 */
public class OverhangFixingManager {

    protected static final Logger logger = LogManager.getLogger(OverhangFixingManager.class);
    private static final boolean DEBUG = false;

    // how many reads should we store in memory before flushing the queue?
    private final int MAX_RECORDS_IN_MEMORY;

    // how many mismatches do we tolerate in the overhangs?
    private final int MAX_MISMATCHES_IN_OVERHANG;

    // how many bases do we tolerate in the overhang before deciding not to clip?
    private final int MAX_BASES_IN_OVERHANG;

    // should we not bother fixing overhangs?
    private final boolean doNotFixOverhangs;

    // header for the reads
    private final SAMFileHeader header;

    // where we ultimately write out our records
    private final SAMFileGATKReadWriter writer;

    // fasta reference reader to check overhanging edges in the exome reference sequence
    private final IndexedFastaSequenceFile referenceReader;

    // the genome loc parser
    private final GenomeLocParser genomeLocParser;

    // the read cache
    private static final int initialCapacity = 5000;
    private final PriorityQueue<SplitRead> waitingReads;

    // the set of current splices to use
    private final Set<Splice> splices = new TreeSet<>(new SpliceComparator());

    protected static final int MAX_SPLICES_TO_KEEP = 1000;


    /**
     *
     * @param header                   header for the reads
     * @param writer                   actual writer
     * @param genomeLocParser          the GenomeLocParser object
     * @param referenceReader          the reference reader
     * @param maxRecordsInMemory       max records to keep in memory
     * @param maxMismatchesInOverhangs max number of mismatches permitted in the overhangs before requiring clipping
     * @param maxBasesInOverhangs      max number of bases permitted in the overhangs before deciding not to clip
     * @param doNotFixOverhangs        if true, don't clip overhangs at all
     */
    public OverhangFixingManager(final SAMFileHeader header,
                                 final SAMFileGATKReadWriter writer,
                                 final GenomeLocParser genomeLocParser,
                                 final IndexedFastaSequenceFile referenceReader,
                                 final int maxRecordsInMemory,
                                 final int maxMismatchesInOverhangs,
                                 final int maxBasesInOverhangs,
                                 final boolean doNotFixOverhangs) {
        this.header = header;
        this.writer = writer;
        this.genomeLocParser = genomeLocParser;
        this.referenceReader = referenceReader;
        this.MAX_RECORDS_IN_MEMORY = maxRecordsInMemory;
        this.MAX_MISMATCHES_IN_OVERHANG = maxMismatchesInOverhangs;
        this.MAX_BASES_IN_OVERHANG = maxBasesInOverhangs;
        this.doNotFixOverhangs = doNotFixOverhangs;
        this.waitingReads = new PriorityQueue<>(initialCapacity, new SplitReadComparator());
    }

    public final int getNReadsInQueue() { return waitingReads.size(); }

    /**
     * For testing purposes only
     *
     * @return the list of reads currently in the queue
     */
    public List<SplitRead> getReadsInQueueForTesting() {
        return new ArrayList<>(waitingReads);
    }

    /**
     * For testing purposes only
     *
     * @return the list of splices currently in the queue
     */
    public List<Splice> getSplicesForTesting() {
        return new ArrayList<>(splices);
    }

    /**
     * Add a new observed split to the list to use
     *
     * @param contig  the contig
     * @param start   the start of the split, inclusive
     * @param end     the end of the split, inclusive
     */
    public void addSplicePosition(final String contig, final int start, final int end) {
        if ( doNotFixOverhangs ) {
            return;
        }

        // is this a new splice?  if not, we are done
        final Splice splice = new Splice(contig, start, end);
        if ( splices.contains(splice) ) {
            return;
        }

        // initialize it with the reference context
        // we don't want to do this until we know for sure that it's a new splice position
        splice.initialize(referenceReader);

        // clear the set of old split positions seen if we hit a new contig
        final boolean sameContig = splices.isEmpty() || splices.iterator().next().loc.getContig().equals(contig);
        if ( !sameContig ) {
            splices.clear();
        }

        // run this position against the existing reads
        for ( final SplitRead read : waitingReads ) {
            fixSplit(read, splice);
        }

        splices.add(splice);

        if ( splices.size() > MAX_SPLICES_TO_KEEP ) {
            cleanSplices();
        }
    }

    /**
     * Add a read to the manager
     *
     * @param read  the read to add
     */
    public void addRead(final GATKRead read) {
        if ( read == null ) {
            throw new IllegalArgumentException("read added to manager is null, which is not allowed");
        }

        // if the new read is on a different contig or we have too many reads, then we need to flush the queue and clear the map
        final boolean tooManyReads = getNReadsInQueue() >= MAX_RECORDS_IN_MEMORY;
        final boolean encounteredNewContig = getNReadsInQueue() > 0 && !waitingReads.peek().read.getContig().equals(read.getContig());

        if ( tooManyReads || encounteredNewContig ) {
            if ( DEBUG ) {
                logger.warn("Flushing queue on " + (tooManyReads ? "too many reads" : ("move to new contig: " + read.getContig() + " from " + waitingReads.peek().read.getContig())) + " at " + read.getStart());
            }

            final int targetQueueSize = encounteredNewContig ? 0 : MAX_RECORDS_IN_MEMORY / 2;

            // write the required number of waiting reads to disk
            while ( getNReadsInQueue() > targetQueueSize ) {
                writer.addRead(waitingReads.poll().read);
            }
        }

        final SplitRead splitRead = new SplitRead(read);

        // fix overhangs, as needed
        for ( final Splice splice : splices) {
            fixSplit(splitRead, splice);
        }

        // add the new read to the queue
        waitingReads.add(splitRead);
    }

    /**
     * Clean up the list of splices
     */
    private void cleanSplices() {
        final int targetQueueSize = splices.size() / 2;
        final Iterator<Splice> iter = splices.iterator();
        for ( int i = 0; i < targetQueueSize; i++ ) {
            iter.next();
            iter.remove();
        }
    }

    /**
     * Try to fix the given read using the given split
     *
     * @param read        the read to fix
     * @param splice      the split (bad region to clip out)
     */
    private void fixSplit(final SplitRead read, final Splice splice) {
        // if the read doesn't even overlap the split position then we can just exit
        if ( read.loc == null || !splice.loc.overlapsP(read.loc) ) {
            return;
        }

        if ( isLeftOverhang(read.loc, splice.loc) ) {
            final int overhang = splice.loc.getStop() - read.loc.getStart() + 1;
            if ( overhangingBasesMismatch(read.read.getBases(), 0, splice.reference, splice.reference.length - overhang, overhang) ) {
                final GATKRead clippedRead = ReadClipper.hardClipByReadCoordinates(read.read, 0, overhang - 1);
                read.setRead(clippedRead);
            }
        }
        else if ( isRightOverhang(read.loc, splice.loc) ) {
            final int overhang = read.loc.getStop() - splice.loc.getStart() + 1;
            if ( overhangingBasesMismatch(read.read.getBases(), read.read.getLength() - overhang, splice.reference, 0, overhang) ) {
                final GATKRead clippedRead = ReadClipper.hardClipByReadCoordinates(read.read, read.read.getLength() - overhang, read.read.getLength() - 1);
                read.setRead(clippedRead);
            }
        }
    }

    /**
     * Is this a proper overhang on the left side of the read?
     *
     * @param readLoc    the read's loc
     * @param spliceLoc   the split's loc
     * @return true if it's a left side overhang
     */
    protected static boolean isLeftOverhang(final GenomeLoc readLoc, final GenomeLoc spliceLoc) {
        return readLoc.getStart() <= spliceLoc.getStop() && readLoc.getStart() > spliceLoc.getStart() && readLoc.getStop() > spliceLoc.getStop();
    }

    /**
     * Is this a proper overhang on the right side of the read?
     *
     * @param readLoc    the read's loc
     * @param spliceLoc   the split's loc
     * @return true if it's a right side overhang
     */
    protected static boolean isRightOverhang(final GenomeLoc readLoc, final GenomeLoc spliceLoc) {
        return readLoc.getStop() >= spliceLoc.getStart() && readLoc.getStop() < spliceLoc.getStop() && readLoc.getStart() < spliceLoc.getStart();
    }

    /**
     * Are there too many mismatches to the reference among the overhanging bases?
     *
     * @param read                  the read bases
     * @param readStartIndex        where to start on the read
     * @param reference             the reference bases
     * @param referenceStartIndex   where to start on the reference
     * @param spanToTest            how many bases to test
     * @return true if too many overhanging bases mismatch, false otherwise
     */
    protected boolean overhangingBasesMismatch(final byte[] read,
                                               final int readStartIndex,
                                               final byte[] reference,
                                               final int referenceStartIndex,
                                               final int spanToTest) {
        // don't process too small a span, too large a span, or a span that is most of a read
        if ( spanToTest < 1 || spanToTest > MAX_BASES_IN_OVERHANG || spanToTest > read.length / 2 ) {
            return false;
        }

        int numMismatchesSeen = 0;
        for ( int i = 0; i < spanToTest; i++ ) {
            if ( read[readStartIndex + i] != reference[referenceStartIndex + i] ) {
                if ( ++numMismatchesSeen > MAX_MISMATCHES_IN_OVERHANG ) {
                    return true;
                }
            }
        }

        // we can still mismatch overall if at least half of the bases mismatch
        return numMismatchesSeen >= ((spanToTest+1)/2);
    }

    /**
     * Close out the manager stream by clearing the read cache
     */
    public void close() {
        // write out all of the remaining reads
        while ( ! waitingReads.isEmpty() ) {
            writer.addRead(waitingReads.poll().read);
        }
    }

    // class to represent the reads with their soft-clip-included GenomeLocs
    public final class SplitRead {

        public GATKRead read;
        public GenomeLoc loc;

        public SplitRead(final GATKRead read) {
            setRead(read);
        }

        public void setRead(final GATKRead read) {
            if ( ! read.isEmpty() ) {
                this.read = read;
                if ( ! read.isUnmapped() ) {
                    loc = genomeLocParser.createGenomeLoc(read.getContig(), ReadUtils.getSoftStart(read), ReadUtils.getSoftEnd(read));
                }
            }
        }
    }

    // class to represent the comparator for the split reads
    private final class SplitReadComparator implements Comparator<SplitRead> {

        private final ReadCoordinateComparator readComparator;

        public SplitReadComparator() {
            readComparator = new ReadCoordinateComparator(header);
        }

        public int compare(final SplitRead read1, final SplitRead read2) {
            return readComparator.compare(read1.read, read2.read);
        }
    }

    // class to represent the split positions
    protected final class Splice {

        public final GenomeLoc loc;
        public byte[] reference;

        public Splice(final String contig, final int start, final int end) {
            loc = genomeLocParser.createGenomeLoc(contig, start, end);
        }

        public void initialize(final IndexedFastaSequenceFile referenceReader) {
            reference = referenceReader.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
        }

        @Override
        public boolean equals(final Object other) {
            return other != null && (other instanceof Splice) && this.loc.equals(((Splice)other).loc);
        }

        @Override
        public int hashCode() {
            return loc.hashCode();
        }
    }

    // class to represent the comparator for the split reads
    private final class SpliceComparator implements Comparator<Splice> {

        public int compare(final Splice position1, final Splice position2) {
            return position1.loc.compareTo(position2.loc);
        }
    }
}
