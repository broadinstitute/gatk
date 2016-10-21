package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import java.util.Iterator;


/** this fake iterator allows us to look at how specific piles of reads are handled */
public class ArtificialReadIterator implements Iterator<GATKRead>, Iterable<GATKRead> {


    protected int currentChromo = 0;
    protected int currentRead = 1;
    protected int totalReadCount = 0;
    protected int unmappedRemaining = 0;
    protected boolean done = false;
    // the next record
    protected GATKRead next = null;
    protected SAMFileHeader header = null;

    // the passed in parameters
    protected final int sChr;
    protected final int eChromosomeCount;
    protected final int rCount;
    protected final int unmappedReadCount;

    // let us know to make a read, we need this to help out the fake sam query iterator
    private boolean initialized = false;

    /**
     * Is this iterator currently open or closed?  Closed iterators can be reused.
     */
    protected boolean open = false;

    /**
     * create the fake iterator, given the mapping of chromosomes and read counts
     *
     * @param startingChr the starting chromosome
     * @param endingChr   the ending chromosome
     * @param readCount   the number of reads in each chromosome
     * @param header      the associated header
     */
    ArtificialReadIterator( int startingChr, int endingChr, int readCount, SAMFileHeader header ) {
        sChr = startingChr;
        eChromosomeCount = (endingChr - startingChr) + 1;
        rCount = readCount;
        this.header = header;
        unmappedReadCount = 0;
        reset();
    }

    protected void reset() {
        this.currentChromo = 0;
        this.currentRead = 1;
        this.totalReadCount = 0;
        this.done = false;
        this.next = null;
        this.initialized = false;
        this.unmappedRemaining = unmappedReadCount;
    }

    /**
     * create the fake iterator, given the mapping of chromosomes and read counts
     *
     * @param startingChr the starting chromosome
     * @param endingChr   the ending chromosome
     * @param readCount   the number of reads in each chromosome
     * @param header      the associated header
     */
    ArtificialReadIterator( int startingChr, int endingChr, int readCount, int unmappedReadCount, SAMFileHeader header ) {
        sChr = startingChr;
        eChromosomeCount = (endingChr - startingChr) + 1;
        rCount = readCount;
        this.header = header;
        this.currentChromo = 0;
        this.unmappedReadCount = unmappedReadCount;
        reset();
    }

    public void close() {
        open = false;
    }

    public SAMFileHeader getHeader() {
        return header;
    }

    @Override
    public boolean hasNext() {
        open = true;

        if (!initialized){
            initialized = true;
            createNextRead();
        }
        return this.next != null;
    }

    protected boolean createNextRead() {
        if (currentRead > rCount) {
            currentChromo++;
            currentRead = 1;
        }
        // check for end condition, have we finished the chromosome listing, and have no unmapped reads
        if (currentChromo >= eChromosomeCount) {
            if (unmappedRemaining < 1) {
                this.next = null;
                return false;
            } else {
                ++totalReadCount;
                this.next = ArtificialReadUtils.createArtificialRead(this.header,
                        String.valueOf(totalReadCount),
                        SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX,
                        SAMRecord.NO_ALIGNMENT_START,
                        50);
                --unmappedRemaining;
                return true;
            }
        }
        ++totalReadCount;
        this.next = ArtificialReadUtils.createArtificialRead(this.header, String.valueOf(totalReadCount), currentChromo, currentRead, 50);
        ++currentRead;
        return true;
    }


    @Override
    public GATKRead next() {
        open = true;

        GATKRead ret = next;
        createNextRead();
        return ret;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("You've tried to remove on a GATKSAMIterator (unsupported), not to mention that this is a fake iterator.");
    }

    /**
     * return this iterator, for the iterable interface
     */
    @Override
    public Iterator<GATKRead> iterator() {
        return this;
    }

    /**
     * some instrumentation methods
     */
    public int readsTaken() {
        return totalReadCount;
    }

    /**
     * peek at the next sam record
     */
    public GATKRead peek() {
        return this.next;
    }
}
