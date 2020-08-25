package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.List;

/**
 * @author aaron
 *
 * allows query calls to the artificial read iterator, which allows you
 * to test out classes that use specific itervals.  The reads returned will
 * all lie in order in the specified interval.
 */
public final class ArtificialReadQueryIterator extends ArtificialReadIterator {

    // get the next position
    protected int finalPos = 0;
    protected int startPos = 0;
    protected int contigIndex = -1;
    protected boolean overlapping = false;
    protected int startingChr = 0;
    protected boolean seeked = false;

    /**
     * create the fake iterator, given the mapping of chromosomes and read counts
     *
     * @param startingChr the starting chromosome
     * @param endingChr   the ending chromosome
     * @param readCount   the number of reads in each chromosome
     * @param header      the associated header
     */
    ArtificialReadQueryIterator( int startingChr, int endingChr, int readCount, int unmappedReadCount, SAMFileHeader header ) {
        super(startingChr, endingChr, readCount, unmappedReadCount, header);
        this.startingChr = startingChr;
    }

    @Override
    protected void reset() {
        this.startPos = 0;
        this.finalPos = 0;
        this.contigIndex = -1;
        // Doesn't make sense to reset the overlapping flag, because we rely on its state later on.
        // TODO: Make this a bit more direct.
        //overlapping = false;
        this.startingChr = 0;
        this.seeked = false;
        super.reset();
    }

    /**
     * query containing - get reads contained by the specified interval
     *
     * @param contig the contig index string
     * @param start  the start position
     * @param stop   the stop position
     */
    public void queryContained( String contig, int start, int stop ) {
        this.overlapping = false;
        initialize(contig, start, stop);
    }

    /**
     * query containing - get reads contained by the specified interval
     *
     * @param contig the contig index string
     * @param start  the start position
     * @param stop   the stop position
     */
    public void queryOverlapping( String contig, int start, int stop ) {
        this.overlapping = true;
        initialize(contig, start, stop);
    }

    public void query( String contig, int start, int stop, boolean contained ) {
        if (contained)
            queryContained(contig, start, stop);
        else
            queryOverlapping(contig, start, stop);
    }


    /**
     * initialize the query iterator
     *
     * @param contig the contig
     * @param start  the start position
     * @param stop   the stop postition
     */
    private void initialize( String contig, int start, int stop ) {
        // throw away data from the previous invocation, if one exists.
        ensureUntouched();
        reset();

        finalPos = stop;
        startPos = start;
        if (finalPos < 0) {
            finalPos = Integer.MAX_VALUE;
        }
        // sanity check that we have the contig
        contigIndex = -1;
        List<SAMSequenceRecord> list = header.getSequenceDictionary().getSequences();
        for (SAMSequenceRecord rec : list) {
            if (rec.getSequenceName().equals(contig)) {
                contigIndex = rec.getSequenceIndex();
            }
        }
        if (contigIndex < 0) { throw new IllegalArgumentException("ArtificialContig" + contig + " doesn't exist"); }
        while (super.hasNext() && ReadUtils.getReferenceIndex(this.peek(), header) < contigIndex) {
            super.next();
        }
        if (!super.hasNext()) {
            throw new GATKException("Unable to find the target chromosome");
        }
        while (super.hasNext() && this.peek().getStart() < start) {
            super.next();
        }
        // sanity check that we have an actual matching read next
        GATKRead rec = this.peek();
        if (!matches(rec)) {
            throw new GATKException("The next read doesn't match");
        }
        // set the seeked variable to true
        seeked = true;
    }

    /**
     * given a read and the query type, check if it matches our regions
     *
     * @param rec the read
     *
     * @return true if it belongs in our region
     */
    public boolean matches( GATKRead rec ) {
        final int recReferenceIndex = ReadUtils.getReferenceIndex(rec, header);

        if (recReferenceIndex != this.contigIndex) {
            return false;
        }
        // if we have an unmapped read, matching the contig is good enough for us
        if (recReferenceIndex < 0) {
            return true;
        }

        if (!overlapping) {
            // if the start or the end are somewhere within our range
            if (( rec.getStart() >= startPos && rec.getEnd() <= finalPos )) {
                return true;
            }
        } else {
            if (( rec.getStart() <= finalPos && rec.getStart() >= startPos ) ||
                    ( rec.getEnd() <= finalPos && rec.getEnd() >= startPos )) {
                return true;
            }
        }
        return false;
    }


    /**
     * override the hasNext, to incorportate our limiting factor
     *
     * @return
     */
    @Override
    public boolean hasNext() {
        boolean res = super.hasNext();
        if (!seeked) {
            return res;
        }
        return res && matches(this.next);
    }

    /** make sure we haven't been used as an iterator yet; this is to miror the MergingSamIterator2 action. */
    public void ensureUntouched() {
        if (open) {
            throw new UnsupportedOperationException("We've already been used as an iterator; you can't query after that");
        }
    }
}
