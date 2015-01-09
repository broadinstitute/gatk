/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.hellbender.utils.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.iterators.GATKSAMIterator;

import java.util.Iterator;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/** this fake iterator allows us to look at how specific piles of reads are handled */
public class ArtificialSAMIterator implements GATKSAMIterator {


    protected int currentChromo = 0;
    protected int currentRead = 1;
    protected int totalReadCount = 0;
    protected int unmappedRemaining = 0;
    protected boolean done = false;
    // the next record
    protected SAMRecord next = null;
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
    ArtificialSAMIterator( int startingChr, int endingChr, int readCount, SAMFileHeader header ) {
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
    ArtificialSAMIterator( int startingChr, int endingChr, int readCount, int unmappedReadCount, SAMFileHeader header ) {
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

    public boolean hasNext() {
        open = true;

        if (!initialized){
            initialized = true;
            createNextRead();
        }
        if (this.next != null) {
            return true;
        }
        return false;
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
                this.next = ArtificialSAMUtils.createArtificialRead(this.header,
                        String.valueOf(totalReadCount),
                        SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX,
                        SAMRecord.NO_ALIGNMENT_START,
                        50);
                --unmappedRemaining;
                return true;
            }
        }
        ++totalReadCount;
        this.next = ArtificialSAMUtils.createArtificialRead(this.header, String.valueOf(totalReadCount), currentChromo, currentRead, 50);
        ++currentRead;
        return true;
    }


    public SAMRecord next() {
        open = true;

        SAMRecord ret = next;
        createNextRead();
        return ret;
    }

    public void remove() {
        throw new UnsupportedOperationException("You've tried to remove on a GATKSAMIterator (unsupported), not to mention that this is a fake iterator.");
    }

    /**
     * return this iterator, for the iterable interface
     */
    public Iterator<SAMRecord> iterator() {
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
    public SAMRecord peek() {
        return this.next;
    }
}
