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

package org.broadinstitute.gatk.utils.contexts;

import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.HasGenomeLocation;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.List;

/**
 * Useful class for forwarding on locusContext data from this iterator
 * 
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:01:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class AlignmentContext implements HasGenomeLocation {
    protected GenomeLoc loc = null;
    protected ReadBackedPileup basePileup = null;
    protected boolean hasPileupBeenDownsampled;

    /**
     * The number of bases we've skipped over in the reference since the last map invocation.
     * Only filled in by RodTraversals right now.  By default, nothing is being skipped, so skippedBases == 0.
     */
    private long skippedBases = 0;

    public AlignmentContext(GenomeLoc loc, ReadBackedPileup basePileup) {
        this(loc, basePileup, 0, false);
    }

    public AlignmentContext(GenomeLoc loc, ReadBackedPileup basePileup, boolean hasPileupBeenDownsampled) {
        this(loc, basePileup, 0, hasPileupBeenDownsampled);
    }

    public AlignmentContext(GenomeLoc loc, ReadBackedPileup basePileup, long skippedBases) {
        this(loc, basePileup, skippedBases, false);
    }

    public AlignmentContext(GenomeLoc loc, ReadBackedPileup basePileup, long skippedBases,boolean hasPileupBeenDownsampled ) {
        if ( loc == null ) throw new ReviewedGATKException("BUG: GenomeLoc in Alignment context is null");
        if ( basePileup == null ) throw new ReviewedGATKException("BUG: ReadBackedPileup in Alignment context is null");
        if ( skippedBases < 0 ) throw new ReviewedGATKException("BUG: skippedBases is -1 in Alignment context");

        this.loc = loc;
        this.basePileup = basePileup;
        this.skippedBases = skippedBases;
        this.hasPileupBeenDownsampled = hasPileupBeenDownsampled;
    }

    /** Returns base pileup over the current genomic location. Deprectated. Use getBasePileup() to make your intentions
     * clear.
     * @return
     */
    @Deprecated
    public ReadBackedPileup getPileup() { return basePileup; }

    /** Returns base pileup over the current genomic location. May return null if this context keeps only
     * extended event (indel) pileup.
     * @return
     */
    public ReadBackedPileup getBasePileup() {
        return basePileup;
    }

    /**
     * Returns true if any reads have been filtered out of the pileup due to excess DoC.
     * @return True if reads have been filtered out.  False otherwise.
     */
    public boolean hasPileupBeenDownsampled() { return hasPileupBeenDownsampled; }

    /**
     * get all of the reads within this context
     * 
     * @return
     */
    @Deprecated
    //todo: unsafe and tailored for current usage only; both pileups can be null or worse, bot can be not null in theory
    public List<GATKSAMRecord> getReads() { return ( basePileup.getReads() ); }

    /**
     * Are there any reads associated with this locus?
     *
     * @return
     */
    public boolean hasReads() {
        return basePileup != null && basePileup.getNumberOfElements() > 0 ;
    }

    /**
     * How many reads cover this locus?
     * @return
     */
    public int size() {
        return basePileup.getNumberOfElements();
    }

    /**
     * get a list of the equivalent positions within in the reads at Pos
     *
     * @return
     */
    @Deprecated
    public List<Integer> getOffsets() {
        return basePileup.getOffsets();
    }

    public String getContig() { return getLocation().getContig(); }
    public long getPosition() { return getLocation().getStart(); }
    public GenomeLoc getLocation() { return loc; }

    public void downsampleToCoverage(int coverage) {
        basePileup = basePileup.getDownsampledPileup(coverage);
        hasPileupBeenDownsampled = true;
    }

    /**
     * Returns the number of bases we've skipped over in the reference since the last map invocation.
     * Only filled in by RodTraversals right now.  A value of 0 indicates that no bases were skipped.
     *
     * @return the number of skipped bases
     */
    public long getSkippedBases() {
        return skippedBases;
    }
}
