package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.HasGenomeLocation;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Useful class for forwarding on locusContext data from this iterator
 */
public final class AlignmentContext implements HasGenomeLocation {
    private Locatable loc = null;
    private ReadPileup basePileup = null;
    private boolean hasPileupBeenDownsampled;

    /**
     * The number of bases we've skipped over in the reference since the last map invocation.
     * Only filled in by RodTraversals right now.  By default, nothing is being skipped, so skippedBases == 0.
     */
    private long skippedBases = 0;

    public AlignmentContext(Locatable loc, ReadPileup basePileup) {
        this(loc, basePileup, 0, false);
    }

    public AlignmentContext(Locatable loc, ReadPileup basePileup, boolean hasPileupBeenDownsampled) {
        this(loc, basePileup, 0, hasPileupBeenDownsampled);
    }

    public AlignmentContext(Locatable loc, ReadPileup basePileup, long skippedBases) {
        this(loc, basePileup, skippedBases, false);
    }

    public AlignmentContext(Locatable loc, ReadPileup basePileup, long skippedBases, boolean hasPileupBeenDownsampled) {
        if ( loc == null ) throw new GATKException("BUG: GenomeLoc in Alignment context is null");
        if ( basePileup == null ) throw new GATKException("BUG: ReadBackedPileup in Alignment context is null");
        if ( skippedBases < 0 ) throw new GATKException("BUG: skippedBases is -1 in Alignment context");

        this.loc = loc;
        this.basePileup = basePileup;
        this.skippedBases = skippedBases;
        this.hasPileupBeenDownsampled = hasPileupBeenDownsampled;
    }

    /**
     * Returns base pileup over the current genomic location. May return null if this context keeps only
     * extended event (indel) pileup.
     * @return
     */
    public ReadPileup getBasePileup() { return basePileup; }


    /**
     * Returns true if any reads have been filtered out of the pileup due to excess DoC.
     * @return True if reads have been filtered out.  False otherwise.
     */
    public boolean hasPileupBeenDownsampled() { return hasPileupBeenDownsampled; }

    /**
     * Are there any reads associated with this locus?
     *
     * @return
     */
    public boolean hasReads() {
        return basePileup != null && basePileup.size() > 0 ;
    }

    /**
     * How many reads cover this locus?
     * @return
     */
    public int size() {
        return basePileup.size();
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
    public Locatable getLocation() { return loc; }

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
