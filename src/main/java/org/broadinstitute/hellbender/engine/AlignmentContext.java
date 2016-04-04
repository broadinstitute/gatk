package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.HasGenomeLocation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

/**
 * Bundles together a pileup and a location.
 */
public final class AlignmentContext implements Locatable, HasGenomeLocation {
    private final Locatable loc;
    private final ReadPileup basePileup;

    public AlignmentContext(final Locatable loc, final ReadPileup basePileup) {
        Utils.nonNull(loc, "BUG: GenomeLoc in Alignment context is null");
        Utils.nonNull(basePileup, "BUG: ReadBackedPileup in Alignment context is null");

        this.loc = loc;
        this.basePileup = basePileup;
    }

	@Override
	public int getStart() {	return getLocation().getStart(); }

	@Override
	public int getEnd() { return getLocation().getEnd(); }

    /**
     * How many reads cover this locus?
     * @return
     */
    public int size() {
        return basePileup.size();
    }

    public String getContig() { return getLocation().getContig(); }

    public long getPosition() { return getLocation().getStart(); }

    public Locatable getLocation() { return loc; }

    public ReadPileup getBasePileup(){
        return basePileup;
    }

}
