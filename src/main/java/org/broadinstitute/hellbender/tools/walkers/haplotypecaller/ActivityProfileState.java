package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * The state of an active region walker's isActive call at a specific locus in the genome
 */
public final class ActivityProfileState {
    private final GenomeLoc loc;
    private final double activeProb;
    private final Type resultState;
    private final Number resultValue;

    public double isActiveProb() {
        return activeProb;
    }

    public Type getResultState() {
        return resultState;
    }

    public Number getResultValue() {
        return resultValue;
    }

    public enum Type {
        NONE,
        HIGH_QUALITY_SOFT_CLIPS
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of activeProb
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param activeProb the probability of being active (between 0 and 1)
     */
    public ActivityProfileState(final GenomeLoc loc, final double activeProb) {
        this(loc, activeProb, Type.NONE, null);
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of activeProb that maintains some
     * information about the result state and value
     *
     * The only state value in use is HIGH_QUALITY_SOFT_CLIPS, and here the value is interpreted as the number
     * of bp affected by the soft clips.
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param activeProb the probability of being active (between 0 and 1)
     */
    public ActivityProfileState(final GenomeLoc loc, final double activeProb, final Type resultState, final Number resultValue) {
        // make sure the location of that activity profile is 1
        if ( loc.size() != 1 ) {
            throw new IllegalArgumentException("Location for an ActivityProfileState must have to size 1 bp but saw " + loc);
        }
        if ( resultValue != null && resultValue.doubleValue() < 0 ) {
            throw new IllegalArgumentException("Result value isn't null and its < 0, which is illegal: " + resultValue);
        }

        this.loc = loc;
        this.activeProb = activeProb;
        this.resultState = resultState;
        this.resultValue = resultValue;
    }

    /**
     * The offset of state w.r.t. our current region's start location
     * @param regionStartLoc the start of the region, as a genome loc
     * @return the position of this profile relative to the start of this region
     */
    public int getOffset(final Locatable regionStartLoc) {
        Utils.nonNull(regionStartLoc);
        return getLoc().getStart() - regionStartLoc.getStart();
    }

    /**
     * Get the genome loc associated with the ActivityProfileState
     * @return the location of this result
     */
    public GenomeLoc getLoc() {
        return loc;
    }

    @Override
    public String toString() {
        return "ActivityProfileState{" +
                "loc=" + loc +
                ", activeProb=" + activeProb +
                ", resultState=" + resultState +
                ", resultValue=" + resultValue +
                '}';
    }
}
