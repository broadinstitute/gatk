package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.EnumMap;
import java.util.Map;
import java.util.Set;

/**
 * Captures the probability that a specific locus in the genome represents an "active" site containing
 * real variation.
 */
public final class ActivityProfileState {
    private final SimpleInterval loc;
    private double activeProb;
    private final Type resultState;
    private final Number resultValue;
    // the kinds of alt allele that the evidence at this locus supports, each with its log odds, if the caller reports them,
    // otherwise empty
    private final Map<AltAlleleType, Double> altAlleleLogOdds;
    // TODO this is used to store the original isActive() prob (before smoothing) fast genotyping likelihood determination score which is used as a heuristic in DRAGEN-GATK for pileup-calling
    private double originalActiveProb = 0.0;

    // When range-checking probabilities, we allow this much tolerance.
    private static final double PROBABILITY_TOLERANCE = 0.01;

    public double isActiveProb() {
        return activeProb;
    }

    /**
     * Set the probability that this site is active.
     *
     * Probabilities should be between 0.0 and 1.0, however this is not currently enforced
     * because the {@link BandPassActivityProfile} can sometimes generate probabilities that
     * slightly exceed 1.0 when moving probability mass around. We intend to fix this by
     * capping at 1.0, but first we must evaluate the effects of capping on the HaplotypeCaller.
     *
     * @param activeProb probability (should be between 0.0 and 1.0) that the site is active
     */
    public void setIsActiveProb( final double activeProb ) {
        this.activeProb = activeProb;
    }

    /**
     * @return The type of the value returned by {@link #getResultValue}
     */
    public Type getResultState() {
        return resultState;
    }

    /**
     * @return Numeric value associated with {@link #getResultState}. If {@link #getResultState} is HIGH_QUALITY_SOFT_CLIPS,
     *         this is the number of bp affected by the soft clips
     */
    public Number getResultValue() {
        return resultValue;
    }

    public double getOriginalActiveProb() {
        return originalActiveProb;
    }

    public void setOriginalActiveProb(double originalActiveProb) {
        this.originalActiveProb = originalActiveProb;
    }

    /**
     * The type of the value returned by {@link #getResultValue}
     */
    public enum Type {
        NONE,
        HIGH_QUALITY_SOFT_CLIPS,
        GERMLINE,
        SOMATIC,
        ARTIFACT
    }

    /**
     * Kinds of alt allele that the evidence at a locus can support: a substitution of one of the four bases, an insertion
     * or deletion immediately after the locus (or a deletion spanning it), or anything else, such as a non-ACGT base or
     * a high-quality soft clip
     */
    public enum AltAlleleType {
        A, C, G, T, INSERTION, DELETION, OTHER;

        /**
         * @return the substituted base if this is a substitution of A, C, G, or T, otherwise 0
         */
        public byte substitutedBase() {
            return ordinal() < 4 ? (byte) name().charAt(0) : 0;
        }
    }

    /**
     * @return the kinds of alt allele that the evidence at this locus supports, if the caller reported them, otherwise empty
     */
    public Set<AltAlleleType> getAltAlleleTypes() {
        return altAlleleLogOdds.keySet();
    }

    /**
     * @return the kinds of alt allele that the evidence at this locus supports, each with the caller's log odds (natural
     * log) that it is present, if the caller reported them, otherwise empty
     */
    public Map<AltAlleleType, Double> getAltAlleleLogOdds() {
        return altAlleleLogOdds;
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of activeProb
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param activeProb the probability of being active (between 0 and 1)
     */
    public ActivityProfileState(final SimpleInterval loc, final double activeProb) {
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
    public ActivityProfileState(final SimpleInterval loc, final double activeProb, final Type resultState, final Number resultValue) {
        this(loc, activeProb, resultState, resultValue, Collections.emptyMap());
    }

    /**
     * As above, also recording the kinds of alt allele that the evidence at this locus supports
     *
     * @param altAlleleLogOdds the kinds of alt allele supported at this locus, which may be none, each with its log odds
     */
    public ActivityProfileState(final SimpleInterval loc, final double activeProb, final Type resultState, final Number resultValue,
                                final Map<AltAlleleType, Double> altAlleleLogOdds) {
        Utils.nonNull(altAlleleLogOdds);
        if ( loc.size() != 1 ) {
            throw new IllegalArgumentException("Location for an ActivityProfileState must have to size 1 bp but saw " + loc);
        }
        if ( resultValue != null && resultValue.doubleValue() < 0 ) {
            throw new IllegalArgumentException("Result value isn't null and its < 0, which is illegal: " + resultValue);
        }

        this.loc = loc;
        setIsActiveProb(activeProb);
        this.resultState = resultState;
        this.resultValue = resultValue;
        this.altAlleleLogOdds = altAlleleLogOdds.isEmpty() ? Collections.emptyMap() : Collections.unmodifiableMap(new EnumMap<>(altAlleleLogOdds));
    }

    /**
     * The offset of state w.r.t. our current region's start location
     * @param regionStartLoc the start of the region, as a Locatable
     * @return the position of this profile relative to the start of this region
     */
    public int getOffset(final Locatable regionStartLoc) {
        Utils.nonNull(regionStartLoc);
        return getLoc().getStart() - regionStartLoc.getStart();
    }

    /**
     * Get the locus associated with the ActivityProfileState
     * @return the locus associated with the ActivityProfileState as a SimpleInterval
     */
    public SimpleInterval getLoc() {
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
