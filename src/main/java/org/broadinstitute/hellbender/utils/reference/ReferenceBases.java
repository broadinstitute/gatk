package org.broadinstitute.hellbender.utils.reference;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.Arrays;

/**
 * ReferenceBases stores the bases of the reference genome for a particular interval.
 * This class requires the bases to be encoded at 8 bits per base.
 */
public final class ReferenceBases implements Serializable {
    private static final long serialVersionUID = 1L;

    private final byte[] bases;
    private final SimpleInterval interval;

    public ReferenceBases( final byte[] bases, final SimpleInterval interval ) {
        Utils.nonNull(bases);
        Utils.nonNull(interval);
        if (interval.size() != bases.length) {
            throw new IllegalArgumentException(
                    "interval must have same length as bases, " + interval + " " + interval.size() + "," + bases.length);
        }
        this.bases = bases;
        this.interval = interval;
    }

    @Override
    public String toString() {
        return "ReferenceBases{" +
                "bases=" + Arrays.toString(bases) +
                ", interval=" + interval +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ReferenceBases that = (ReferenceBases) o;

        if (!Arrays.equals(getBases(), that.getBases())) return false;
        return getInterval().equals(that.getInterval());

    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(getBases());
        result = 31 * result + getInterval().hashCode();
        return result;
    }

    /**
     * Returns the bases.
     * @return never {@code null}, this is a direct reference to this object bases, so the caller must refrain from modifying them
     *  as it will alter the state of this object.
     */
    public byte[] getBases() {
        return bases;
    }

    /**
     * Returns a copy of a subrange of the bases given the start and end positions (1-based closed).
     *
     * <p>
     *     This is a more efficient short-hand for {@code this.getSubset(new SimpleInterval(this.getContig(), start, end)).getBases() }.
     * </p>
     * @param start the sub-interval start.
     * @param end the sub-interval end.
     * @IllegalArgumentException if the input start-end interval is not completelly contained in this reference-bases own interval or it
     * is an empty interval (start > end).
     * @return the returned array is a new byte array so the caller is free to modify it.
     */
    public byte[] getSubsetBases(final int start, final int end) {
        final int intervalStart = interval.getStart();
        final int intervalEnd = interval.getEnd();
        Utils.validateArg(start >= intervalStart, "the start of the subset must be within this reference bases interval");
        Utils.validateArg(end <= intervalEnd, "the end of the subset must be within this reference bases interval");
        Utils.validateArg(start <= end, "the start must be less or equal to the end");
        return Arrays.copyOfRange(bases, start - intervalStart, bases.length  + end - intervalEnd);
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    /**
     * getSubset returns only the bases of the interval passed in.
     * @param subsetInterval, the subset to be returned
     * @return the subset of ReferenceBases
     */
    public ReferenceBases getSubset(SimpleInterval subsetInterval) {
        if (!this.interval.contains(subsetInterval)) {
            throw new GATKException("Reference doesn't match input interval (asked for "+subsetInterval.toString()+" but we have "+this.interval+")");
        }
        int start = subsetInterval.getStart() - this.interval.getStart();
        int end = subsetInterval.getEnd() - this.interval.getStart();
        return new ReferenceBases(Arrays.copyOfRange(this.bases, start, end + 1), subsetInterval);
    }
}
