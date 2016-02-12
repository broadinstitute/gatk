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

    public byte[] getBases() {
        return bases;
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
