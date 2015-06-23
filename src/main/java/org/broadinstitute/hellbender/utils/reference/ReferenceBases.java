package org.broadinstitute.hellbender.utils.reference;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;
import java.util.Arrays;

public class ReferenceBases implements Serializable {
    private static final long serialVersionUID = 1L;

    private final byte[] bases;
    private final SimpleInterval interval;

    public ReferenceBases( final byte[] bases, final SimpleInterval interval ) {
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

    public ReferenceBases getSubset(SimpleInterval interval) {
        if (!this.interval.contains(interval)) {
            throw new GATKException("Reference doesn't match read...");
        }
        int start = interval.getStart() - this.interval.getStart();
        int end = interval.getEnd() - this.interval.getStart();
        return new ReferenceBases(Arrays.copyOfRange(this.bases, start, end + 1/* ??? */), interval);
    }
}
