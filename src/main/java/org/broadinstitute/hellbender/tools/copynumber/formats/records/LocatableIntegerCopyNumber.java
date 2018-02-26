package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * This class represents an instance of {@link IntegerCopyNumberState} associated with a {@link Locatable} genomic
 * interval.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class LocatableIntegerCopyNumber implements Locatable {
    private final SimpleInterval interval;
    private final IntegerCopyNumberState integerCopyNumberState;

    public LocatableIntegerCopyNumber(final SimpleInterval interval,
                                      final IntegerCopyNumberState integerCopyNumberState) {
        this.interval = Utils.nonNull(interval);
        this.integerCopyNumberState = Utils.nonNull(integerCopyNumberState);
    }

    public IntegerCopyNumberState getIntegerCopyNumberState() {
        return integerCopyNumberState;
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    @Override
    public boolean equals(final Object other) {
        if (this == other) {
            return true;
        }
        if (other == null || getClass() != other.getClass()) {
            return false;
        }
        LocatableIntegerCopyNumber that = (LocatableIntegerCopyNumber) other;
        return (interval.equals(that.interval) &&
                integerCopyNumberState.equals(that.integerCopyNumberState));
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + integerCopyNumberState.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return "LocatableIntegerCopyNumber{" +
                "interval=" + interval +
                ", integerCopyNumberState=" + integerCopyNumberState +
                '}';
    }
}
