package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

/**
 * Represents a count at an interval.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class SimpleCount implements Locatable {

    private final SimpleInterval interval;
    private final int count;

    public SimpleCount(final SimpleInterval interval,
                       final int count) {
        this.interval = Utils.nonNull(interval);
        this.count = ParamUtils.isPositiveOrZero(count, "Can't construct SimpleCount with negative count.");
    }

    @Override
    public String getContig() { return interval.getContig(); }

    @Override
    public int getStart() { return interval.getStart(); }

    @Override
    public int getEnd() { return interval.getEnd(); }

    public SimpleInterval getInterval() { return interval; }

    public int getCount() { return count; }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final SimpleCount that = (SimpleCount) o;
        return count == that.count && interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + count;
        return result;
    }

    @Override
    public String toString() {
        return "SimpleCount{" +
                "interval=" + interval +
                ", count=" + count +
                '}';
    }
}
