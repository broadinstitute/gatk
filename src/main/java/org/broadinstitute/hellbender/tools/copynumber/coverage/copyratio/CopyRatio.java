package org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

public class CopyRatio implements Locatable {
    private final SimpleInterval interval;
    private final double log2CopyRatioValue;

    public CopyRatio(final SimpleInterval interval,
                     final double log2CopyRatioValue) {
        Utils.nonNull(interval);
        this.interval = interval;
        this.log2CopyRatioValue = log2CopyRatioValue;
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

    public SimpleInterval getInterval() {
        return interval;
    }

    public double getLog2CopyRatioValue() {
        return log2CopyRatioValue;
    }

    /**
     * The midpoint is used to characterize the interval for the purposes of determining overlaps
     * so that each copy-ratio interval will be uniquely contained in a single segment.
     */
    public SimpleInterval getMidpoint() {
        final int midPoint = (getStart() + getEnd()) / 2;
        return new SimpleInterval(interval.getContig(), midPoint, midPoint);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final CopyRatio copyRatio = (CopyRatio) o;
        return Double.compare(copyRatio.log2CopyRatioValue, log2CopyRatioValue) == 0 && interval.equals(copyRatio.interval);
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = interval.hashCode();
        temp = Double.doubleToLongBits(log2CopyRatioValue);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    @Override
    public String toString() {
        return "CopyRatio{" +
                "interval=" + interval +
                ", log2CopyRatioValue=" + log2CopyRatioValue +
                '}';
    }
}
