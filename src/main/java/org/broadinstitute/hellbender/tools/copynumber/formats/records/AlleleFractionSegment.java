package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;

public class AlleleFractionSegment implements Locatable {
    private final SimpleInterval interval;
    private final int numPoints;

    public AlleleFractionSegment(final SimpleInterval interval,
                                 final int numPoints) {
        Utils.nonNull(interval);
        ParamUtils.isPositiveOrZero(numPoints, "Number of points must be non-negative.");
        this.interval = interval;
        this.numPoints = numPoints;
    }

    public AlleleFractionSegment(final SimpleInterval interval,
                                 final List<AllelicCount> allelicCounts) {
        Utils.nonNull(interval);
        Utils.nonNull(allelicCounts);
        this.interval = interval;
        numPoints = allelicCounts.size();
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

    public int getNumPoints() {
        return numPoints;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final AlleleFractionSegment that = (AlleleFractionSegment) o;
        return numPoints == that.numPoints &&
                interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        int result;
        result = interval.hashCode();
        result = 31 * result + numPoints;
        return result;
    }

    @Override
    public String toString() {
        return "AlleleFractionSegment{" +
                "interval=" + interval +
                ", numPoints=" + numPoints +
                '}';
    }
}
