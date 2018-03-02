package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

/**
 * A genotyped integer copy-number segment.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class IntegerCopyNumberSegment implements Locatable {
    private final SimpleInterval interval;
    private final IntegerCopyNumberState callIntegerCopyNumberState;
    private final IntegerCopyNumberState baselineIntegerCopyNumberState;
    private final int numPoints;
    private final double qualitySomeCalled;
    private final double qualityAllCalled;
    private final double qualityStart;
    private final double qualityEnd;

    public IntegerCopyNumberSegment(final SimpleInterval interval,
                                    final IntegerCopyNumberState callIntegerCopyNumberState,
                                    final IntegerCopyNumberState baselineIntegerCopyNumberState,
                                    final int numPoints,
                                    final double qualitySomeCalled,
                                    final double qualityAllCalled,
                                    final double qualityStart,
                                    final double qualityEnd) {
        this.interval = Utils.nonNull(interval, "The interval for the segment must be non-null.");
        this.callIntegerCopyNumberState = Utils.nonNull(callIntegerCopyNumberState,
                "The call integer copy-number state for the segment must be non-null.");
        this.baselineIntegerCopyNumberState = Utils.nonNull(baselineIntegerCopyNumberState,
                "The baseline integer copy-number state for the segment must be non-null.");
        this.numPoints = ParamUtils.isPositive(numPoints,
                "Number of points in the segment must be positive.");
        this.qualitySomeCalled = ParamUtils.isPositiveOrZero(qualitySomeCalled,
                "The phred-scaled quality of \"some points called\" must be non-negative.");
        this.qualityAllCalled = ParamUtils.isPositiveOrZero(qualityAllCalled,
                "The phred-scaled quality of \"all points called\" must be non-negative.");
        this.qualityStart = ParamUtils.isPositiveOrZero(qualityStart,
                "The phred-scaled quality of \"segment start\" must be non-negative.");
        this.qualityEnd = ParamUtils.isPositiveOrZero(qualityEnd,
                "The phred-scaled quality of \"segment end\" must be non-negative");
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

    public IntegerCopyNumberState getCallIntegerCopyNumberState() {
        return callIntegerCopyNumberState;
    }

    public IntegerCopyNumberState getBaselineIntegerCopyNumberState() {
        return baselineIntegerCopyNumberState;
    }

    public int getNumPoints() {
        return numPoints;
    }

    public double getQualitySomeCalled() {
        return qualitySomeCalled;
    }

    public double getQualityAllCalled() {
        return qualityAllCalled;
    }

    public double getQualityStart() {
        return qualityStart;
    }

    public double getQualityEnd() {
        return qualityEnd;
    }

    @Override
    public boolean equals(final Object other) {
        if (this == other) {
            return true;
        }
        if (other == null || getClass() != other.getClass()) {
            return false;
        }
        final IntegerCopyNumberSegment that = (IntegerCopyNumberSegment) other;
        return (numPoints == that.numPoints &&
                Double.compare(that.qualitySomeCalled, qualitySomeCalled) == 0 &&
                Double.compare(that.qualityAllCalled, qualityAllCalled) == 0 &&
                Double.compare(that.qualityStart, qualityStart) == 0 &&
                Double.compare(that.qualityEnd, qualityEnd) == 0 &&
                interval.equals(that.interval) &&
                callIntegerCopyNumberState.equals(that.callIntegerCopyNumberState) &&
                baselineIntegerCopyNumberState.equals(that.baselineIntegerCopyNumberState));
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = interval.hashCode();
        result = 31 * result + callIntegerCopyNumberState.hashCode();
        result = 31 * result + baselineIntegerCopyNumberState.hashCode();
        result = 31 * result + numPoints;
        temp = Double.doubleToLongBits(qualitySomeCalled);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(qualityAllCalled);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(qualityStart);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(qualityEnd);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    @Override
    public String toString() {
        return "IntegerCopyNumberSegment{" +
                "interval=" + interval +
                ", callIntegerCopyNumberState=" + callIntegerCopyNumberState +
                ", baselineIntegerCopyNumberState=" + baselineIntegerCopyNumberState +
                ", numPoints=" + numPoints +
                ", qualitySomeCalled=" + qualitySomeCalled +
                ", qualityAllCalled=" + qualityAllCalled +
                ", qualityStart=" + qualityStart +
                ", qualityEnd=" + qualityEnd +
                '}';
    }
}
