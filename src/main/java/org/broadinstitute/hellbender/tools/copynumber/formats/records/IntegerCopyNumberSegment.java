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
    private final int numSpanningIntervals;
    private final double someQuality;
    private final double exactQuality;
    private final double startQuality;
    private final double endQuality;

    public IntegerCopyNumberSegment(final SimpleInterval interval,
                                    final IntegerCopyNumberState callIntegerCopyNumberState,
                                    final IntegerCopyNumberState baselineIntegerCopyNumberState,
                                    final int numSpanningIntervals,
                                    final double someQuality,
                                    final double exactQuality,
                                    final double startQuality,
                                    final double endQuality) {
        this.interval = Utils.nonNull(interval, "The interval for the segment must be non-null.");
        this.callIntegerCopyNumberState = Utils.nonNull(callIntegerCopyNumberState,
                "The call integer copy-number state for the segment must be non-null.");
        this.baselineIntegerCopyNumberState = Utils.nonNull(baselineIntegerCopyNumberState,
                "The baseline integer copy-number state for the segment must be non-null.");
        this.numSpanningIntervals = ParamUtils.isPositive(numSpanningIntervals,
                "Number of intervals spanned by the integer copy-number segment must be positive.");
        this.someQuality = ParamUtils.isPositiveOrZero(someQuality, "Some quality must be non-negative.");
        this.exactQuality = ParamUtils.isPositiveOrZero(exactQuality, "Exact quality must be non-negative.");
        this.startQuality = ParamUtils.isPositiveOrZero(startQuality, "Left breakpoint quality must be non-negative.");
        this.endQuality = ParamUtils.isPositiveOrZero(endQuality, "Right breakpoint quality must be non-negative");
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

    public int getNumSpanningIntervals() {
        return numSpanningIntervals;
    }

    public double getSomeQuality() {
        return someQuality;
    }

    public double getExactQuality() {
        return exactQuality;
    }

    public double getStartQuality() {
        return startQuality;
    }

    public double getEndQuality() {
        return endQuality;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        IntegerCopyNumberSegment that = (IntegerCopyNumberSegment) o;

        if (numSpanningIntervals != that.numSpanningIntervals) return false;
        if (Double.compare(that.someQuality, someQuality) != 0) return false;
        if (Double.compare(that.exactQuality, exactQuality) != 0) return false;
        if (Double.compare(that.startQuality, startQuality) != 0) return false;
        if (Double.compare(that.endQuality, endQuality) != 0) return false;
        if (!interval.equals(that.interval)) return false;
        if (!callIntegerCopyNumberState.equals(that.callIntegerCopyNumberState)) return false;
        return baselineIntegerCopyNumberState.equals(that.baselineIntegerCopyNumberState);
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = interval.hashCode();
        result = 31 * result + callIntegerCopyNumberState.hashCode();
        result = 31 * result + baselineIntegerCopyNumberState.hashCode();
        result = 31 * result + numSpanningIntervals;
        temp = Double.doubleToLongBits(someQuality);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(exactQuality);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(startQuality);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(endQuality);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    @Override
    public String toString() {
        return "IntegerCopyNumberSegment{" +
                "interval=" + interval +
                ", callIntegerCopyNumberState=" + callIntegerCopyNumberState +
                ", baselineIntegerCopyNumberState=" + baselineIntegerCopyNumberState +
                ", numSpanningIntervals=" + numSpanningIntervals +
                ", someQuality=" + someQuality +
                ", exactQuality=" + exactQuality +
                ", startQuality=" + startQuality +
                ", endQuality=" + endQuality +
                '}';
    }
}
