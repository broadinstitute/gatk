package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;

/**
 * Represent segments.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 *
 * @param <C> call type.
 */
public class Segment<C> implements Locatable {

    protected SimpleInterval interval;
    protected C call;
    protected double mean;
    protected long targetCount;

    /**
     * Creates a segment given the list of targets it overlaps.
     * @param targets the targets the segment will contain.
     * @param mean the mean coverage across the targets.
     * @param call the call for that segment.
     * @throws IllegalArgumentException if {@code targets} or call is {@code null}, {@code targets} contains targets that do
     *  not belong to the same contig.
     */
    public Segment(final List<? extends Target> targets, final double mean, final C call) {
        this(IntervalUtils.getSpanningInterval(targets), Utils.nonNull(targets).size(), mean, call);
    }

    /**
     * Creates a segment given its interval and the number of targets within the interval.
     * @param interval the segment interval.
     * @param targetCount the number of targets in the segment.
     * @param mean the mean coverage across targets.
     * @param call the call for the segment.
     * @throws IllegalArgumentException if {@code interval} is {@code null}, {@code targetCount} is negative or
     * {@code call} is {@code null}.
     */
    public Segment(final SimpleInterval interval, final long targetCount, final double mean, final C call) {
        this.interval = Utils.nonNull(interval, "the input interval cannot be null");
        this.targetCount = ParamUtils.isPositiveOrZero(targetCount, "Number of original probes must be positive or zero.");
        this.call = Utils.nonNull(call, "the call cannot be null");
        this.mean = mean;
    }

    /**
     * Returns the segment call.
     *
     * @return never {@code null}.
     */
    public final C getCall() {
        return this.call;
    }

    /**
     * Returns the mean across the segment.
     *
     * @return a valid mean value between 0 and {@link Double#MAX_VALUE}.
     */
    public final double getMean() {
        return this.mean;
    }

    /**
     * Returns the segment coordinate interval.
     *
     * @return a valid interval.
     */
    public final SimpleInterval getInterval() {
        return this.interval;
    }

    /**
     * Returns the number of targets in the segment.
     *
     * @return never {@code null}.
     */
    public final long getTargetCount() {
        return this.targetCount;
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
}
