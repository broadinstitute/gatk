package org.broadinstitute.hellbender.utils.hmm.segmentation;

import org.broadinstitute.hellbender.tools.exome.Segment;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;
import java.util.List;

/**
 * Represents a sequence of hidden states with the same call as a segment.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class HiddenStateSegment<S, T extends Target> extends Segment<S> implements Serializable {

    private static final long serialVersionUID = 1135086084964809775L;

    private final double exactQuality;

    private final double someQuality;

    private final double startQuality;

    private final double endQuality;

    private final double eventQuality;

    private final double stdev;

    /**
     * Returns the "exact quality" (EQ): 1 - posterior probability of the hidden state to be the call state in
     * all the target in the segment in Phred scale.
     *
     * @return any {@code double} value.
     */
    public double getExactQuality() {
        return exactQuality;
    }

    /**
     * Returns the "some quality" (SQ): 1 - posterior probability of the hidden state
     * to be the same as the segment call state in at least one target enclosed in the segment in Phred scale.
     * @return any {@code double} value.
     */
    public double getSomeQuality() {
        return someQuality;
    }

    /**
     * Returns the "start quality" (LQ): 1 - posterior probability that the actual event described by
     * the segment starts at the first target enclosed in the segment in Phred scale.
     * <p>
     * More concretely, this quality is equal 1 - posterior probability that the hidden state is not the segment
     * call state in the target just before the segment and the hidden state is the segment call state
     * in the first target inside the segment in Phred scale.
     * </p>
     * @return any {@code double} value.
     */
    public double getStartQuality() {
        return startQuality;
    }

    /**
     * Returns the "end quality" (RQ): 1 - posterior probability that the actual event described by
     * the segment finishes at the last target enclosed in the segment in Phred scale.
     * <p>
     * More concretely, this quality is equal 1 - posterior probability that the hidden state is not the segment
     * call state in the target that follows the segment and is the segment's call the hidden state is the last target
     * in the segment in Phred scale.
     * </p>
     * @return any {@code double} value.
     */
    public double getEndQuality() {
        return endQuality;
    }

    /**
     * Returns the "event quality" (NDP): 1 - posterior probability that the hidden
     * state is different from the the "reference state" across the whole segment
     * in Phred scale.
     *
     * The reference state is determined by the context. For a 3-state HMM, the reference
     * state is the copy-neutral state.
     *
     * @return any {@code double} value.
     */
    public double getEventQuality() {
        return eventQuality;
    }

    /**
     * Returns the standard deviation of emission model data (defined in context) across targets in
     * the segment.
     */
    public double getStdev() {
        return stdev;
    }

    /**
     * Creates a new segment given all its properties and overlapped targets.
     *
     * @param targets the overlapped targets.
     * @param mean the mean coverage value accross targets.
     * @param call the segment call.
     * @param exactQuality the phred scaled EQ score.
     * @param someQuality the phred scaled SQ score.
     * @param startQuality the phred scaled LQ score.
     * @param stopQuality the phred scaled SQ score.
     * @param eventQuality the phred scaled NDQ score.
     */
    public HiddenStateSegment(final List<T> targets, final double mean, final double stdev,
                              final S call,
                              final double exactQuality, final double someQuality, final double startQuality,
                              final double stopQuality, final double eventQuality) {
        super(targets, mean, call);
        this.stdev = stdev;
        this.exactQuality = exactQuality;
        this.someQuality = someQuality;
        this.startQuality = startQuality;
        endQuality = stopQuality;
        this.eventQuality = eventQuality;
    }


    /**
     * Creates a new segment given all its properties and overlapped targets.
     *
     * @param interval the overlapped targets.
     * @param mean the mean coverage value accross targets.
     * @param call the segment call.
     * @param exactQuality the phred scaled EQ score.
     * @param someQuality the phred scaled SQ score.
     * @param startQuality the phred scaled LQ score.
     * @param stopQuality the phred scaled SQ score.
     * @param eventQuality the phred scaled NDQ score.
     */
    public HiddenStateSegment(final SimpleInterval interval, final int targetCount, final double mean, final double stdev,
                              final S call,
                              final double exactQuality, final double someQuality, final double startQuality,
                              final double stopQuality, final double eventQuality) {
        super(interval, targetCount, mean, call);
        this.stdev = stdev;
        this.exactQuality = exactQuality;
        this.someQuality = someQuality;
        this.startQuality = startQuality;
        endQuality = stopQuality;
        this.eventQuality = eventQuality;
    }

    public String toString() {
        return String.format("{interval: %s, targetCount: %d, call: %s, mean: %.2f, stdev: %.2f, EQ: %.2f, SQ: %.2f, LQ: %.2f, RQ: %.2f, NDQ: %.2f}",
                interval, targetCount, call, mean,
                stdev, exactQuality, someQuality, startQuality, endQuality, eventQuality);
    }

}
