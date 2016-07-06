package org.broadinstitute.hellbender.tools.picard.interval;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * @author mccowan
 */
public final class IntervalListScatterer {

    public enum Mode {
        /**
         * A simple scatter approach in which all output intervals have size equal to the total base count of the source list divide by the
         * scatter count (except for possible variance in the final interval list).
         */
        INTERVAL_SUBDIVISION,
        /**
         * A scatter approach that differs from {@link Mode#INTERVAL_SUBDIVISION} in a few ways.
         * <ol>
         * <li>No interval will be subdivided, and consequently, the requested scatter count is an upper bound of scatter count, not a
         * guarantee as to how many {@link IntervalList}s will be produced (e.g., if scatterCount = 10 but there is only one input interval,
         * only 1 interval list will be emitted).</li>
         * <li>When an interval would otherwise be split, it is instead deferred to the next scatter list.</li>
         * <li>The "target width" of each scatter list may be wider than what is computed for {@link Mode#INTERVAL_SUBDIVISION}.
         * Specifically, if the widest interval in the source interval list is larger than what would otherwise be the target width, that
         * interval's width is used.<br/><br/>The reasoning for this is that this approach produces more consistently-sized interval lists,
         * which is one of the objectives of scattering.</li>
         * </ol>
         */
        BALANCING_WITHOUT_INTERVAL_SUBDIVISION
    }

    private final Mode mode;

    public IntervalListScatterer(final Mode mode) {this.mode = mode;}

    private int deduceIdealSplitLength(final IntervalList uniquedList, final int scatterCount) {
        final int splitWidth = Math.max((int) Math.floor(uniquedList.getBaseCount() / (1.0 * scatterCount)), 1);
        switch (mode) {
            case INTERVAL_SUBDIVISION:
                return splitWidth;
            case BALANCING_WITHOUT_INTERVAL_SUBDIVISION:
                final int widestIntervalLength = Collections.max(uniquedList.getIntervals(), new Comparator<Interval>() {
                    @Override
                    public int compare(final Interval o1, final Interval o2) {
                        return Integer.valueOf(o1.length()).compareTo(o2.length());
                    }
                }).length();

                // There is no purpose to splitting more granularly than the widest interval, so do not.
                return Math.max(widestIntervalLength, splitWidth);
            default:
                throw new IllegalStateException();
        }
    }

    public List<IntervalList> scatter(final IntervalList uniquedIntervalList, final int scatterCount) {
        return scatter(uniquedIntervalList, scatterCount, false);
    }

    public List<IntervalList> scatter(final IntervalList sourceIntervalList, final int scatterCount, final boolean isUniqued) {
        Utils.validateArg(scatterCount >= 1, "scatterCount < 1");

        final IntervalList uniquedList = isUniqued ? sourceIntervalList : sourceIntervalList.uniqued();
        final long idealSplitLength = deduceIdealSplitLength(uniquedList, scatterCount);

        final List<IntervalList> accumulatedIntervalLists = new ArrayList<>();

        IntervalList runningIntervalList = new IntervalList(uniquedList.getHeader());
        final ArrayDeque<Interval> intervalQueue = new ArrayDeque<>(uniquedList.getIntervals());

        while (!intervalQueue.isEmpty() && accumulatedIntervalLists.size() < scatterCount - 1) {
            final Interval interval = intervalQueue.pollFirst();
            final long projectedSize = runningIntervalList.getBaseCount() + interval.length();
            if (projectedSize <= idealSplitLength) {
                runningIntervalList.add(interval);
            } else {
                switch (mode) {
                    case INTERVAL_SUBDIVISION:
                        final int amountToConsume = (int) (idealSplitLength - runningIntervalList.getBaseCount());
                        final Interval left = new Interval(
                                interval.getContig(),
                                interval.getStart(),
                                interval.getStart() + amountToConsume - 1,
                                interval.isNegativeStrand(),
                                interval.getName()
                        );
                        final Interval right = new Interval(
                                interval.getContig(),
                                interval.getStart() + amountToConsume,
                                interval.getEnd(),
                                interval.isNegativeStrand(),
                                interval.getName()
                        );
                        runningIntervalList.add(left);

                        // Push back the excess back onto our queue for reconsideration.
                        intervalQueue.addFirst(right);
                        break;

                    case BALANCING_WITHOUT_INTERVAL_SUBDIVISION:
                        if (runningIntervalList.getIntervals().isEmpty()) {
                            runningIntervalList.add(interval);
                        } else {
                            // Push this interval into the next scatter; re-inject it into the queue, then advance the scatter.
                            intervalQueue.addFirst(interval);
                            accumulatedIntervalLists.add(runningIntervalList.uniqued());
                            runningIntervalList = new IntervalList(uniquedList.getHeader());
                        }
                        break;
                }
            }

            if (runningIntervalList.getBaseCount() >= idealSplitLength) {
                accumulatedIntervalLists.add(runningIntervalList.uniqued());
                runningIntervalList = new IntervalList(uniquedList.getHeader());
            }
        }

        // Flush the remaining intervals into the last split.
        while (!intervalQueue.isEmpty()) {
            runningIntervalList.add(intervalQueue.pollFirst());
        }
        if (!runningIntervalList.getIntervals().isEmpty()) {
            accumulatedIntervalLists.add(runningIntervalList.uniqued());
        }

        return accumulatedIntervalLists;
    }
}
