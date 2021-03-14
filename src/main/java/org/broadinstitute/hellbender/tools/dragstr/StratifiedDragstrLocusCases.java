package org.broadinstitute.hellbender.tools.dragstr;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

/**
 * Collection of {@link DragstrLocusCase} stratified by period and repeat length.
 * <p>
 *  Cases belonging to every possible combination can be added individually, and
 *   added or retrieved as a group using a {@link DragstrLocusCases} instance.
 * </p>
 */
final class StratifiedDragstrLocusCases {

    private int size;
    DragstrLocusCases[][] perPeriodAndRepeat;

    StratifiedDragstrLocusCases(final int maxPeriod, final int maxRepeats) {
        perPeriodAndRepeat = new DragstrLocusCases[maxPeriod][maxRepeats];
        for (int i = 0; i < maxPeriod; i++) {
            for (int j = 0; j < maxRepeats; j++) {
                perPeriodAndRepeat[i][j] = new DragstrLocusCases(i + 1, j + 1);
            }
        }
    }

    public static StratifiedDragstrLocusCases make(final int maxPeriod, final int maxRepeats) {
        ParamUtils.isPositive(maxPeriod, "max-period must be greater than 0");
        ParamUtils.isPositive(maxRepeats, "max-repeats must be greater than 0");
        return new StratifiedDragstrLocusCases(maxPeriod, maxRepeats);
    }

    public StratifiedDragstrLocusCases addAll(final StratifiedDragstrLocusCases other) {
        Utils.validate(perPeriodAndRepeat.length == other.perPeriodAndRepeat.length, "incompatible dimensions");
        for (int i = 0; i < perPeriodAndRepeat.length; i++) {
            Utils.validate(perPeriodAndRepeat[i].length == other.perPeriodAndRepeat[i].length, "invalid dimensions");
            for (int j = 0; j < perPeriodAndRepeat[i].length; j++) {
                perPeriodAndRepeat[i][j].addAll(other.perPeriodAndRepeat[i][j]);
            }
        }
        size += other.size;
        return this;
    }

    public int size() {
        return size;
    }

    public static StratifiedDragstrLocusCases merge(final StratifiedDragstrLocusCases... cazes) {
        if (cazes.length == 0) {
            return new StratifiedDragstrLocusCases(DragstrHyperParameters.DEFAULT_MAX_PERIOD, DragstrHyperParameters.DEFAULT_MAX_REPEAT_LENGTH);
        } else {
            final StratifiedDragstrLocusCases result = new StratifiedDragstrLocusCases(cazes[0].perPeriodAndRepeat.length, cazes[0].perPeriodAndRepeat[0].length);
            for (final StratifiedDragstrLocusCases col : cazes) {
                result.addAll(col);
            }
            return result;
        }
    }


    public DragstrLocusCases get(final int period, final int repeats) {
        final int periodIndex = Utils.validIndex(period - 1, perPeriodAndRepeat.length, "period is out of range");
        final int repeatIndex = Math.min(ParamUtils.isPositive(repeats, "repeats must be greater than 0 ") - 1,
                         perPeriodAndRepeat[periodIndex].length - 1);
        return perPeriodAndRepeat[periodIndex][repeatIndex];
    }

    /**
     * Adds all the cases in the input collection into this one segratting them
     * based on period and repeat length.
     * @param in the collection of cases to add in.
     */
    public void addAll(final DragstrLocusCases in) {
        final int periodIndex = in.getPeriod() - 1;
        final DragstrLocusCases cases = perPeriodAndRepeat[periodIndex]
                [Math.min(perPeriodAndRepeat[periodIndex].length - 1, in.getRepeatLength() - 1)];
        cases.addAll(in);
        size += in.size();
    }

    /**
     * Adds a single case.
     * @param caze the case to add.
     */
    public void add(final DragstrLocusCase caze) {
        final DragstrLocusCases cases = perPeriodAndRepeat[caze.getPeriod() - 1]
                [Math.min(perPeriodAndRepeat[0].length - 1, caze.getRepeatLength() - 1)];
        cases.add(caze);
        size++;
    }

    /**
     * Returns the subset of cases that pass a set of standard filters
     * based on depth, minimum MQ and maximum number of
     * supplementary alignments.
     *
     * <p>
     *  Changes on the returned collection won't have any effect on this one and <i>vice-versa</i>.
     * </p>
     */
    StratifiedDragstrLocusCases qualifyingOnly(final int minDepth, final int minMQ, final int maxSup) {
        final StratifiedDragstrLocusCases result = new StratifiedDragstrLocusCases(perPeriodAndRepeat.length, perPeriodAndRepeat[0].length);
        for (final DragstrLocusCases[] periodCases : perPeriodAndRepeat) {
            for (final DragstrLocusCases cases : periodCases) {
                result.addAll(cases.qualifyingOnly(minDepth, minMQ, maxSup));
            }
        }
        return result;
    }
}
