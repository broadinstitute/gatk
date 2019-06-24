package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import org.broadinstitute.hellbender.utils.MathUtils;

/**
 * Created by valentin on 12/4/16.
 */
class STREvaluationCounter {
    private int[] values = new int[STREvaluationClass.values().length];

    STREvaluationCounter() {
    }

    public void add(final STREvaluationClass count) {
        values[count.ordinal()]++;
    }

    public int get(final STREvaluationClass evalClass) {
        return values[evalClass.ordinal()];
    }

    public long total() {
        return MathUtils.sum(values);
    }
}
