package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

/**
 * General heterogeneous ploidy model.
 */
public final class HeterogeneousPloidyModel implements PloidyModel {

    private final SampleList sampleList;

    private final int[] ploidies;

    private final int ploidySum;

    private final boolean isHomogeneous;

    public HeterogeneousPloidyModel(final SampleList sampleList, final int[] ploidies) {
        Utils.nonNull(sampleList, "the sample list cannot be null");
        Utils.nonNull(ploidies, "the ploidies cannot be null");
        if (sampleList.numberOfSamples() != ploidies.length) {
            throw new IllegalArgumentException("sample-list and ploidy array length must match");
        }

        this.ploidies = ploidies.clone();

        int ploidySum = 0;
        for (int i = 0; i < ploidies.length; i++) {
            final int p = this.ploidies[i];
            if (p < 0) {
                throw new IllegalArgumentException("no ploidy can be less than 0");
            }
            ploidySum += p;
        }
        this.ploidySum = ploidySum;
        isHomogeneous = ploidies.length == 0 || ploidies.length * this.ploidies[0] == ploidySum;
        this.sampleList = sampleList;
    }

    @Override
    public int samplePloidy(final int sampleIndex) {
        Utils.validIndex(sampleIndex, ploidies.length);
        return ploidies[sampleIndex];
    }

    @Override
    public boolean isHomogeneous() {
        return isHomogeneous;
    }

    @Override
    public int totalPloidy() {
        return ploidySum;
    }

    @Override
    public int numberOfSamples() {
        return ploidies.length;
    }

    @Override
    public int indexOfSample(final String sample) {
        Utils.nonNull(sample);
        return sampleList.indexOfSample(sample);
    }

    @Override
    public String getSample(final int sampleIndex) {
        Utils.validIndex(sampleIndex, numberOfSamples());
        return sampleList.getSample(sampleIndex);
    }
}