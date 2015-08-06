package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.genotyper.SampleList;

/**
 * General heterogeneous ploidy model.
 *
 * <p>
 *     Currenly only avaialable for testing but will be promoted at some point and have its own unit test.
 * </p>
 */
public final class HeterogeneousPloidyModel implements PloidyModel {

    private final SampleList sampleList;

    private final int[] ploidies;

    private final int ploidySum;

    private final boolean isHomogeneous;

    public HeterogeneousPloidyModel(final SampleList sampleList, final int[] ploidies) {
        if (sampleList == null) {
            throw new IllegalArgumentException("the sample list cannot be null");
        }
        if (ploidies == null) {
            throw new IllegalArgumentException("the ploidies cannot be null");
        }
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
        if (sampleIndex < 0 || sampleIndex > ploidies.length) {
            throw new IllegalArgumentException("invalid sample index: " + sampleIndex);
        }
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
        return sampleList.indexOfSample(sample);
    }

    @Override
    public String getSample(final int sampleIndex) {
        return sampleList.getSample(sampleIndex);
    }
}