package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

import java.util.Arrays;

/**
 * General heterogeneous ploidy model.
 */
public final class HeterogeneousPloidyModel implements PloidyModel {

    private final SampleList sampleList;

    private final int[] ploidies;

    private final int ploidySum;

    private final boolean isHomogeneous;

    public HeterogeneousPloidyModel(final SampleList sampleList, final int[] ploidies) {
        this.sampleList = Utils.nonNull(sampleList, "the sample list cannot be null");
        this.ploidies = Utils.nonNull(ploidies, "the ploidies cannot be null").clone();
        Utils.validateArg(sampleList.numberOfSamples() == ploidies.length, "sample-list and ploidy array length must match");
        Arrays.stream(ploidies).forEach(p -> Utils.validateArg(p >= 0, "no ploidy can be less than 0"));
        ploidySum = (int) MathUtils.sum(ploidies);
        isHomogeneous = ploidies.length == 0 || Arrays.stream(ploidies).allMatch(p -> p == ploidies[0]);
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