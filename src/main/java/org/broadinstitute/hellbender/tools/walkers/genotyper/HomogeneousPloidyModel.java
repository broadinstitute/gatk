package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

/**
* {@link PloidyModel} implementation tailored to work with a homogeneous constant ploidy
* across samples and positions.
*/
public final class HomogeneousPloidyModel implements PloidyModel {

    private SampleList sampleList;

    private final int ploidy;

    /**
     * Constructs a homogeneous ploidy model given the sample list and ploidy.
     *
     * @param sampleList the sample list.
     * @param ploidy the common ploidy for all samples in {@code samples}.
     *
     * @throws IllegalArgumentException if {@code sampleList} is {@code null},
     *    or ploidy is 0 or less.
     */
    public HomogeneousPloidyModel(final SampleList sampleList, final int ploidy) {
        Utils.nonNull(sampleList);
        Utils.validateArg(ploidy > 0, "does not support negative ploidy");

        this.ploidy = ploidy;
        this.sampleList = sampleList;
    }

    @Override
    public int numberOfSamples() {
        return sampleList.numberOfSamples();
    }

    @Override
    public String getSample(final int index) {
        Utils.validIndex(index, numberOfSamples());
        return sampleList.getSample(index);
    }

    @Override
    public int indexOfSample(final String sample) {
        Utils.nonNull(sample);
        return sampleList.indexOfSample(sample);
    }

    @Override
    public int samplePloidy(final int sampleIndex) {
        Utils.validIndex(sampleIndex, numberOfSamples());
        return ploidy;
    }

    @Override
    public boolean isHomogeneous() {
        return true;
    }

    @Override
    public int totalPloidy() {
        return ploidy * sampleList.numberOfSamples();
    }

}
