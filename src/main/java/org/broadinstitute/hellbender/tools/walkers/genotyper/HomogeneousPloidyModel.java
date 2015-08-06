package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.genotyper.SampleList;

/**
* {@link PloidyModel} implementation tailored to work with a homogeneous constant ploidy
* across samples and positions.
*
* @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
*/
public final class HomogeneousPloidyModel implements PloidyModel, SampleList {

    private SampleList sampleList;

    private final int ploidy;

    /**
     * Constructs a homogeneous ploidy model given the sample list and ploidy.
     *
     * @param samples the sample list.
     * @param ploidy the common ploidy for all samples in {@code samples}.
     *
     * @throws IllegalArgumentException if {@code samples} is {@code null},
     *    or ploidy is 0 or less.
     */
    public HomogeneousPloidyModel(final SampleList samples, final int ploidy) {
        if (ploidy <= 0) {
            throw new IllegalArgumentException("does not support negative ploidy");
        }
        this.ploidy = ploidy;

        sampleList = samples;
    }

    @Override
    public int numberOfSamples() {
        return sampleList.numberOfSamples();
    }

    @Override
    public String getSample(final int index) {
        return sampleList.getSample(index);
    }

    @Override
    public int indexOfSample(final String sample) {
        return sampleList.indexOfSample(sample);
    }

    @Override
    public int samplePloidy(final int sampleIndex) {
        checkSampleIndex(sampleIndex);
        return ploidy;
    }

    private void checkSampleIndex(final int sampleIndex) {
        if (sampleIndex < 0) {
            throw new IllegalArgumentException("the sample index cannot be negative: " + sampleIndex);
        }
        if (sampleIndex >= sampleList.numberOfSamples()) {
            throw new IllegalArgumentException("the sample index is equal or larger than the sample count: " + sampleIndex + " >= " + sampleList.numberOfSamples());
        }
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
