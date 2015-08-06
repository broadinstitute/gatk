package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.genotyper.SampleList;

/**
 * Information about the number of chromosome per sample at a given location.
 *
 */
public interface PloidyModel extends SampleList {

    /**
     * Return the assumed ploidy for a sample given its index.
     *
     * @param sampleIndex target sample index.
     * @return 0 or greater.
     */
    public int samplePloidy(final int sampleIndex);

    /**
     * Checks whether the ploidy is homogeneous across all samples.
     *
     * @return {@code true} if all samples has the same ploidy.
     */
    public boolean isHomogeneous();

    /**
     * Sum of all ploidy across all samples.
     * <p>
     *     It must match the sum of all ploidies across samples.
     * </p>
     *
     * @return 0 or greater.
     */
    public int totalPloidy();
}
