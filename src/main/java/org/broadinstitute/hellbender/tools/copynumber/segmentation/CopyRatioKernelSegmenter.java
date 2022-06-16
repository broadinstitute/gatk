package org.broadinstitute.hellbender.tools.copynumber.segmentation;

import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.List;

/**
 * Segments copy-ratio data using kernel segmentation.  Segments do not span chromosomes.
 *
 * Refactored to be a thin wrapper around the {@link MultisampleMultidimensionalKernelSegmenter} so that
 * preexisting test code could be reused for testing that class; do not use this class in production code.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioKernelSegmenter {
    private static final double DUMMY_ALLELE_FRACTION_KERNEL_VARIANCE = 1.;
    private static final double DUMMY_KERNEL_SCALING_ALLELE_FRACTION = 1.;

    private final MultisampleMultidimensionalKernelSegmenter segmenter;

    /**
     * @param denoisedCopyRatios  in log2 space
     */
    public CopyRatioKernelSegmenter(final CopyRatioCollection denoisedCopyRatios) {
        Utils.nonNull(denoisedCopyRatios);
        segmenter = new MultisampleMultidimensionalKernelSegmenter(
                Collections.singletonList(denoisedCopyRatios),
                Collections.singletonList(new AllelicCountCollection(denoisedCopyRatios.getMetadata(), Collections.emptyList())));
    }

    /**
     * @param kernelVariance    variance of the Gaussian kernel; if zero, a linear kernel is used instead
     */
    public SimpleIntervalCollection findSegmentation(final int maxNumSegmentsPerChromosome,
                                                     final double kernelVariance,
                                                     final int kernelApproximationDimension,
                                                     final List<Integer> windowSizes,
                                                     final double numChangepointsPenaltyLinearFactor,
                                                     final double numChangepointsPenaltyLogLinearFactor) {
        return segmenter.findSegmentation(
                maxNumSegmentsPerChromosome,
                kernelVariance,
                DUMMY_ALLELE_FRACTION_KERNEL_VARIANCE,
                DUMMY_KERNEL_SCALING_ALLELE_FRACTION,
                kernelApproximationDimension,
                windowSizes,
                numChangepointsPenaltyLinearFactor,
                numChangepointsPenaltyLogLinearFactor);
    }
}
