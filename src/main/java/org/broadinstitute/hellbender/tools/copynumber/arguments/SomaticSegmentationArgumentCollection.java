package org.broadinstitute.hellbender.tools.copynumber.arguments;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SomaticSegmentationArgumentCollection implements Serializable {
    public static final long serialVersionUID = 1L;

    //segmentation argument names
    public static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME = "maximum-number-of-segments-per-chromosome";
    public static final String KERNEL_VARIANCE_COPY_RATIO_LONG_NAME = "kernel-variance-copy-ratio";
    public static final String KERNEL_VARIANCE_ALLELE_FRACTION_LONG_NAME = "kernel-variance-allele-fraction";
    public static final String KERNEL_SCALING_ALLELE_FRACTION_LONG_NAME = "kernel-scaling-allele-fraction";
    public static final String KERNEL_APPROXIMATION_DIMENSION_LONG_NAME = "kernel-approximation-dimension";
    public static final String WINDOW_SIZE_LONG_NAME = "window-size";
    public static final String NUMBER_OF_CHANGEPOINTS_PENALTY_FACTOR_LONG_NAME = "number-of-changepoints-penalty-factor";

    @Argument(
            doc = "Maximum number of segments allowed per chromosome.",
            fullName = MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME,
            minValue = 1,
            optional = true
    )
    public int maxNumSegmentsPerChromosome = 1000;

    @Argument(
            doc = "Variance of Gaussian kernel for copy-ratio segmentation, if performed.  If zero, a linear kernel will be used.",
            fullName = KERNEL_VARIANCE_COPY_RATIO_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    public double kernelVarianceCopyRatio = 0.;

    @Argument(
            doc = "Variance of Gaussian kernel for allele-fraction segmentation, if performed.  If zero, a linear kernel will be used.",
            fullName = KERNEL_VARIANCE_ALLELE_FRACTION_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    public double kernelVarianceAlleleFraction = 0.025;

    @Argument(
            doc = "Relative scaling S of the kernel K_AF for allele-fraction segmentation to the kernel K_CR for copy-ratio segmentation.  " +
                    "If multidimensional segmentation is performed, the total kernel used will be K_CR + S * K_AF.",
            fullName = KERNEL_SCALING_ALLELE_FRACTION_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    public double kernelScalingAlleleFraction = 1.0;

    @Argument(
            doc = "Dimension of the kernel approximation.  A subsample containing this number of data points " +
                    "will be used to construct the approximation for each chromosome.  " +
                    "If the total number of data points in a chromosome is greater " +
                    "than this number, then all data points in the chromosome will be used.  " +
                    "Time complexity scales quadratically and space complexity scales linearly with this parameter.",
            fullName = KERNEL_APPROXIMATION_DIMENSION_LONG_NAME,
            minValue = 1,
            optional = true
    )
    public int kernelApproximationDimension = 100;

    @Argument(
            doc = "Window sizes to use for calculating local changepoint costs.  " +
                    "For each window size, the cost for each data point to be a changepoint will be calculated " +
                    "assuming that the point demarcates two adjacent segments of that size.  " +
                    "Including small (large) window sizes will increase sensitivity to small (large) events.  " +
                    "Duplicate values will be ignored.",
            fullName = WINDOW_SIZE_LONG_NAME,
            minValue = 1,
            optional = true
    )
    public List<Integer> windowSizes = new ArrayList<>(Arrays.asList(8, 16, 32, 64, 128, 256));

    @Argument(
            doc = "Factor A for the penalty on the number of changepoints per chromosome for segmentation.  " +
                    "Adds a penalty of the form A * C * [1 + log (N / C)], " +
                    "where C is the number of changepoints in the chromosome, " +
                    "to the cost function for each chromosome.  " +
                    "Must be non-negative.",
            fullName = NUMBER_OF_CHANGEPOINTS_PENALTY_FACTOR_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    public double numChangepointsPenaltyFactor = 1.;
}
