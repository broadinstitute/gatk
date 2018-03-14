package org.broadinstitute.hellbender.tools.copynumber.segmentation;

import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Segments copy-ratio data using kernel segmentation.  Segments do not span chromosomes.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioKernelSegmenter {
    private static final Logger logger = LogManager.getLogger(CopyRatioKernelSegmenter.class);

    private static final int MIN_NUM_POINTS_REQUIRED_PER_CHROMOSOME = 10;

    //Gaussian kernel for a specified variance; if variance is zero, use a linear kernel
    private static final Function<Double, BiFunction<Double, Double, Double>> KERNEL =
            variance -> variance == 0.
                    ? (x, y) -> x * y
                    : (x, y) -> FastMath.exp(-(x - y) * (x - y) / (2. * variance));

    private final CopyRatioCollection denoisedCopyRatios;
    private final Map<String, List<CopyRatio>> denoisedCopyRatiosPerChromosome;    //in log2 space

    /**
     * @param denoisedCopyRatios  in log2 space
     */
    public CopyRatioKernelSegmenter(final CopyRatioCollection denoisedCopyRatios) {
        Utils.nonNull(denoisedCopyRatios);
        this.denoisedCopyRatios = denoisedCopyRatios;
        denoisedCopyRatiosPerChromosome = denoisedCopyRatios.getRecords().stream()
                .collect(Collectors.groupingBy(
                        CopyRatio::getContig,
                        LinkedHashMap::new,
                        Collectors.mapping(Function.identity(), Collectors.toList())));
    }

    /**
     * Segments the internally held {@link CopyRatioCollection} using a separate {@link KernelSegmenter} for each chromosome.
     * @param kernelVariance    variance of the Gaussian kernel; if zero, a linear kernel is used instead
     */
    public CopyRatioSegmentCollection findSegmentation(final int maxNumChangepointsPerChromosome,
                                                       final double kernelVariance,
                                                       final int kernelApproximationDimension,
                                                       final List<Integer> windowSizes,
                                                       final double numChangepointsPenaltyLinearFactor,
                                                       final double numChangepointsPenaltyLogLinearFactor) {
        ParamUtils.isPositiveOrZero(maxNumChangepointsPerChromosome, "Maximum number of changepoints must be non-negative.");
        ParamUtils.isPositiveOrZero(kernelVariance, "Variance of Gaussian kernel must be non-negative (if zero, a linear kernel will be used).");
        ParamUtils.isPositive(kernelApproximationDimension, "Dimension of kernel approximation must be positive.");
        Utils.validateArg(windowSizes.stream().allMatch(ws -> ws > 0), "Window sizes must all be positive.");
        Utils.validateArg(new HashSet<>(windowSizes).size() == windowSizes.size(), "Window sizes must all be unique.");
        ParamUtils.isPositiveOrZero(numChangepointsPenaltyLinearFactor,
                "Linear factor for the penalty on the number of changepoints per chromosome must be non-negative.");
        ParamUtils.isPositiveOrZero(numChangepointsPenaltyLogLinearFactor,
                "Log-linear factor for the penalty on the number of changepoints per chromosome must be non-negative.");

        logger.info(String.format("Finding changepoints in %d data points and %d chromosomes...",
                denoisedCopyRatios.size(), denoisedCopyRatiosPerChromosome.size()));

        //loop over chromosomes, find changepoints, and create copy-ratio segments
        final List<CopyRatioSegment> segments = new ArrayList<>();
        for (final String chromosome : denoisedCopyRatiosPerChromosome.keySet()) {
            final List<CopyRatio> denoisedCopyRatiosInChromosome = denoisedCopyRatiosPerChromosome.get(chromosome);
            final int numDenoisedCopyRatiosInChromosome = denoisedCopyRatiosInChromosome.size();
            logger.info(String.format("Finding changepoints in %d data points in chromosome %s...",
                    numDenoisedCopyRatiosInChromosome, chromosome));

            if (numDenoisedCopyRatiosInChromosome < MIN_NUM_POINTS_REQUIRED_PER_CHROMOSOME) {
                logger.warn(String.format("Number of points in chromosome %s (%d) is less than that required (%d), skipping segmentation...",
                        chromosome, numDenoisedCopyRatiosInChromosome, MIN_NUM_POINTS_REQUIRED_PER_CHROMOSOME));
                final int start = denoisedCopyRatiosPerChromosome.get(chromosome).get(0).getStart();
                final int end = denoisedCopyRatiosPerChromosome.get(chromosome).get(numDenoisedCopyRatiosInChromosome - 1).getEnd();
                segments.add(new CopyRatioSegment(
                        new SimpleInterval(chromosome, start, end), denoisedCopyRatiosInChromosome));
                continue;
            }

            final List<Double> denoisedLog2CopyRatioValuesInChromosome = denoisedCopyRatiosInChromosome.stream()
                    .map(CopyRatio::getLog2CopyRatioValue)
                    .collect(Collectors.toList());
            final List<Integer> changepoints = new ArrayList<>(new KernelSegmenter<>(denoisedLog2CopyRatioValuesInChromosome)
                .findChangepoints(maxNumChangepointsPerChromosome, KERNEL.apply(kernelVariance), kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, KernelSegmenter.ChangepointSortOrder.INDEX));

            if (!changepoints.contains(numDenoisedCopyRatiosInChromosome)) {
                changepoints.add(numDenoisedCopyRatiosInChromosome - 1);
            }
            int previousChangepoint = -1;
            for (final int changepoint : changepoints) {
                final int start = denoisedCopyRatiosPerChromosome.get(chromosome).get(previousChangepoint + 1).getStart();
                final int end = denoisedCopyRatiosPerChromosome.get(chromosome).get(changepoint).getEnd();
                final List<CopyRatio> denoisedCopyRatiosInSegment = denoisedCopyRatiosInChromosome.subList(
                        previousChangepoint + 1, changepoint + 1);
                segments.add(new CopyRatioSegment(
                        new SimpleInterval(chromosome, start, end),
                        denoisedCopyRatiosInSegment));
                previousChangepoint = changepoint;
            }
        }
        logger.info(String.format("Found %d segments in %d chromosomes.", segments.size(), denoisedCopyRatiosPerChromosome.keySet().size()));
        return new CopyRatioSegmentCollection(denoisedCopyRatios.getMetadata(), segments);
    }
}
