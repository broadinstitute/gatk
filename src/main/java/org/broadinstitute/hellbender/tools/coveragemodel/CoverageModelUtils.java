package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Utility methods related to the target coverage model.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelUtils {
    /**
     * Estimate mean read depth $d$ according to a Poisson model for the distribution of read counts.
     *
     * The Poisson model is:
     *
     *     n_t \sim \mathrm{Poisson}(d c_t m_t p_t)
     *
     * where $n_t$ is the read count, $c_t$ is copy ratio, $p_t$ is ploidy, and $m_t$ is
     * the multiplicative bias (the estimation method is described in CNV-methods.pdf).
     *
     * Note: the mask array {@param mask} consists of 0s or 1s, and targets with 0 mask will be
     * effectively dropped out.
     *
     * Note: the product $c_t m_t p_t$ must be positive on all targets. Otherwise, they will be automatically
     * masked out (required for stability).
     *
     * TODO
     *
     * Note: This method prone to over-estimation because the Gaussian approximation of the Poisson distribution,
     * which is currently used in CoverageModelUtils.estimateReadDepthFromPoissonModel breaks down
     * for outlier read counts (if n >> depth). For the time being, we use medians which is a robust
     * statistic. In the future, we must find a better approximation for the Poisson distribution which
     * is both analytically tractable and robust.
     *
     * @param readCounts read counts
     * @param multiplicativeBiases estimate of total multiplicative bias
     * @param ploidies nominal germline ploidies of targets
     * @param copyRatio estimate of copy ratio
     * @param mask target masks (only 0 and 1 is expected)
     * @return read depth per base per homolog
     */
    public static double estimateReadDepthFromPoissonModel(final double[] readCounts,
                                                           final double[] multiplicativeBiases,
                                                           final int[] ploidies,
                                                           final double[] copyRatio,
                                                           final int[] mask) {
        Utils.nonNull(readCounts, "Read count array must be non-null");
        Utils.nonNull(multiplicativeBiases, "Multiplicative bias array must be non-null");
        Utils.nonNull(ploidies, "Ploidy array must be non-null");
        Utils.nonNull(copyRatio, "Copy ratio array must be non-null");
        Utils.nonNull(mask, "Mask array must be non-null");
        final int numTargets = readCounts.length;
        Utils.validateArg(numTargets > 0, "At least one target is required");
        Utils.validateArg(multiplicativeBiases.length == numTargets, "Multiplicative bias array must have the" +
                " same length as read counts array");
        Utils.validateArg(ploidies.length == numTargets, "Ploidies array must have the same length as read counts array");
        Utils.validateArg(copyRatio.length == numTargets, "Copy ratio array must have the same length as read counts array");
        Utils.validateArg(mask.length == numTargets, "Mask array must have the same length as read counts array");

        /* calculate total multiplicative factor */
        final double[] lambda = IntStream.range(0, numTargets)
                .mapToDouble(i -> multiplicativeBiases[i] * ploidies[i] * copyRatio[i])
                .toArray();

        /* mask out targets that have exactly zero Poisson factor to avoid division by zero */
        final int[] updatedMask = mask.clone();
        IntStream.range(0, numTargets)
                .filter(i -> lambda[i] == 0.0)
                .forEach(i -> { updatedMask[i] = 0; lambda[i] = 1; });
        Utils.validateArg(MathUtils.allMatch(updatedMask, m -> m == 0 || m == 1), "Target mask may only have 0 and 1 values");

        final long maskSum = Arrays.stream(updatedMask).mapToLong(m -> (long)m).sum();
        Utils.validateArg(maskSum > 0, "All of the targets are masked out, can not continue");

        final double readCountsSqrDivLambdaAve = IntStream.range(0, numTargets)
                .mapToDouble(i -> updatedMask[i] * readCounts[i] * readCounts[i] / lambda[i]).sum() / maskSum;
        final double lambdaAve = IntStream.range(0, numTargets)
                .mapToDouble(i -> updatedMask[i] * lambda[i]).sum() / maskSum;
        return (FastMath.sqrt(4 * readCountsSqrDivLambdaAve * lambdaAve + 1) - 1) / (2 * lambdaAve);
    }
}
