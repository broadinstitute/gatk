package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
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
     * Estimate mean read depth density $\rho$ according to a Poisson model for the distribution of read counts.
     *
     * The Poisson model is:
     *
     *     n_t \sim \mathrm{Poisson}(\rho l_t c_t m_t p_t)
     *
     * where $n_t$ is the read count, $l_t$ is target length, $c_t$ is copy ratio, $p_t$ is ploidy, and $m_t$ is
     * the multiplicative bias (the estimation method is described in CNV-methods.pdf).
     *
     * Note: the mask array is 0 or 1, and targets with 0 mask will be effectively dropped out.
     *
     * Note: the product $l_t c_t m_t p_t$ must be positive on all targets. Otherwise, they will be automatically
     * masked out (required for stability).
     *
     * @param readCounts read counts
     * @param targetLengths target lengths
     * @param multiplicativeBiases estimate of total multiplicative bias
     * @param ploidies nominal germline ploidies of targets
     * @param copyRatio estimate of copy ratio
     * @param mask target masks (only 0 and 1 is expected)
     * @return read depth per base per homolog
     */
    public static double estimateReadDepthDensityFromPoissonModel(final double[] readCounts,
                                                                  final int[] targetLengths,
                                                                  final double[] multiplicativeBiases,
                                                                  final int[] ploidies,
                                                                  final double[] copyRatio,
                                                                  final int[] mask) {
        Utils.nonNull(readCounts, "Read count array must be non-null");
        Utils.nonNull(targetLengths, "Target length array must be non-null");
        Utils.nonNull(multiplicativeBiases, "Multiplicative bias array must be non-null");
        Utils.nonNull(ploidies, "Ploidy array must be non-null");
        Utils.nonNull(copyRatio, "Copy ratio array must be non-null");
        Utils.nonNull(mask, "Mask array must be non-null");
        final int numTargets = readCounts.length;
        if (numTargets == 0) {
            throw new UserException.BadInput("At least one target is required");
        }
        if (targetLengths.length != numTargets) {
            throw new UserException.BadInput("Target lengths array must have the same length as read counts array");
        }
        if (multiplicativeBiases.length != numTargets) {
            throw new UserException.BadInput("Multiplicative bias array must have the same length as read counts array");
        }
        if (ploidies.length != numTargets) {
            throw new UserException.BadInput("Ploidies array must have the same length as read counts array");
        }
        if (copyRatio.length != numTargets) {
            throw new UserException.BadInput("Copy ratio array must have the same length as read counts array");
        }
        if (mask.length != numTargets) {
            throw new UserException.BadInput("Mask array must have the same length as read counts array");
        }

        /* calculate total multiplicative factor */
        final double[] lambda = IntStream.range(0, numTargets)
                .mapToDouble(i -> targetLengths[i] * multiplicativeBiases[i] * ploidies[i] * copyRatio[i])
                .toArray();

        /* mask out targets that have exactly zero Poisson factor to avoid division by zero */
        final int[] updatedMask = mask.clone();
        IntStream.range(0, numTargets)
                .filter(i -> lambda[i] == 0.0)
                .forEach(i -> { updatedMask[i] = 0; lambda[i] = 1; });
        if (Arrays.stream(updatedMask).filter(m -> m != 0 && m != 1).count() > 0) {
            throw new UserException.BadInput("Target mask may only have 0 and 1 values");
        }
        final long maskSum = Arrays.stream(updatedMask).mapToLong(m -> (long)m).sum();
        if (maskSum == 0) {
            throw new UserException.BadInput("All of the targets are masked out, can not continue");
        }

        final double readCountsSqrDivLambdaAve = IntStream.range(0, numTargets)
                .mapToDouble(i -> updatedMask[i] * readCounts[i] * readCounts[i] / lambda[i]).sum() / maskSum;
        final double lambdaAve = IntStream.range(0, numTargets)
                .mapToDouble(i -> updatedMask[i] * lambda[i]).sum() / maskSum;
        return (FastMath.sqrt(4 * readCountsSqrDivLambdaAve * lambdaAve + 1) - 1) / (2 * lambdaAve);
    }
}
