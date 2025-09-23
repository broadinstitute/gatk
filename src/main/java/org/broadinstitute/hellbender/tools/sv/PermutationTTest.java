package org.broadinstitute.hellbender.tools.sv;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public class PermutationTTest {

    protected static final NormalDistribution STD_NORMAL_DISTRIBUTION = new NormalDistribution();

    // Enum for alternative hypotheses
    public enum Alternative {
        TWO_SIDED,
        LESS,
        GREATER
    }

    // Result class to hold test results
    public record PermTestResult(double statistic, double pValue, Alternative alternative) {
    }

    /**
     * Main permutation test function - equivalent to permTS.default in R
     */
    public static PermTestResult test(final double[] x, final double[] y, final Alternative alternative) {
        Utils.nonNull(x);
        Utils.nonNull(y);
        Utils.nonNull(alternative);

        // Combine data
        final double[] W = new double[x.length + y.length];
        final int[] Z = new int[x.length + y.length];
        System.arraycopy(x, 0, W, 0, x.length);
        System.arraycopy(y, 0, W, x.length, y.length);
        Arrays.fill(Z, 0, x.length, 1);
        Arrays.fill(Z, x.length, Z.length, 0);

        final int n = W.length;
        final int nFirstGroup = x.length;

        // Calculate T0 (sum of scores * grp)
        final double t0 = MathUtils.sum(W, 0, nFirstGroup);

        // Calculate means
        final double meanScores = Arrays.stream(W).average().orElse(0.0);
        final double meanGrp = n == 0 ? 0 : nFirstGroup / (double) n;

        // Calculate SSE for scores
        final double sseScores = Arrays.stream(W).map(a -> Math.pow(a - meanScores, 2)).sum();

        // Calculate SSE for group
        final double sseGroup = nFirstGroup * Math.pow(1 - meanGrp, 2) + (n - nFirstGroup) * Math.pow(-meanGrp, 2);

        // Calculate Z statistic
        final double statZ = Math.sqrt(n - 1) * (t0 - n * meanScores * meanGrp) /
                Math.sqrt(sseScores * sseGroup);

        // Calculate p-values using standard normal distribution
        final double pLess = STD_NORMAL_DISTRIBUTION.cumulativeProbability(statZ);
        final double pGreater = 1 - STD_NORMAL_DISTRIBUTION.cumulativeProbability(statZ);

        final double pValue = switch (alternative) {
            case TWO_SIDED -> Math.min(1.0, 2 * Math.min(pLess, pGreater));
            case GREATER -> pGreater;
            case LESS -> pLess;
            default -> throw new IllegalArgumentException("Unknown alternative");
        };
        return new PermTestResult(statZ, pValue, alternative);
    }
}