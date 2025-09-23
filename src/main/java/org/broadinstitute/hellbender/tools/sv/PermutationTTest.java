package org.broadinstitute.hellbender.tools.sv;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.IntStream;

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

        return twosamplePclt(W, Z, alternative);
    }

    public static PermTestResult twosamplePclt(final double[] scores, final int[] groups, final Alternative alternative) {
        final int n = scores.length;

        // Find the groups and determine which is "group 1"
        final Set<Integer> uniqueGroups = new HashSet<>();
        for (int g : groups) {
            uniqueGroups.add(g);
        }
        final Integer[] groupArray = uniqueGroups.toArray(new Integer[0]);
        Arrays.sort(groupArray);

        // The second group in sorted order becomes "Grp1"
        final int group1 = groupArray.length > 1 ? groupArray[1] : groupArray[0];

        // Create binary group indicator (grp)
        final int[] groupIndicators = new int[n];
        for (int i = 0; i < n; i++) {
            groupIndicators[i] = groups[i] == group1 ? 1 : 0;
        }

        // Calculate T0 (sum of scores * grp)
        final double t0 = IntStream.range(0, n)
                .mapToDouble(i -> scores[i] * groupIndicators[i])
                .sum();

        // Calculate means
        final double meanScores = Arrays.stream(scores).average().orElse(0.0);
        final double meanGrp = Arrays.stream(groupIndicators).average().orElse(0.0);

        // Calculate SSE for scores
        final double sseScores = Arrays.stream(scores)
                .map(x -> Math.pow(x - meanScores, 2))
                .sum();

        // Calculate SSE for group
        final double sseGroup = Arrays.stream(groupIndicators)
                .mapToDouble(x -> Math.pow(x - meanGrp, 2))
                .sum();

        // Calculate Z statistic
        final double Z = Math.sqrt(n - 1) * (t0 - n * meanScores * meanGrp) /
                Math.sqrt(sseScores * sseGroup);

        // Calculate p-values using standard normal distribution
        final double pLess = STD_NORMAL_DISTRIBUTION.cumulativeProbability(Z);
        final double pGreater = 1 - STD_NORMAL_DISTRIBUTION.cumulativeProbability(Z);

        final double pValue = switch (alternative) {
            case TWO_SIDED -> Math.min(1.0, 2 * Math.min(pLess, pGreater));
            case GREATER -> pGreater;
            case LESS -> pLess;
            default -> throw new IllegalArgumentException("Unknown alternative");
        };
        return new PermTestResult(Z, pValue, alternative);
    }
}