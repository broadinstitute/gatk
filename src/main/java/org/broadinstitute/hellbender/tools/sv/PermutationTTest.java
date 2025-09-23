package org.broadinstitute.hellbender.tools.sv;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.IntStream;

public class PermutationTTest {

    private static final NormalDistribution STD_NORMAL_DISTRIBUTION = new NormalDistribution();
    private static final ChiSquaredDistribution CHI_SQUARED_DISTRIBUTION = new ChiSquaredDistribution(1);


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
    public static PermTestResult permTSDefault(double[] x, double[] y, Alternative alternative) {

        // Validate input
        Utils.nonNull(x);
        Utils.nonNull(y);

        // Combine data
        double[] W = new double[x.length + y.length];
        int[] Z = new int[x.length + y.length];
        System.arraycopy(x, 0, W, 0, x.length);
        System.arraycopy(y, 0, W, x.length, y.length);
        Arrays.fill(Z, 0, x.length, 1);
        Arrays.fill(Z, x.length, Z.length, 0);

        // Execute appropriate method
        TwoSamplePCLT.Result mout = TwoSamplePCLT.twosamplePclt(W, Z);

        // Extract p-value based on alternative
        double pValue;
        switch (alternative) {
            case TWO_SIDED:
                pValue = mout.pValues.get(Alternative.TWO_SIDED);
                break;
            case GREATER:
                pValue = mout.pValues.get(Alternative.GREATER);
                break;
            case LESS:
                pValue = mout.pValues.get(Alternative.LESS);
                break;
            default:
                throw new IllegalArgumentException("Unknown alternative");
        }

        return new PermTestResult(mout.Z, pValue, alternative);
    }

    public class TwoSamplePCLT {

        public static class Result {
            public Map<Alternative, Double> pValues;
            public double Z;

            public Result(Map<Alternative, Double> pValues, double Z) {
                this.pValues = pValues;
                this.Z = Z;
            }
        }

        public static Result twosamplePclt(double[] scores, int[] group) {
            int n = scores.length;

            // Create frequency table equivalent
            Map<Integer, Map<Double, Integer>> tab = new HashMap<>();
            for (int i = 0; i < n; i++) {
                tab.computeIfAbsent(group[i], k -> new HashMap<>())
                        .merge(scores[i], 1, Integer::sum);
            }

            // Find the groups and determine which is "group 1"
            Set<Integer> uniqueGroups = new HashSet<>();
            for (int g : group) {
                uniqueGroups.add(g);
            }
            Integer[] groupArray = uniqueGroups.toArray(new Integer[0]);
            Arrays.sort(groupArray);

            // The second group in sorted order becomes "Grp1" (equivalent to dimnames(tab)[[1]][2])
            int Grp1 = groupArray.length > 1 ? groupArray[1] : groupArray[0];

            // Create binary group indicator (grp)
            int[] grp = new int[n];
            for (int i = 0; i < n; i++) {
                grp[i] = group[i] == Grp1 ? 1 : 0;
            }

            // Calculate T0 (sum of scores * grp)
            double T0 = IntStream.range(0, n)
                    .mapToDouble(i -> scores[i] * grp[i])
                    .sum();

            // Calculate means
            double meanScores = Arrays.stream(scores).average().orElse(0.0);
            double meanGrp = Arrays.stream(grp).average().orElse(0.0);

            // Calculate SSE for scores
            double SSE_scores = Arrays.stream(scores)
                    .map(x -> Math.pow(x - meanScores, 2))
                    .sum();

            // Calculate SSE for group
            double SSE_grp = Arrays.stream(grp)
                    .mapToDouble(x -> Math.pow(x - meanGrp, 2))
                    .sum();

            // Calculate Z statistic
            double Z = Math.sqrt(n - 1) * (T0 - n * meanScores * meanGrp) /
                    Math.sqrt(SSE_scores * SSE_grp);

            // Calculate p-values using standard normal distribution
            double p_lte = STD_NORMAL_DISTRIBUTION.cumulativeProbability(Z);
            double p_gte = 1 - STD_NORMAL_DISTRIBUTION.cumulativeProbability(Z);
            double p_twosided = Math.min(1.0, 2 * Math.min(p_lte, p_gte));

            // Create result map
            Map<Alternative, Double> pValues = new LinkedHashMap<>();
            pValues.put(Alternative.TWO_SIDED, p_twosided);
            pValues.put(Alternative.LESS, p_lte);
            pValues.put(Alternative.GREATER, p_gte);

            return new Result(pValues, Z);
        }
    }
}