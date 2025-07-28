package org.broadinstitute.hellbender.tools.sv;

import java.util.*;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public class PermutationTTest {

    // Enum for alternative hypotheses
    public enum Alternative {
        TWO_SIDED("two.sided"),
        LESS("less"),
        GREATER("greater"),
        TWO_SIDED_ABS("two.sidedAbs");

        private final String value;
        Alternative(String value) { this.value = value; }
        public String getValue() { return value; }

        public static Alternative fromString(String text) {
            for (Alternative alt : Alternative.values()) {
                if (alt.value.equalsIgnoreCase(text)) {
                    return alt;
                }
            }
            throw new IllegalArgumentException("Unknown alternative: " + text);
        }
    }

    // Enum for test methods
    public enum Method {
        PCLT("pclt"),
        EXACT_MC("exact.mc"),
        EXACT_NETWORK("exact.network"),
        EXACT_CE("exact.ce");

        private final String value;
        Method(String value) { this.value = value; }
        public String getValue() { return value; }

        public static Method fromString(String text) {
            for (Method method : Method.values()) {
                if (method.value.equalsIgnoreCase(text)) {
                    return method;
                }
            }
            throw new IllegalArgumentException("Unknown method: " + text);
        }
    }

    // Control parameters class
    public static class PermControl {
        public double cm = 1e6;
        public int nmc = 999;
        public Integer seed = 1234321;
        public int digits = 12;
        public double pConfLevel = 0.99;
        public boolean setSEED = true;
        public String tsmethod = "central";

        public PermControl() {}

        public PermControl(double cm, int nmc, Integer seed, int digits,
                           double pConfLevel, boolean setSEED, String tsmethod) {
            this.cm = cm;
            this.nmc = nmc;
            this.seed = seed;
            this.digits = digits;
            this.pConfLevel = pConfLevel;
            this.setSEED = setSEED;
            this.tsmethod = tsmethod;
        }
    }

    // Result class to hold test results
    public static class PermTestResult {
        public double statistic;
        public double estimate;
        public Double parameter;
        public double pValue;
        public double nullValue;
        public Alternative alternative;
        public String method;
        public String dataName;
        public Map<String, Double> pValues;
        public Map<String, Double> pConfInt;
        public Integer nmc;
        public String testClass;

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("Permutation Test Result:\n");
            sb.append("Method: ").append(method).append("\n");
            sb.append("Data: ").append(dataName).append("\n");
            sb.append("Statistic: ").append(statistic).append("\n");
            sb.append("P-value: ").append(pValue).append("\n");
            sb.append("Alternative: ").append(alternative.getValue()).append("\n");
            sb.append("Estimate: ").append(estimate).append("\n");
            return sb.toString();
        }
    }

    // Internal result class for method outputs
    private static class MethodOutput {
        public Map<String, Double> pValues;
        public Double Z;
        public Map<String, Double> pConfInt;

        public MethodOutput(Map<String, Double> pValues, Double Z, Map<String, Double> pConfInt) {
            this.pValues = pValues;
            this.Z = Z;
            this.pConfInt = pConfInt;
        }
    }

    /**
     * Main permutation test function - equivalent to permTS.default in R
     */
    public static PermTestResult permTSDefault(double[] x, double[] y,
                                               Alternative alternative,
                                               Boolean exact,
                                               Method method,
                                               PermControl control) {

        // Handle default control
        if (control == null) control = new PermControl();

        // Handle deprecated two.sidedAbs alternative
        if (alternative == Alternative.TWO_SIDED_ABS) {
            System.out.println("Warning: alternative='two.sidedAbs' may be deprecated in the future, " +
                    "use alternative='two.sided' and control with tsmethod='abs'");
            alternative = Alternative.TWO_SIDED;
            control.tsmethod = "abs";
        }

        // Validate tsmethod for two-sided tests
        if (!("central".equals(control.tsmethod) || "abs".equals(control.tsmethod)) &&
                alternative == Alternative.TWO_SIDED) {
            throw new IllegalArgumentException("only tsmethod='central' and tsmethod='abs' allowed");
        }

        // Adjust alternative based on tsmethod
        if ("abs".equals(control.tsmethod) && alternative == Alternative.TWO_SIDED) {
            alternative = Alternative.TWO_SIDED_ABS;
        } else if ("central".equals(control.tsmethod) && alternative == Alternative.TWO_SIDED) {
            alternative = Alternative.TWO_SIDED;
        }

        // Validate input
        if (x == null || y == null) {
            throw new IllegalArgumentException("x and y must be numeric arrays");
        }

        // Combine data
        double[] W = new double[x.length + y.length];
        int[] Z = new int[x.length + y.length];
        System.arraycopy(x, 0, W, 0, x.length);
        System.arraycopy(y, 0, W, x.length, y.length);
        Arrays.fill(Z, 0, x.length, 1);
        Arrays.fill(Z, x.length, Z.length, 0);

        // Determine method if not specified
        if (method == null) {
            method = methodRuleTS1(W, Z, exact);
        }

        // Validate method
        if (!(method == Method.PCLT || method == Method.EXACT_MC ||
                method == Method.EXACT_NETWORK || method == Method.EXACT_CE)) {
            throw new IllegalArgumentException("method not one of: 'pclt', 'exact.mc', 'exact.network', 'exact.ce'");
        }

        // Execute appropriate method
        TwoSamplePCLT.Result mout;
        switch (method) {
            case PCLT:
                mout = TwoSamplePCLT.twosamplePclt(W, Z);
                break;
            default:
                throw new IllegalArgumentException("Unsupported method");
        }

        // Extract p-value based on alternative
        double pValue;
        switch (alternative) {
            case TWO_SIDED:
                pValue = mout.pValues.get("p.twosided");
                break;
            case GREATER:
                pValue = mout.pValues.get("p.gte");
                break;
            case LESS:
                pValue = mout.pValues.get("p.lte");
                break;
            case TWO_SIDED_ABS:
                pValue = mout.pValues.get("p.twosidedAbs");
                break;
            default:
                throw new IllegalArgumentException("Unknown alternative");
        }

        // Determine method description
        String methodDesc;
        switch (method) {
            case EXACT_NETWORK:
                methodDesc = "Exact Permutation Test (network algorithm)";
                break;
            case PCLT:
                methodDesc = "Permutation Test using Asymptotic Approximation";
                break;
            case EXACT_MC:
                methodDesc = "Exact Permutation Test Estimated by Monte Carlo";
                break;
            case EXACT_CE:
                methodDesc = "Exact Permutation Test (complete enumeration)";
                break;
            default:
                methodDesc = "Unknown method";
        }

        // Create data names (simplified version)
        String xname = "GROUP 1";
        String yname = "GROUP 2";
        String dataName = xname + " and " + yname;

        // Calculate estimates
        double xMean = Arrays.stream(x).average().orElse(0.0);
        double yMean = Arrays.stream(y).average().orElse(0.0);
        double estimate = xMean - yMean;

        // Build result
        PermTestResult result = new PermTestResult();
        result.statistic = mout.Z;
        result.estimate = estimate;
        result.parameter = null;
        result.pValue = pValue;
        result.nullValue = 0.0;
        result.alternative = alternative;
        result.method = methodDesc;
        result.dataName = dataName;
        result.pValues = mout.pValues;
        result.pConfInt = null;
        result.nmc = (method == Method.EXACT_MC) ? control.nmc : null;
        result.testClass = (method == Method.EXACT_MC) ? "mchtest" : "htest";

        return result;
    }

    /**
     * Simplified method rule - placeholder implementation
     */
    private static Method methodRuleTS1(double[] W, int[] Z, Boolean exact) {
        // Simple heuristic: use exact.mc for small samples, pclt for large
        if (exact != null && exact) {
            return Method.EXACT_MC;
        } else if (exact != null && !exact) {
            return Method.PCLT;
        } else {
            return W.length < 100 ? Method.EXACT_MC : Method.PCLT;
        }
    }

    /**
     * Placeholder implementations for the statistical methods
     * In a real implementation, these would contain the actual algorithms
     */
    private static MethodOutput twosamplePclt(double[] W, int[] Z) {
        // Placeholder for asymptotic approximation
        Map<String, Double> pValues = new HashMap<>();
        pValues.put("p.twosided", 0.05);
        pValues.put("p.gte", 0.025);
        pValues.put("p.lte", 0.975);
        pValues.put("p.twosidedAbs", 0.05);

        double statistic = calculateTestStatistic(W, Z);

        return new MethodOutput(pValues, statistic, new HashMap<>());
    }

    private static MethodOutput twosampleExactNetwork(double[] W, int[] Z, int digits) {
        // Placeholder for network algorithm
        return twosamplePclt(W, Z); // Simplified
    }

    private static MethodOutput twosampleExactCe(double[] W, int[] Z, double cm, int digits) {
        // Placeholder for complete enumeration
        return twosamplePclt(W, Z); // Simplified
    }

    private static MethodOutput twosampleExactMc(double[] W, int[] Z, Alternative alternative,
                                                 int nmc, Integer seed, int digits,
                                                 double pConfLevel, boolean setSEED) {
        // Placeholder for Monte Carlo exact test
        if (setSEED && seed != null) {
            Random random = new Random(seed);
        }

        // In real implementation, this would perform Monte Carlo permutation testing
        Map<String, Double> pValues = new HashMap<>();
        pValues.put("p.twosided", 0.05);
        pValues.put("p.gte", 0.025);
        pValues.put("p.lte", 0.975);
        pValues.put("p.twosidedAbs", 0.05);

        double statistic = calculateTestStatistic(W, Z);

        Map<String, Double> confInt = new HashMap<>();
        confInt.put("lower", -1.96);
        confInt.put("upper", 1.96);

        return new MethodOutput(pValues, statistic, confInt);
    }

    /**
     * Calculate test statistic (simplified implementation)
     */
    private static double calculateTestStatistic(double[] W, int[] Z) {
        double sum1 = 0.0, sum0 = 0.0;
        int n1 = 0, n0 = 0;

        for (int i = 0; i < W.length; i++) {
            if (Z[i] == 1) {
                sum1 += W[i];
                n1++;
            } else {
                sum0 += W[i];
                n0++;
            }
        }

        double mean1 = n1 > 0 ? sum1 / n1 : 0.0;
        double mean0 = n0 > 0 ? sum0 / n0 : 0.0;

        return mean1 - mean0;
    }

    // Convenience methods with default parameters
    public static PermTestResult permTSDefault(double[] x, double[] y) {
        return permTSDefault(x, y, Alternative.TWO_SIDED, null, null, null);
    }

    public static PermTestResult permTSDefault(double[] x, double[] y, Alternative alternative) {
        return permTSDefault(x, y, alternative, null, null, null);
    }

    public static PermTestResult permTSDefault(double[] x, double[] y, Alternative alternative, Method method) {
        return permTSDefault(x, y, alternative, null, method, null);
    }

    public class TwoSamplePCLT {

        public static class Result {
            public Map<String, Double> pValues;
            public double Z;

            public Result(Map<String, Double> pValues, double Z) {
                this.pValues = pValues;
                this.Z = Z;
            }

            @Override
            public String toString() {
                StringBuilder sb = new StringBuilder();
                sb.append("Z = ").append(Z).append("\n");
                sb.append("P-values:\n");
                for (Map.Entry<String, Double> entry : pValues.entrySet()) {
                    sb.append("  ").append(entry.getKey()).append(" = ").append(entry.getValue()).append("\n");
                }
                return sb.toString();
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
            double p_lte = normalCDF(Z);
            double p_gte = 1 - normalCDF(Z);
            double p_twosidedAbs = 1 - chiSquareCDF(Z * Z, 1);
            double p_twosided = Math.min(1.0, 2 * Math.min(p_lte, p_gte));

            // Create result map
            Map<String, Double> pValues = new LinkedHashMap<>();
            pValues.put("p.twosided", p_twosided);
            pValues.put("p.lte", p_lte);
            pValues.put("p.gte", p_gte);
            pValues.put("p.twosidedAbs", p_twosidedAbs);

            return new Result(pValues, Z);
        }

        // Approximation of standard normal CDF using error function
        private static double normalCDF(double x) {
            return 0.5 * (1 + erf(x / Math.sqrt(2)));
        }

        // Error function approximation (Abramowitz and Stegun)
        private static double erf(double x) {
            double a1 =  0.254829592;
            double a2 = -0.284496736;
            double a3 =  1.421413741;
            double a4 = -1.453152027;
            double a5 =  1.061405429;
            double p  =  0.3275911;

            int sign = x < 0 ? -1 : 1;
            x = Math.abs(x);

            double t = 1.0 / (1.0 + p * x);
            double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

            return sign * y;
        }

        // Chi-square CDF approximation for df=1
        private static double chiSquareCDF(double x, int df) {
            if (df == 1) {
                // For df=1, chi-square CDF = 2*Φ(sqrt(x)) - 1, where Φ is standard normal CDF
                return 2 * normalCDF(Math.sqrt(x)) - 1;
            }
            throw new IllegalArgumentException("Only df=1 is supported in this implementation");
        }
    }
}