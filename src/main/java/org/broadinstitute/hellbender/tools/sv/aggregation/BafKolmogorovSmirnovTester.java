package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.hipparchus.stat.inference.KolmogorovSmirnovTest;

import java.util.*;

/**
 * Calculates metrics for detecting copy number variants from B-allele (BAF) evidence. Het SNP calls are aggregated over
 * the variant interval and BAFs in carrier samples are compared to those in controls using a Kolmogorov-Smirnov test.
 *
 * <p>This matches the v1.1 svtk KS2sample implementation: no per-locus frequency filtering is applied.</p>
 */

public class BafKolmogorovSmirnovTester {

    // Cap at 50 to approximate scipy.stats.ks_2samp effective precision (~10^-300 -> -log10 ~300)
    // and avoid distorting the BAF1 RF feature space with extreme outliers
    private static final double MAX_QUAL = 50;

    private final int minBafCount;
    private final KolmogorovSmirnovTest kolmogorovSmirnovTest;

    // For N carrier BAFs and M control BAFs, use KS monte carlo approximation if N*M exceeds this threshold
    private final static long KS_APPROX_MIN_LENGTH_PRODUCT = 10000L;

    /**
     * @param minBafCount       min BAFs needed in both cases and controls
     * @param seed              PRNG seed for KS test approximation
     */
    public BafKolmogorovSmirnovTester(final int minBafCount, final long seed) {
        Utils.validateArg(minBafCount > 0, "minBafCount must be positive");
        this.minBafCount = minBafCount;
        this.kolmogorovSmirnovTest = new KolmogorovSmirnovTest(seed);
    }

    /**
     * @param record            query record
     * @param evidence          BAF evidence over the given variant
     * @param carrierSamples    samples putatively carrying the variant
     * @return                  KS statistic and p-value, or null if the test fails
     */
    public KSTestResult test(final SVCallRecord record, final List<BafEvidence> evidence, final Set<String> carrierSamples) {
        if (!record.isSimpleCNV() || evidence == null || evidence.isEmpty()) {
            return null;
        }
        // subset to variant interval
        final List<BafEvidence> bafList = new ArrayList<>(evidence.size());
        for (final BafEvidence baf : evidence) {
            if (baf.getStart() >= record.getPositionA() && baf.getStart() < record.getPositionB()) {
                bafList.add(baf);
            }
        }
        return calculate(bafList, carrierSamples);
    }

    /**
     * Adds KS stat and -log10(p-value) to record attributes. Returns a shallow copy of the new record.
     */
    public SVCallRecord applyToRecord(final SVCallRecord record, final KSTestResult result) {
        Utils.nonNull(record);
        if (result == null) {
            return record;
        }
        final Map<String, Object> attributes = new HashMap<>(record.getAttributes());
        final double q = Math.min(-Math.log10(result.getP()), MAX_QUAL);
        attributes.put(GATKSVVCFConstants.BAF_KS_STAT_ATTRIBUTE, result.getStat());
        attributes.put(GATKSVVCFConstants.BAF_KS_Q_ATTRIBUTE, q);
        return SVCallRecordUtils.copyCallWithNewAttributes(record, attributes);
    }


    private KSTestResult calculate(final List<BafEvidence> evidence, final Set<String> carrierSamples) {
        // No per-locus frequency filtering (matching v1.1)
        double[] carrierBaf = evidence.stream().filter(baf -> carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).toArray();
        if (carrierBaf.length < minBafCount) {
            return null;
        }
        double[] nullBaf = evidence.stream().filter(baf -> !carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).toArray();
        if (nullBaf.length < minBafCount) {
            return null;
        }
        final double stat;
        final double p;
        if (carrierBaf.length >= 2 && nullBaf.length >= 2) {
            // Standard hipparchus two-sample KS test
            stat = kolmogorovSmirnovTest.kolmogorovSmirnovStatistic(carrierBaf, nullBaf);
            p = calculateP(stat, carrierBaf.length, nullBaf.length);
        } else {
            // At least one sample has exactly 1 element; hipparchus requires n >= 2,
            // so compute the KS statistic and exact p-value manually to match scipy ks_2samp.
            final double[] result = singletonKsTest(carrierBaf, nullBaf);
            stat = result[0];
            p = result[1];
        }

        return new KSTestResult(stat, p);
    }

    /**
     * Computes the two-sample KS statistic and exact p-value when at least one sample has size 1.
     * Matches scipy.stats.ks_2samp behavior: the statistic is the standard two-sample KS stat
     * (well-defined for any n >= 1), and the exact p-value is computed by enumerating all (n+m choose 1)
     * possible single-element assignments from the combined sample.
     *
     * @param sample1 first sample (length >= 1)
     * @param sample2 second sample (length >= 1)
     * @return double[]{statistic, pValue}
     */
    static double[] singletonKsTest(final double[] sample1, final double[] sample2) {
        // Normalize so the singleton is always singletonArr, the other is otherArr
        // KS test is symmetric so this doesn't affect results
        final double[] singletonArr;
        final double[] otherArr;
        if (sample1.length <= sample2.length) {
            singletonArr = sample1;
            otherArr = sample2;
        } else {
            singletonArr = sample2;
            otherArr = sample1;
        }
        Utils.validate(singletonArr.length == 1,
                "singletonKsTest requires at least one sample of size 1, got " + singletonArr.length);

        final double x = singletonArr[0];
        final int m = otherArr.length;

        // KS statistic: max|F_singleton(t) - F_other(t)|
        // F_singleton is a step function: 0 for t < x, 1 for t >= x
        // The max difference occurs either just below x (dMinus) or at x (dPlus)
        int countLess = 0;
        int countLessOrEqual = 0;
        for (final double v : otherArr) {
            if (v < x) {
                countLess++;
            }
            if (v <= x) {
                countLessOrEqual++;
            }
        }
        final double dMinus = (double) countLess / m;
        final double dPlus = 1.0 - (double) countLessOrEqual / m;
        final double stat = Math.max(dMinus, dPlus);

        // Exact p-value: enumerate all (m+1) possible assignments of which element is the singleton
        // from the combined sample. Count the fraction that produces a stat >= observed.
        final double[] combined = new double[m + 1];
        System.arraycopy(otherArr, 0, combined, 0, m);
        combined[m] = x;
        Arrays.sort(combined);

        int countGe = 0;
        for (int pick = 0; pick <= m; pick++) {
            int kLess = 0;
            int kLessOrEqual = 0;
            for (int j = 0; j <= m; j++) {
                if (j == pick) {
                    continue;
                }
                if (combined[j] < combined[pick]) {
                    kLess++;
                }
                if (combined[j] <= combined[pick]) {
                    kLessOrEqual++;
                }
            }
            final double pickDMinus = (double) kLess / m;
            final double pickDPlus = 1.0 - (double) kLessOrEqual / m;
            final double pickStat = Math.max(pickDMinus, pickDPlus);
            if (pickStat >= stat - 1e-10) {
                countGe++;
            }
        }

        final double pValue = (double) countGe / (m + 1);
        return new double[]{stat, pValue};
    }

    public double calculateP(final double stat, final int n, final int m) {
        final long lengthProduct = ((long) n) * m;
        if (lengthProduct < KS_APPROX_MIN_LENGTH_PRODUCT) {
            return kolmogorovSmirnovTest.exactP(stat, n, m, false);
        }
        return kolmogorovSmirnovTest.approximateP(stat, n, m);
    }

    public static final class KSTestResult {
        private final Double stat;
        private final Double p;

        public KSTestResult(final Double stat, final Double p) {
            this.stat = stat;
            this.p = p;
        }

        /**
         * @return KS test statistic
         */
        public Double getStat() {
            return stat;
        }

        /**
         * @return KS test p-value
         */
        public Double getP() {
            return p;
        }
    }
}
