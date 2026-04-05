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

    private static final double MAX_QUAL = 999;

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
        final double stat = kolmogorovSmirnovTest.kolmogorovSmirnovStatistic(carrierBaf, nullBaf);
        final double p = calculateP(stat, carrierBaf.length, nullBaf.length);

        return new KSTestResult(stat, p);
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
