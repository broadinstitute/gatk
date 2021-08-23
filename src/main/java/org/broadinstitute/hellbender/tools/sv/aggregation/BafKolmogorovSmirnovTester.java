package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public class BafKolmogorovSmirnovTester {

    private final KolmogorovSmirnovTest KS_TEST =  new KolmogorovSmirnovTest();

    private final int minSnpCarriers;
    private final int minBafCount;

    public BafKolmogorovSmirnovTester(final int minSnpCarriers, final int minBafCount) {
        Utils.validateArg(minSnpCarriers > 0, "minSnpCarriers must be positive");
        Utils.validateArg(minBafCount > 0, "minBafCount must be positive");
        this.minSnpCarriers = minSnpCarriers;
        this.minBafCount = minBafCount;
    }

    public KSTestResult test(final SVCallRecord record, final List<BafEvidence> evidence, final Set<String> carrierSamples) {
        if (!record.isSimpleCNV() || evidence == null || evidence.isEmpty()) {
            return null;
        }

        final List<BafEvidence> bafList = new ArrayList<>(evidence.size());
        for (final BafEvidence baf : evidence) {
            if (baf.getStart() >= record.getPositionA() && baf.getStart() < record.getPositionB()) {
                bafList.add(baf);
            }
        }
        return calculate(bafList, carrierSamples);
    }

    public SVCallRecord applyToRecord(final SVCallRecord record, final KSTestResult result) {
        Utils.nonNull(record);
        if (result == null) {
            return record;
        }
        final Map<String, Object> attributes = new HashMap<>();
        final Integer q = EvidenceStatUtils.probToQual(result.getP(), (byte) 99);
        attributes.put(GATKSVVCFConstants.BAF_KS_STAT_ATTRIBUTE, result.getStat());
        attributes.put(GATKSVVCFConstants.BAF_KS_Q_ATTRIBUTE, q);
        return SVCallRecordUtils.copyCallWithNewAttributes(record, attributes);
    }

    private KSTestResult calculate(final List<BafEvidence> evidence, final Set<String> carrierSamples) {
        final List<BafEvidence> frequencyFilteredEvidence = new ArrayList<>();
        final Iterator<BafEvidence> iter = evidence.iterator();
        final List<BafEvidence> buffer = new ArrayList<>();
        int pos = -1;
        while (iter.hasNext()) {
            final BafEvidence baf = iter.next();
            if (baf.getStart() != pos) {
                if (buffer.size() >= minSnpCarriers) {
                    frequencyFilteredEvidence.addAll(buffer);
                }
                buffer.clear();
                pos = baf.getStart();
            }
            buffer.add(baf);
        }
        if (buffer.size() >= minSnpCarriers) {
            frequencyFilteredEvidence.addAll(buffer);
        }

        double[] carrierBaf = frequencyFilteredEvidence.stream().filter(baf -> carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).map(x -> Math.min(x, 1.0 - x)).toArray();
        if (carrierBaf.length < minBafCount) {
            return null;
        }
        double[] nullBaf = frequencyFilteredEvidence.stream().filter(baf -> !carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).map(x -> Math.min(x, 1.0 - x)).toArray();
        if (nullBaf.length < minBafCount) {
            return null;
        }
        final double stat = KS_TEST.kolmogorovSmirnovStatistic(carrierBaf, nullBaf);
        final double p = calculateP(stat, carrierBaf.length, nullBaf.length);
        return new KSTestResult(stat, p);
        //final double p = KS_TEST.monteCarloP(stat, carrierBaf.length, nullBaf.length, true, 10000);
        //final double p = kolmogorovSmirnovTest(stat, carrierBaf, nullBaf);
        //return BafResult.createDuplicationResult(p, stat);
        //final MannWhitneyU test = new MannWhitneyU();
        //final MannWhitneyU.Result stat = test.test(carrierBaf, nullBaf, MannWhitneyU.TestType.TWO_SIDED);
        //return BafResult.createDuplicationResult(stat.getP());
        //final MannWhitneyUTest test = new MannWhitneyUTest();
        //final double u = (long) carrierBaf.length * nullBaf.length - test.mannWhitneyU(nullBaf, carrierBaf);

        //final double p = calculateAsymptoticPValue(u, carrierBaf.length, nullBaf.length);
        //return BafResult.createDuplicationResult(p, u);
        /*
        final double nullMedian = MEDIAN.evaluate(nullBaf);
        final double carrierMedian = MEDIAN.evaluate(carrierBaf);
        final double nullStd = STDEV.evaluate(nullBaf);
        final double carrierStd = STDEV.evaluate(carrierBaf);
        return (nullMedian - carrierMedian) / Math.sqrt(nullStd * nullStd + carrierStd * carrierStd);
        */
    }

    public double calculateP(final double stat, final int n, final int m) {
        final long lengthProduct = (long) n * m;
        if (lengthProduct < 100) {
            return KS_TEST.exactP(stat, n, m, true);
        }
        return KS_TEST.approximateP(stat, n, m);
    }

    public static final class KSTestResult {
        private final Double stat;
        private final Double p;

        private KSTestResult(final Double stat, final Double p) {
            this.stat = stat;
            this.p = p;
        }

        public Double getStat() {
            return stat;
        }

        public Double getP() {
            return p;
        }
    }

    /*
    private static double[] shrinkArray(final double[] arr, final int size) {
        if (arr.length <= size) {
            return arr;
        }
        final double[] sorted = Arrays.copyOf(arr, arr.length);
        Arrays.sort(sorted);
        final int k = arr.length / size;
        int j = k / 2;
        final double[] out = new double[size];
        for (int i = 0; i < size; i++) {
            out[i] = sorted[j];
            j += k;
        }
        return out;
    }

    private static double[] downsampleArray(final double[] arr, final int size) {
        for (int i = arr.length - 1; i >= arr.length - size; i--) {
            final double tmp = arr[i];
            final int j = RAND.nextInt(arr.length);
            arr[i] = arr[j];
            arr[j] = tmp;
        }
        return Arrays.copyOf(arr, size);
    }
     */
}
