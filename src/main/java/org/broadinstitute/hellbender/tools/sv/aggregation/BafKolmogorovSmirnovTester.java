package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.hipparchus.stat.inference.KolmogorovSmirnovTest;

import java.util.*;

public class BafKolmogorovSmirnovTester {

    private final int minSnpCarriers;
    private final int minBafCount;
    private final KolmogorovSmirnovTest kolmogorovSmirnovTest;

    public BafKolmogorovSmirnovTester(final int minSnpCarriers, final int minBafCount, final long seed) {
        Utils.validateArg(minSnpCarriers > 0, "minSnpCarriers must be positive");
        Utils.validateArg(minBafCount > 0, "minBafCount must be positive");
        this.minSnpCarriers = minSnpCarriers;
        this.minBafCount = minBafCount;
        this.kolmogorovSmirnovTest = new KolmogorovSmirnovTest(seed);
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
        final Map<String, Object> attributes = new HashMap<>(record.getAttributes());
        final Double q = EvidenceStatUtils.probToQual(result.getP(), (byte) 99);
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

        double[] carrierBaf = frequencyFilteredEvidence.stream().filter(baf -> carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).toArray();
        if (carrierBaf.length < minBafCount) {
            return null;
        }
        double[] nullBaf = frequencyFilteredEvidence.stream().filter(baf -> !carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).toArray();
        if (nullBaf.length < minBafCount) {
            return null;
        }
        final double stat = kolmogorovSmirnovTest.kolmogorovSmirnovStatistic(carrierBaf, nullBaf);
        final double p = calculateP(stat, carrierBaf.length, nullBaf.length);

        return new KSTestResult(stat, p);
    }

    public double calculateP(final double stat, final int n, final int m) {
        final long lengthProduct = ((long) n) * m;
        if (lengthProduct < 10000) {
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

        public Double getStat() {
            return stat;
        }

        public Double getP() {
            return p;
        }
    }
}
