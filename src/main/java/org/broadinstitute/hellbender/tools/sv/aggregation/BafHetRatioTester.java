package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class BafHetRatioTester {

    private static final Median MEDIAN = new Median();

    private final double pSnp;
    private final double pMaxHomozygous;

    public BafHetRatioTester(final double pSnp, final double pMaxHomozygous) {
        Utils.validateArg(pSnp > 0 && pSnp <= 1, "pSnp must be a probability on (0, 1]");
        Utils.validateArg(pMaxHomozygous > 0 && pMaxHomozygous < 1, "pMaxHomozygous must be a probability on (0, 1)");
        this.pSnp = pSnp;
        this.pMaxHomozygous = pMaxHomozygous;
    }

    public Double test(final SVCallRecord record, final List<BafEvidence> evidence, final Set<String> allSamples,
                       final Set<String> carrierSamples, final int flankSize) {
        if (!record.isSimpleCNV() || evidence == null || evidence.isEmpty()) {
            return null;
        }

        final Map<String, HetSnpStats> sampleStats = allSamples.stream().collect(Collectors.toMap(s -> s, s -> new HetSnpStats()));
        for (final BafEvidence baf : evidence) {
            if (baf.getStart() < record.getPositionA()) {
                sampleStats.get(baf.getSample()).upstreamHets++;
            } else if (baf.getStart() >= record.getPositionB()) {
                sampleStats.get(baf.getSample()).downstreamHets++;
            } else {
                sampleStats.get(baf.getSample()).containedHets++;
            }
        }

        return calculate(record.getLength(), sampleStats, carrierSamples, flankSize);
    }

    public SVCallRecord applyToRecord(final SVCallRecord record, final Double result) {
        Utils.nonNull(record);
        if (result == null) {
            return record;
        }
        final Map<String, Object> attributes = new HashMap<>(record.getAttributes());
        attributes.put(GATKSVVCFConstants.BAF_HET_RATIO_ATTRIBUTE, result);
        return SVCallRecordUtils.copyCallWithNewAttributes(record, attributes);
    }

    private Double calculate(final int length,
                             final Map<String, HetSnpStats> sampleStats,
                             final Set<String> carrierSamples,
                             final int flankSize) {
        //final double pSnp = 0.001; //Math.min(50. / length, 0.0005);
        //int totalInnerCount = 0;
        final List<Double> nullLogRatios = new ArrayList<>();
        final List<Double> carrierLogRatios = new ArrayList<>();

        final BinomialDistribution binomialDistributionFlank = new BinomialDistribution(flankSize, pSnp);
        final BinomialDistribution binomialDistributionInner = new BinomialDistribution(length, pSnp);

        for (final Map.Entry<String, HetSnpStats> entry : sampleStats.entrySet()) {
            final String sample = entry.getKey();
            final HetSnpStats stats = entry.getValue();
            final double pFlank = binomialDistributionFlank.cumulativeProbability(Math.min(stats.upstreamHets, stats.downstreamHets));
            final double pInner = binomialDistributionInner.cumulativeProbability(stats.containedHets);
            if (!(pInner < pMaxHomozygous && pFlank < pMaxHomozygous)) {
                stats.logRatio = Math.log(stats.containedHets + 1.);
                if (carrierSamples.contains(sample)) {
                    carrierLogRatios.add(stats.logRatio);
                } else {
                    nullLogRatios.add(stats.logRatio);
                }
            }
            //totalInnerCount += stats.innerHetCount;
        }

        if (carrierLogRatios.isEmpty() || nullLogRatios.isEmpty()) {
            return null;
        }
        final double medianCarrier = MEDIAN.evaluate(carrierLogRatios.stream().mapToDouble(Double::doubleValue).toArray());
        final double medianNull = MEDIAN.evaluate(nullLogRatios.stream().mapToDouble(Double::doubleValue).toArray());
        return medianNull - medianCarrier;

        /*
        if (nullRatios.size() <= 10 || totalInnerCount < 10) {
            return null;
        }
        final double minNullRatio = nullRatios.stream().mapToDouble(d -> d).min().getAsDouble();
        final double maxNullRatio = nullRatios.stream().mapToDouble(d -> d).max().getAsDouble();
        if (maxNullRatio - minNullRatio < 0.0001) {
            return null;
        }
        final List<SampleStats> carrierStats = carrierSamples.stream().map(sampleStats::get).filter(s -> s.deletionRatio != null).collect(Collectors.toList());
        if (carrierStats.isEmpty() || carrierStats.size() > nullRatios.size()) {
            return null;
        }
        final double[] carrierRatiosArr = carrierSamples.stream().map(sampleStats::get).filter(s -> s.deletionRatio != null).mapToDouble(s -> s.deletionRatio).toArray();
        return median.evaluate(carrierRatiosArr);
         */
    }

    private static final class HetSnpStats {
        public int upstreamHets = 0;
        public int containedHets = 0;
        public int downstreamHets = 0;
        public Double logRatio = null;
    }

    /*
    protected boolean isRegionOfHomozygosity(final int beforeHetCount, final int innerHetCount, final int afterHetCount, final double length, final double threshold) {
        return innerHetCount < threshold * length &&
                (beforeHetCount  < threshold * length || afterHetCount < threshold * length);
    }

    protected double calculateDeletionRatio(final int beforeHetCount, final int innerHetCount, final int afterHetCount, final double length, final double threshold) {
        final int flankHetCount = Math.min(beforeHetCount, afterHetCount);
        return Math.log10(innerHetCount + length * threshold) - Math.log10(flankHetCount + length * threshold);
    }
     */
}
