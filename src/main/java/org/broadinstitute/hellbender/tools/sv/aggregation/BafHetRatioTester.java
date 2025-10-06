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

/**
 * Measures quality of a CNV with a metric of the ratio of median heterozygous SNPs in carriers to that in controls.
 */

public class BafHetRatioTester {

    private static final Median MEDIAN = new Median();

    private final double pSnp;
    private final double pMaxHomozygous;

    /**
     * @param pSnp              prior probability of SNP at any locus
     * @param pMaxHomozygous    probability threshold for detecting regions of homozygosity
     */
    public BafHetRatioTester(final double pSnp, final double pMaxHomozygous) {
        Utils.validateArg(pSnp > 0 && pSnp <= 1, "pSnp must be a probability on (0, 1]");
        Utils.validateArg(pMaxHomozygous > 0 && pMaxHomozygous < 1, "pMaxHomozygous must be a probability on (0, 1)");
        this.pSnp = pSnp;
        this.pMaxHomozygous = pMaxHomozygous;
    }

    /**
     * @param record            query DEL/DUP record
     * @param evidence          BAF evidence associated with this record, including flanking regions
     * @param allSamples        all samples in the record
     * @param carrierSamples    carrier samples from the record
     * @param flankSize         flank size, for detecting ROH
     * @return                  log ratio of median het SNP count in carriers to controls
     */
    public Double test(final SVCallRecord record, final List<BafEvidence> evidence, final Set<String> allSamples,
                       final Set<String> carrierSamples, final int flankSize) {
        Utils.nonNull(record);
        Utils.nonNull(allSamples);
        Utils.nonNull(carrierSamples);
        Utils.validateArg(flankSize >= 0, "flankSize cannot be negative");
        Utils.validateArg(allSamples.size() >= carrierSamples.size() && allSamples.containsAll(carrierSamples), "Sample set must contain all carrier samples");

        if (!record.isSimpleCNV() || evidence == null || evidence.isEmpty()) {
            return null;
        }

        // Count het SNPs in upstream/downstream flanks and across the variant itself
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

    /**
     * Annotates record with het SNP ratio
     */
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
        final List<Double> nullLogCounts = new ArrayList<>();
        final List<Double> carrierLogCounts = new ArrayList<>();
        final BinomialDistribution binomialDistributionFlank = new BinomialDistribution(flankSize, pSnp);
        final BinomialDistribution binomialDistributionInner = new BinomialDistribution(length, pSnp);
        for (final Map.Entry<String, HetSnpStats> entry : sampleStats.entrySet()) {
            final String sample = entry.getKey();
            final HetSnpStats stats = entry.getValue();
            // binomial p-value on observed het counts
            final double pFlank = binomialDistributionFlank.cumulativeProbability(Math.min(stats.upstreamHets, stats.downstreamHets));
            final double pInner = binomialDistributionInner.cumulativeProbability(stats.containedHets);
            if (!(pInner < pMaxHomozygous && pFlank < pMaxHomozygous)) {
                // Not region of homozygosity
                stats.logHetCount = Math.log(stats.containedHets + 1.);
                if (carrierSamples.contains(sample)) {
                    carrierLogCounts.add(stats.logHetCount);
                } else {
                    nullLogCounts.add(stats.logHetCount);
                }
            }
        }
        if (carrierLogCounts.isEmpty() || nullLogCounts.isEmpty()) {
            return null;
        }
        // median het counts in carriers and controls
        final double medianCarrier = MEDIAN.evaluate(carrierLogCounts.stream().mapToDouble(Double::doubleValue).toArray());
        final double medianNull = MEDIAN.evaluate(nullLogCounts.stream().mapToDouble(Double::doubleValue).toArray());
        // return log ratio
        return medianNull - medianCarrier;
    }

    private static final class HetSnpStats {
        public int upstreamHets = 0;
        public int containedHets = 0;
        public int downstreamHets = 0;
        public Double logHetCount = null;
    }
}
