package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Measures quality of a CNV deletion with two metrics:
 * <ul>
 *   <li>{@code BAF_HET_RATIO}: The mean of 10^(-Ratio) across non-ROH carrier samples, where
 *       Ratio = log10((M + epsilon*L) / (min(F,E) + epsilon*L)). For a true deletion, carriers
 *       have fewer hets inside the SV (M is small) so Ratio is very negative and 10^(-Ratio) is large.</li>
 *   <li>{@code BAF_DEL_LOGLIK}: The negated average log-likelihood of carrier log-ratios under
 *       a 3-component Gaussian mixture model fit on control (non-carrier, non-ROH) log-ratios,
 *       approximating the v1.1 BayesianGaussianMixture implementation.
 *       Higher values indicate carriers are more unlikely under the null distribution.</li>
 * </ul>
 *
 * <p>This implementation matches the v1.1 svtk baf-test methodology (BAFpysam.py DeletionTest class).</p>
 */
public class BafHetRatioTester {

    /** Minimum number of non-ROH control log-ratios required to fit the GMM */
    static final int MIN_TRAINABLE_CONTROLS = 11;
    /** Minimum range of control log-ratios required for GMM fitting */
    static final double MIN_TRAINABLE_RANGE = 0.0001;
    /** Number of GMM components (matching v1.1) */
    static final int GMM_N_COMPONENTS = 3;
    /** Minimum number of het SNPs inside the SV region across all samples */
    static final int MIN_TOTAL_SNPS = 10;
    /** Default pseudocount rate for ratio computation */
    static final double DEFAULT_EPSILON = 0.0005;

    private final int seed;

    /**
     * @param seed  random seed for GMM initialization
     */
    public BafHetRatioTester(final int seed) {
        this.seed = seed;
    }

    /**
     * Computes BAF deletion metrics for a DEL record.
     *
     * @param record            query DEL record
     * @param evidence          BAF evidence including flanking regions
     * @param allSamples        all samples
     * @param carrierSamples    carrier samples
     * @param svLength          length of the SV
     * @return                  result containing BAF_HET_RATIO and BAF_DEL_LOGLIK, or null
     */
    public BafDelResult test(final SVCallRecord record, final List<BafEvidence> evidence,
                             final Set<String> allSamples, final Set<String> carrierSamples,
                             final int svLength) {
        Utils.nonNull(record);
        Utils.nonNull(allSamples);
        Utils.nonNull(carrierSamples);
        Utils.validateArg(allSamples.size() >= carrierSamples.size() && allSamples.containsAll(carrierSamples),
                "Sample set must contain all carrier samples");

        if (!record.isSimpleCNV() || evidence == null || evidence.isEmpty()) {
            return null;
        }

        // Count het SNPs per sample in upstream flank, inside SV, and downstream flank
        final Map<String, int[]> hetCounts = new LinkedHashMap<>();
        for (final String s : allSamples) {
            hetCounts.put(s, new int[3]); // [before, inside, after]
        }
        for (final BafEvidence baf : evidence) {
            final int[] counts = hetCounts.get(baf.getSample());
            if (counts == null) continue;
            if (baf.getStart() < record.getPositionA()) {
                counts[0]++;
            } else if (baf.getStart() >= record.getPositionB()) {
                counts[2]++;
            } else {
                counts[1]++;
            }
        }

        return calculate(hetCounts, carrierSamples, svLength);
    }

    /**
     * Annotates record with BAF_HET_RATIO and BAF_DEL_LOGLIK.
     */
    public SVCallRecord applyToRecord(final SVCallRecord record, final BafDelResult result) {
        Utils.nonNull(record);
        if (result == null) {
            return record;
        }
        final Map<String, Object> attributes = new HashMap<>(record.getAttributes());
        attributes.put(GATKSVVCFConstants.BAF_HET_RATIO_ATTRIBUTE, result.hetRatio);
        if (!Double.isNaN(result.delLoglik)) {
            attributes.put(GATKSVVCFConstants.BAF_DEL_LOGLIK_ATTRIBUTE, result.delLoglik);
        }
        return SVCallRecordUtils.copyCallWithNewAttributes(record, attributes);
    }

    BafDelResult calculate(final Map<String, int[]> hetCounts,
                                   final Set<String> carrierSamples,
                                   final int svLength) {
        // Effective length capped at 1M (matching v1.1)
        final int effectiveLength = Math.min(svLength, 1000000);
        final double epsilon = Math.min(50.0 / effectiveLength, DEFAULT_EPSILON);

        // Compute per-sample log-ratio and filter ROH
        final List<Double> controlRatios = new ArrayList<>();
        final List<Double> carrierRatios = new ArrayList<>();
        int totalSnps = 0;

        for (final Map.Entry<String, int[]> entry : hetCounts.entrySet()) {
            final String sample = entry.getKey();
            final int[] counts = entry.getValue();
            final int F = counts[0]; // before
            final int M = counts[1]; // inside
            final int E = counts[2]; // after
            totalSnps += M;

            // ROH detection (matching v1.1 Deltest function)
            final double threshold = epsilon;
            if ((F / (double) effectiveLength < threshold && M / (double) effectiveLength < threshold)
                    || (E / (double) effectiveLength < threshold && M / (double) effectiveLength < threshold)) {
                // ROH - skip this sample
                continue;
            }

            // Compute log-ratio
            final double flank = Math.min(F, E);
            final double ratio = Math.log10((M + epsilon * effectiveLength) / (flank + epsilon * effectiveLength));

            if (carrierSamples.contains(sample)) {
                carrierRatios.add(ratio);
            } else {
                controlRatios.add(ratio);
            }
        }

        if (carrierRatios.isEmpty()) {
            return null;
        }

        // Compute BAF_HET_RATIO = mean(10^(-Ratio)) for carriers (matching v1.1 BAF_snp_ratio)
        final double hetRatio = carrierRatios.stream()
                .mapToDouble(r -> Math.pow(10, -r))
                .average()
                .orElse(Double.NaN);

        // Check if carriers exceed controls or too few SNPs (matching v1.1 error conditions)
        if (carrierRatios.size() > controlRatios.size() || totalSnps < MIN_TOTAL_SNPS) {
            return new BafDelResult(hetRatio, Double.NaN);
        }

        // Compute BAF_DEL_LOGLIK using GMM
        final double delLoglik = computeGmmLoglik(controlRatios, carrierRatios);

        return new BafDelResult(hetRatio, delLoglik);
    }

    /**
    * Fits a 3-component GMM on control log-ratios and scores carrier log-ratios.
     * Returns the negated average log-likelihood (higher = more anomalous).
     */
    double computeGmmLoglik(final List<Double> controlRatios, final List<Double> carrierRatios) {
        if (controlRatios.size() < MIN_TRAINABLE_CONTROLS) {
            return Double.NaN;
        }
        final double[] controlArray = controlRatios.stream().mapToDouble(Double::doubleValue).toArray();
        final double min = Arrays.stream(controlArray).min().orElse(0);
        final double max = Arrays.stream(controlArray).max().orElse(0);
        if (max - min < MIN_TRAINABLE_RANGE) {
            return Double.NaN;
        }

        // Fit GMM
        final SimpleGMM1D gmm = new SimpleGMM1D(GMM_N_COMPONENTS, 20, 200, seed);
        gmm.fit(controlArray);

        // Score carrier samples: average log-likelihood under the null model
        final double[] carrierArray = carrierRatios.stream().mapToDouble(Double::doubleValue).toArray();
        final double avgLogLikelihood = gmm.score(carrierArray);

        // Negate: higher BAF_DEL_LOGLIK means carriers are more unlikely under null
        return -avgLogLikelihood;
    }

    /**
     * Result record containing both BAF metrics for deletions.
     */
    public record BafDelResult(double hetRatio, double delLoglik) {
    }
}
