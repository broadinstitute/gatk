package org.broadinstitute.hellbender.tools.sv;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.aggregation.DiscordantPairEvidenceTester;
import org.broadinstitute.hellbender.tools.sv.aggregation.EvidenceStatUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.DoubleStream;

public class DiscordantPairEvidenceGenotyper {

    private final int minSize;
    private final int minCount;
    private final double targetCoverage;
    private final int maxQuality;
    private final Map<String, Double> sampleCoverageMap;
    private final Map<String, FirstPassResult> firstPassCounts;
    private Double hetMedian = null;
    private Double hetMad = null;
    private Double hetCutoff = null;
    private Double normalization = null;
    private Double normalizationVariant = null;
    private final List<Double> hetCounts;
    private final List<Double> homCounts;
    private boolean firstPassMade = false;
    private boolean secondPassMade = false;

    private static final Median MEDIAN = new Median();
    private static final NormalDistribution Z_DISTRIBUTION = new NormalDistribution();
    private static final StandardDeviation STD_DEV = new StandardDeviation();

    public DiscordantPairEvidenceGenotyper(final Map<String,Double> sampleCoverageMap, final double qualityCutoff, final int minSize, final double targetCoverage, final int maxQuality) {
        this.sampleCoverageMap = Utils.nonNull(sampleCoverageMap);
        this.minCount = computeCountCutoff(qualityCutoff);
        Utils.validateArg(maxQuality > 0, "Maximum quality must be greater than zero");
        this.maxQuality = maxQuality;
        this.minSize = minSize;
        this.targetCoverage = targetCoverage;
        this.firstPassCounts = new HashMap<>();
        this.hetCounts = new ArrayList<>();
        this.homCounts = new ArrayList<>();
    }

    private static int computeCountCutoff(final double qualityCutoff) {
        int i = 1;
        while (true) {
            final PoissonDistribution dist = new PoissonDistribution(i);
            final double p = dist.cumulativeProbability(0);
            if (p == 0) {
                throw new IllegalArgumentException("Precision error - quality cutoff " + qualityCutoff + " is too high");
            }
            final double qual = QualityUtils.errorProbToQual(p);
            if (qual > qualityCutoff) {
                return Math.max(i - 1, 1);
            }
            i++;
        }
    }

    public void addFirstPass(final SVCallRecord record, final List<DiscordantPairEvidence> evidence) {
        Utils.nonNull(evidence);
        Utils.validate(!firstPassMade, "First pass already made");
        // TODO check all samples are in coverage map
        final Map<String, Double> counts = normalizeCounts(evidence);
        firstPassCounts.put(record.getId(), new FirstPassResult(counts));
    }

    private Map<String, Double> normalizeCounts(final List<DiscordantPairEvidence> evidence) {
        final Map<String, Integer> countsMap = DiscordantPairEvidenceTester.countEvidence(evidence);
        final Map<String, Double> counts = new HashMap<>();
        for (final String sample : countsMap.keySet()) {
            final int count = countsMap.getOrDefault(sample, 0);
            final double sampleCoverage = sampleCoverageMap.getOrDefault(sample, 0.);
            final double normalizedCount = EvidenceStatUtils.getNormalizedCount(count, sampleCoverage, targetCoverage);
            if (normalizedCount >= minCount) {
                counts.put(sample, normalizedCount);
            }
        }
        return counts;
    }

    public void finalizeFirstPass() {
        Utils.validate(!firstPassCounts.isEmpty(), "No discordant pair counts after first pass");
        Utils.validate(!firstPassMade, "First pass has already been made");
        final double[] counts = firstPassCounts.values().stream().map(FirstPassResult::getCounts).flatMapToDouble(DoubleStream::of).toArray();
        hetMedian = MEDIAN.evaluate(counts);
        final double[] deviations = DoubleStream.of(counts).map(d -> Math.abs(d - hetMedian)).toArray();
        hetMad = MEDIAN.evaluate(deviations);
        hetCutoff = hetMedian + 1.645 * hetMad;
        firstPassMade = true;
    }

    public void addSecondPass(final SVCallRecord record, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype,
                              final List<String> depthGenotypeSamples) {
        Utils.nonNull(record);
        Utils.nonNull(depthGenotype);
        Utils.nonNull(depthGenotypeSamples);
        Utils.validate(firstPassMade, "First pass must be made");
        Utils.validate(!secondPassMade, "First pass has already been made");
        if (!firstPassCounts.containsKey(record.getId())) {
            // Skip if this variant wasn't added in the first pass
            return;
        }
        final FirstPassResult firstPassResult = firstPassCounts.get(record.getId());
        final int[] copyStates = depthGenotype.copyStates();
        Utils.validateArg(copyStates.length == depthGenotypeSamples.size(), "Copy states array and sample list must be the same size");
        for (int i = 0; i < copyStates.length; i++) {
            final String sample = depthGenotypeSamples.get(i);
            final Double count = firstPassResult.getCount(sample);
            if (count != null && !(count < hetCutoff && (copyStates[i] == 0 || copyStates[i] >= 4))) {
                if (copyStates[i] == 0 || copyStates[i] == 4) {
                    homCounts.add(count);
                } else if (copyStates[i] == 1 || copyStates[i] == 3) {
                    hetCounts.add(count);
                }
            }
        }
    }

    public DiscordantPairGenotypeParameters finalizeSecondPass() {
        Utils.validate(!homCounts.isEmpty(), "No discordant pair counts after second pass");
        Utils.validate(!hetCounts.isEmpty(), "No discordant pair counts after second pass");
        Utils.validate(!secondPassMade, "Second pass has already been made");
        final double homMedian = MEDIAN.evaluate(homCounts.stream().mapToDouble(Double::valueOf).toArray());
        final double sdHet = STD_DEV.evaluate(hetCounts.stream().mapToDouble(Double::valueOf).toArray());
        normalization = maxQuality / (- 10. * Math.log10(1 - Z_DISTRIBUTION.cumulativeProbability(0.5 * homMedian / sdHet)) - 3.0103);
        final PoissonDistribution poissonDistribution = new PoissonDistribution(0.5 * homMedian);
        normalizationVariant = maxQuality / (-10. * Math.log10(poissonDistribution.cumulativeProbability(0)));
        secondPassMade = true;
        return new DiscordantPairGenotypeParameters(minCount, homMedian, sdHet);
    }

    public DiscordantPairGenotypeResult genotype(final SVCallRecord record, final List<DiscordantPairEvidence> evidence,
                                                       final DiscordantPairGenotypeParameters parameters, final List<String> samples) {
        final Map<String, Double> counts = normalizeCounts(evidence);
        final int[] genotypes = new int[samples.size()];
        final int[] genotypeQuals = new int[samples.size()];
        final List<Integer> nonRefQuals = new ArrayList<>(samples.size());
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final double discordantPairCount = counts.getOrDefault(sample, 0.);
            if (discordantPairCount < minCount) {
                genotypes[i] = 0;
            } else if (discordantPairCount <= parameters.medianHom - parameters.sdHet) {
                genotypes[i] = 1;
            } else {
                genotypes[i] = (int) ((discordantPairCount / (parameters.medianHom * 0.5)) + 0.5);
            }
            if (genotypes[i] == 0) {
                if (discordantPairCount == 0) {
                    genotypeQuals[i] = maxQuality;
                } else {
                    final PoissonDistribution dist = new PoissonDistribution(discordantPairCount);
                    genotypeQuals[i] = (int) Math.round(- 10. * Math.log10(1. - dist.cumulativeProbability(0)) * normalization);
                }
            } else {
                final double z0 = discordantPairCount - genotypes[i] * parameters.medianHom * 0.5;
                final double z1 = z0 + parameters.medianHom * 0.5;
                final double z2 = z0 - parameters.medianHom * 0.5;
                final double q0 = getGenotypeLikelihood(z0, parameters.sdHet);
                final double q1 = getGenotypeLikelihood(z1, parameters.sdHet);
                final double q2 = getGenotypeLikelihood(z2, parameters.sdHet);
                genotypeQuals[i] = (int) Math.round((Math.min(q1, q2) - q0) * normalization);
            }
            genotypeQuals[i] = Math.max(Math.min(genotypeQuals[i], maxQuality), 1);
            if (genotypes[i] != 0) {
                nonRefQuals.add(genotypeQuals[i]);
            }
        }
        final double medianCarrierQual = MEDIAN.evaluate(nonRefQuals.stream().mapToDouble(Double::valueOf).toArray());
        final PoissonDistribution variantQualDist = new PoissonDistribution(medianCarrierQual);
        final double variantQual = - 10. * Math.log10(variantQualDist.cumulativeProbability(0));
        return new DiscordantPairGenotypeResult(genotypes, genotypeQuals, variantQual);
    }

    private static double getGenotypeLikelihood(final double z, final double sdHet) {
        return -10. * Math.log10(1. - Z_DISTRIBUTION.cumulativeProbability(Math.abs(z) / sdHet));
    }

    public boolean trainableRecord(final SVCallRecord record, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype) {
        if (record.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.DEL && record.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            return false;
        }
        if (record.getLength() == null || record.getLength() < minSize) {
            return false;
        }
        if (record.isDepthOnly()) {
            return false;
        }
        if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
            boolean foundOne = false;
            boolean foundTwo = false;
            for (final int c : depthGenotype.copyStates()) {
                if (c > 2) {
                    return false;
                }
                if (c == 1) {
                    foundOne = true;
                }
                if (c == 2) {
                    foundTwo = true;
                }
            }
            if (!(foundOne && foundTwo)) {
                return false;
            }
        }
        if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            if (Arrays.stream(depthGenotype.copyStates()).anyMatch(s -> s < 2 || s > 4)) {
                return false;
            }
            boolean foundThree = false;
            boolean foundTwo = false;
            for (final int c : depthGenotype.copyStates()) {
                if (c < 2 || c > 4) {
                    return false;
                }
                if (c == 3) {
                    foundThree = true;
                }
                if (c == 2) {
                    foundTwo = true;
                }
            }
            if (!(foundThree && foundTwo)) {
                return false;
            }
        }
        return true;
    }

    private static final class FirstPassResult {
        private final Map<String, Double> counts;
        public FirstPassResult(final Map<String, Double> counts) {
            this.counts = counts;
        }

        public Double getCount(final String sample) {
            return counts.get(sample);
        }

        public double[] getCounts() {
            return counts.values().stream().mapToDouble(Double::valueOf).toArray();
        }
    }

    public record DiscordantPairGenotypeResult(int[] genotypes, int[] genotypeQuals, double variantQual) {}
    public record DiscordantPairGenotypeParameters(double minCount, double medianHom, double sdHet) {}
}
