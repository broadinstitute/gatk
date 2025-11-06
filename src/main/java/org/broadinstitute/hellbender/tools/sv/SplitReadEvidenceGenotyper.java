package org.broadinstitute.hellbender.tools.sv;

import autovalue.shaded.com.google.common.collect.Sets;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.aggregation.EvidenceStatUtils;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.DoubleStream;

public class SplitReadEvidenceGenotyper {

    private final int minSize;
    private final int countCutoff;
    private final double targetCoverage;
    private final int maxQuality;
    private final Map<String, Double> sampleCoverageMap;
    private final Map<String, FirstPassResult> firstPassCounts;
    private Double hetMedian = null;
    private Double hetMad = null;
    private Double hetCutoff = null;
    private final List<Double> hetCounts;
    private final List<Double> homCounts;
    private boolean firstPassMade = false;
    private boolean secondPassMade = false;

    private static final Median MEDIAN = new Median();
    private static final NormalDistribution Z_DISTRIBUTION = new NormalDistribution();
    private static final StandardDeviation STD_DEV = new StandardDeviation();

    private static final Set<Integer> HET_COPY_STATES = Set.of(1, 3);
    private static final Set<Integer> HOM_COPY_STATES = Set.of(0, 4);

    public SplitReadEvidenceGenotyper(final Map<String,Double> sampleCoverageMap, final double qualityCutoff, final int minSize, final double targetCoverage, final int maxQuality) {
        this.sampleCoverageMap = Utils.nonNull(sampleCoverageMap);
        this.countCutoff = computeCountCutoff(qualityCutoff);
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

    public void addFirstPass(final SVCallRecord record, final List<SplitReadEvidence> startEvidence,
                             final List<SplitReadEvidence> endEvidence, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype, final List<String> depthGenotypeSamples) {
        Utils.nonNull(startEvidence);
        Utils.nonNull(endEvidence);
        Utils.validate(!firstPassMade, "First pass already made");
        // TODO check all samples are in coverage map
        final int minCount = Math.max(countCutoff / 2, 1);
        final Map<String, Double> startCounts = normalizeCounts(startEvidence);
        final Map<String, Double> endCounts = normalizeCounts(endEvidence);
        if (hasBothSideSupport(startCounts, endCounts, minCount)) {
            firstPassCounts.put(record.getId(), new FirstPassResult(startCounts, endCounts, depthGenotype, depthGenotypeSamples));
        }
    }

    // Assumes samples with 0 counts are not present in the key set
    private static boolean hasBothSideSupport(final Map<String, Double> startCounts,
                                              final Map<String, Double> endCounts, final double threshold) {
        for (final Map.Entry<String, Double> entry : startCounts.entrySet()) {
            if (entry.getValue() >= threshold) {
                if (endCounts.getOrDefault(entry.getKey(), 0.0) >= threshold) {
                    return true;
                }
            }
        }
        return false;
    }

    // Assumes samples with 0 counts are not present in the key set
    private static int countBothSideSupport(final Map<String, Double> startCounts,
                                            final Map<String, Double> endCounts, final double threshold) {
        int count = 0;
        for (final Map.Entry<String, Double> entry : startCounts.entrySet()) {
            if (entry.getValue() > threshold) {
                if (endCounts.getOrDefault(entry.getKey(), 0.0) > threshold) {
                    count++;
                }
            }
        }
        return count;
    }

    // Assumes samples with 0 counts are not present in the key set
    private static int countSummedSupport(final Map<String, Double> startCounts,
                                          final Map<String, Double> endCounts, final double threshold) {
        int count = 0;
        final Set<String> samples = Sets.union(startCounts.keySet(), endCounts.keySet());
        for (final String sample : samples) {
            final double sum = startCounts.getOrDefault(sample, 0.0) + endCounts.getOrDefault(sample, 0.0);
            if (sum > threshold) {
                count++;
            }
        }
        return count;
    }

    private Map<String, Double> normalizeCounts(final List<SplitReadEvidence> evidence) {
        final Map<String, Double> counts = new HashMap<>();
        for (final SplitReadEvidence e : evidence) {
            final double sampleCoverage = sampleCoverageMap.getOrDefault(e.getSample(), 0.);
            final double normalizedCount = EvidenceStatUtils.getNormalizedCount(e.getCount(), sampleCoverage, targetCoverage);
            if (normalizedCount > 0) {
                counts.put(e.getSample(), normalizedCount);
            }
        }
        return counts;
    }

    public void finalizeFirstPass() {
        Utils.validate(!firstPassCounts.isEmpty(), "No split read counts after first pass");
        Utils.validate(!firstPassMade, "First pass has already been made");
        final double[] hetCounts = firstPassCounts.values().stream().map(c -> c.getCounts(HET_COPY_STATES)).flatMapToDouble(DoubleStream::of).toArray();
        hetMedian = MEDIAN.evaluate(hetCounts);
        final double[] deviations = DoubleStream.of(hetCounts).map(d -> Math.abs(d - hetMedian)).toArray();
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
        if (firstPassResult == null) {
            throw new IllegalArgumentException("Record " + record.getId() + " wasn't added in the first pass");
        }
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

    public SplitReadGenotypeParameters finalizeSecondPass() {
        Utils.validate(!homCounts.isEmpty(), "No split read counts after second pass");
        Utils.validate(!hetCounts.isEmpty(), "No split read counts after second pass");
        Utils.validate(!secondPassMade, "Second pass has already been made");
        final double homMedian = MEDIAN.evaluate(homCounts.stream().mapToDouble(Double::valueOf).toArray());
        final double hetMedian = MEDIAN.evaluate(hetCounts.stream().mapToDouble(Double::valueOf).toArray());
        final double sdHet = 1.645 * MEDIAN.evaluate(hetCounts.stream().mapToDouble(d -> Math.abs(d - hetMedian)).toArray());
        secondPassMade = true;
        return new SplitReadGenotypeParameters(countCutoff, homMedian, sdHet);
    }

    public SplitReadGenotypeResult genotype(final SVCallRecord record,
                                            final List<SplitReadEvidence> startEvidence,
                                            final List<SplitReadEvidence> endEvidence,
                                            final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype,
                                            final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult discordantPairGenotype,
                                            final SplitReadGenotypeParameters parameters,
                                            final List<String> samples) {
        final Map<String, Double> startCounts = normalizeCounts(startEvidence);
        final Map<String, Double> endCounts = normalizeCounts(endEvidence);
        final int[] genotypes = new int[samples.size()];
        final double[] countSum = new double[samples.size()];
        final int[] genotypeQuals = new int[samples.size()];
        final List<Integer> nonRefQuals = new ArrayList<>(samples.size());
        int nonRefCount = 0;
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final double startCount = startCounts.getOrDefault(sample, 0.);
            final double endCount = endCounts.getOrDefault(sample, 0.);
            countSum[i] = startCount + endCount;if (countSum[i] < parameters.minCount) {
                genotypes[i] = 0;
            } else if (countSum[i] <= parameters.medianHom - parameters.sdHet) {
                genotypes[i] = 1;
            } else {
                genotypes[i] = (int) ((countSum[i] / (parameters.medianHom * 0.5)) + 0.5);
            }
            if (genotypes[i] != 0) {
                ++nonRefCount;
            }
        }
        final boolean largeEnough = record.getLength() != null && record.getLength() >= minSize;
        final boolean hasSplitReadEvidence = record.getEvidence().contains(GATKSVVCFConstants.EvidenceTypes.SR);
        final int minCount = Math.max(countCutoff / 2, 1);
        final int twoSidedPassCount = countBothSideSupport(startCounts, endCounts, minCount);
        final int bothsideNonZeroCount = countBothSideSupport(startCounts, endCounts, 0);
        final int samplesOverOneCount = countSummedSupport(startCounts, endCounts, 1);
        final int rareMin = 0;
        final int rareMax = Math.max(samples.size() / 100, 2);
        final int commonMin = rareMax;
        final int commonMax = samples.size();

        final List<Double> recoverSinglePass = new ArrayList<>();
        final List<Double> recoverSingleFail = new ArrayList<>();
        final List<Double> recoverBothPass = new ArrayList<>();
        final List<Double> recoverBothFail = new ArrayList<>();
        final double backgroundRatio = twoSidedPassCount / (double) bothsideNonZeroCount;
        final double genotypeRatio = nonRefCount / (double) samplesOverOneCount;
        for (int i = 0; i < samples.size(); i++) {
            // TODO this is a bug I think - should check for non-ref genotype not GQ>0
            final boolean nonRefDiscordantPair = discordantPairGenotype != null && discordantPairGenotype.genotypeQuals()[i] > 0;
            final boolean nonRefDepth = depthGenotype != null && depthGenotype.copyStates()[i] != 2;
            final boolean pass = nonRefCount > 0 && largeEnough && hasSplitReadEvidence && (nonRefDiscordantPair || nonRefDepth); // TODO should || be && ?

            if (record.getId().equals("all_samples_manta_chr20_00000008") && samples.get(i).equals("HG03085")) {
                int x = 0;
            }
            final double startCount = startCounts.getOrDefault(samples.get(i), 0.0);
            final double endCount = endCounts.getOrDefault(samples.get(i), 0.0);
            if (bothsideNonZeroCount > 0 && twoSidedPassCount > 0) {
                // compute ratio with number of fully supported samples
                // recover.bothsides.txt
                if (genotypes[i] > 0) {
                    if (pass) {
                        // recover.both.txt
                        recoverBothPass.add(backgroundRatio);
                    } else {
                        // recover.both.fail.txt
                        recoverBothFail.add(backgroundRatio);
                    }
                }
            }
            if (samplesOverOneCount > 0 && nonRefCount > 0) {
                // One-sided support
                if (genotypes[i] > 0) {
                    // recover.txt
                    if (pass) {
                        // recover.single.txt
                        recoverSinglePass.add(genotypeRatio);
                    } else {
                        // recover.single.fail.txt
                        recoverSingleFail.add(genotypeRatio);
                    }
                }
            }
        }

        if (record.getId().equals("all_samples_wham_chr20_0000000e")) {
            int x = 0;
        }

        // Quality scoring
        final double normalization = 1; // TODO
        final double normalizationVariant = 1; // TODO
        for (int i = 0; i < samples.size(); i++) {
            if (genotypes[i] == 0) {
                if (countSum[i] == 0) {
                    genotypeQuals[i] = maxQuality;
                } else {
                    final PoissonDistribution dist = new PoissonDistribution(countSum[i]);
                    genotypeQuals[i] = (int) Math.round(-10. * Math.log10(1. - dist.cumulativeProbability(0)) * normalization);
                }
            } else {
                final double z0 = countSum[i] - genotypes[i] * parameters.medianHom * 0.5;
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
        final double variantQual = -10. * Math.log10(variantQualDist.cumulativeProbability(0)) * normalizationVariant;
        return new SplitReadGenotypeResult(genotypes, genotypeQuals, variantQual);
    }

    private static double getGenotypeLikelihood(final double z, final double sdHet) {
        return -10. * Math.log10(1. - Z_DISTRIBUTION.cumulativeProbability(Math.abs(z) / sdHet));
    }

    public boolean trainableRecord(final SVCallRecord record,
                                   final DiscordantPairEvidenceGenotyper discordantPairGenotyper,
                                   final SVStratificationEngine exclusionEngine) {
        if (!discordantPairGenotyper.isTrainingRecord(record)) {
            return false;
        }
        if (record.getLength() == null || record.getLength() < minSize) {
            return false;
        }
        if (exclusionEngine != null && !exclusionEngine.getMatches(record, 0, 1, 1).isEmpty()) {
            return false;
        }
        if (!record.getEvidence().contains(GATKSVVCFConstants.EvidenceTypes.SR)) {
            return false;
        }
        return true;
    }

    private static final class FirstPassResult {
        private final Map<String, Double> counts;
        private final Map<String, Integer> depthGenotypes;
        public FirstPassResult(final Map<String, Double> startCounts, final Map<String, Double> endCounts, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype, final List<String> depthGenotypeSamples) {
            this.counts = new HashMap<>();
            final Set<String> keys = Sets.union(startCounts.keySet(), endCounts.keySet());
            for (final String key : keys) {
                final double startCount = startCounts.getOrDefault(key, 0.);
                final double endCount = endCounts.getOrDefault(key, 0.);
                counts.put(key, startCount + endCount);
            }
            this.depthGenotypes = new HashMap<>();
            for (int i = 0; i < depthGenotypeSamples.size(); i++) {
                final String sample = depthGenotypeSamples.get(i);
                if (counts.containsKey(sample)) {
                    depthGenotypes.put(sample, depthGenotype.copyStates()[i]);
                }
            }
        }

        public Double getCount(final String sample) {
            return counts.get(sample);
        }

        public double[] getCounts(final Set<Integer> validDepthStates) {
            final List<Double> countsList = new ArrayList<>(counts.size());
            for (final Map.Entry<String, Double> entry : counts.entrySet()) {
                if (depthGenotypes.containsKey(entry.getKey())) {
                    final int copyState = depthGenotypes.get(entry.getKey());
                    if (validDepthStates.contains(copyState)) {
                        countsList.add(entry.getValue());
                    }
                } else {
                    throw new IllegalArgumentException("Sample '" + entry.getKey() + "' does not exist in depth genotype");
                }
            }
            return countsList.stream().mapToDouble(Double::valueOf).toArray();
        }
    }

    public record SplitReadGenotypeResult(int[] genotypes, int[] genotypeQuals, double variantQual) {}
    public record SplitReadGenotypeParameters(double minCount, double medianHom, double sdHet) {}
}
