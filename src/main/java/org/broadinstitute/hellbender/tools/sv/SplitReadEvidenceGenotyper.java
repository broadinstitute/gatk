package org.broadinstitute.hellbender.tools.sv;

import autovalue.shaded.com.google.common.collect.Sets;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.aggregation.EvidenceStatUtils;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.*;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public class SplitReadEvidenceGenotyper {

    private final Integer minSize;
    private final int trainingCountCutoff;
    private final double targetCoverage;
    private final int maxQuality;
    private final Map<String, Double> sampleCoverageMap;
    private final Map<String, FirstPassResult> firstPassCounts = new HashMap<>();
    private Double hetMedian = null;
    private Double hetMad = null;
    private Double hetCutoff = null;
    private final List<Double> hetCounts = new ArrayList<>();
    private final List<Double> homCounts = new ArrayList<>();
    private boolean firstPassMade = false;
    private boolean secondPassMade = false;
    private boolean thirdPassMade = false;

    private final int rareMin;
    private final int rareMax;
    private final int commonMin;
    private final int commonMax;

    private final List<RecoverResult> recoverSinglePass = new ArrayList<>();
    private final List<RecoverResult> recoverSingleFail = new ArrayList<>();
    private final List<RecoverResult> recoverBothPass = new ArrayList<>();
    private final List<RecoverResult> recoverBothFail = new ArrayList<>();

    private static final Median MEDIAN = new Median();
    private static final NormalDistribution Z_DISTRIBUTION = new NormalDistribution();

    private static final Set<Integer> HET_COPY_STATES = Set.of(1, 3);

    public SplitReadEvidenceGenotyper(final Map<String,Double> sampleCoverageMap, final int numSamples, final double qualityCutoff, final Integer minSize, final double targetCoverage, final int maxQuality) {
        this.sampleCoverageMap = Utils.nonNull(sampleCoverageMap);
        this.trainingCountCutoff = computeCountCutoff(qualityCutoff);
        Utils.validateArg(maxQuality > 0, "Maximum quality must be greater than zero");
        this.maxQuality = maxQuality;
        this.minSize = minSize;
        this.targetCoverage = targetCoverage;
        this.rareMin = 0;
        this.rareMax = Math.max(numSamples / 100, 2);
        this.commonMin = rareMax;
        this.commonMax = numSamples;
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
        final int minCount = Math.max(trainingCountCutoff / 2, 1);
        final Map<String, Double> startCounts = normalizeCounts(startEvidence);
        final Map<String, Double> endCounts = normalizeCounts(endEvidence);
        final Set<String> bothSideSupportSamples = hasBothSideSupport(startCounts, endCounts, minCount);
        if (!bothSideSupportSamples.isEmpty()) {
            firstPassCounts.put(record.getId(), new FirstPassResult(bothSideSupportSamples, startCounts, endCounts, depthGenotype, depthGenotypeSamples));
        }
    }

    // Assumes samples with 0 counts are not present in the key set
    private static Set<String> hasBothSideSupport(final Map<String, Double> startCounts,
                                              final Map<String, Double> endCounts, final double threshold) {
        final Set<String> result = new HashSet<>();
        for (final Map.Entry<String, Double> entry : startCounts.entrySet()) {
            if (entry.getValue() >= threshold) {
                if (endCounts.getOrDefault(entry.getKey(), 0.0) >= threshold) {
                    result.add(entry.getKey());
                }
            }
        }
        return result;
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

    private record RecoverResult(String variant, String sample, int count, double frac) {}
    public record CutoffResult(double fracSingle, double fracBoth, int countPass, int countFail, int freqMin, int freqMax) {}

    /**
     * For a sorted array of frac values (ascending), count the number of elements >= threshold
     * using binary search. Returns the count of elements in the suffix starting from the
     * leftmost element >= threshold.
     */
    private static int countAboveThreshold(final double[] sortedFracs, final double threshold) {
        if (sortedFracs.length == 0) {
            return 0;
        }
        // Binary search for leftmost index where sortedFracs[idx] >= threshold
        int lo = 0;
        int hi = sortedFracs.length;
        while (lo < hi) {
            final int mid = (lo + hi) >>> 1;
            if (sortedFracs[mid] < threshold) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        return sortedFracs.length - lo;
    }

    /**
     * Pre-filter a recovery list by count range, sort the resulting frac values,
     * then compute the count of elements with frac >= cutoff for each cutoff value
     * using binary search. This replaces N linear scans with one sort + N binary searches.
     */
    private static int[] precomputeCountsForCutoffs(final List<RecoverResult> list, final double[] cutoffs,
                                                     final int freqMin, final int freqMax) {
        final double[] filteredFracs = list.stream()
                .filter(r -> r.count > freqMin && r.count <= freqMax)
                .mapToDouble(RecoverResult::frac)
                .sorted()
                .toArray();
        final int[] counts = new int[cutoffs.length];
        for (int i = 0; i < cutoffs.length; i++) {
            counts[i] = countAboveThreshold(filteredFracs, cutoffs[i]);
        }
        return counts;
    }

    public SplitReadGenotypeFrequencyCutoffs finalizeThirdPass() {
        Utils.validate(!thirdPassMade, "Third pass has already been made");
        final double[] cutoffs = IntStream.rangeClosed(0, 10).mapToDouble(i -> i * 0.1).toArray();
        final CutoffResult commonResult = cutoffOptimization(cutoffs, commonMin, commonMax);
        final CutoffResult rareResult = cutoffOptimization(cutoffs, rareMin, rareMax);
        thirdPassMade = true;
        // Free recovery lists that are no longer needed
        recoverSinglePass.clear();
        recoverSingleFail.clear();
        recoverBothPass.clear();
        recoverBothFail.clear();
        return new SplitReadGenotypeFrequencyCutoffs(rareResult, commonResult);
    }

    /**
     * Optimized cutoff search: pre-sort each recovery list's frac values (filtered by count range)
     * once, then use binary search for each cutoff threshold. This reduces the complexity from
     * O(121 * M) to O(M * log(M) + 44 * log(M)) where M is the recovery list size.
     */
    private CutoffResult cutoffOptimization(final double[] cutoffs, final int freqMin, final int freqMax) {
        // Precompute counts for each cutoff value via sort + binary search
        final int[] singlePassCounts = precomputeCountsForCutoffs(recoverSinglePass, cutoffs, freqMin, freqMax);
        final int[] singleFailCounts = precomputeCountsForCutoffs(recoverSingleFail, cutoffs, freqMin, freqMax);
        final int[] bothPassCounts = precomputeCountsForCutoffs(recoverBothPass, cutoffs, freqMin, freqMax);
        final int[] bothFailCounts = precomputeCountsForCutoffs(recoverBothFail, cutoffs, freqMin, freqMax);

        final List<CutoffResult> combine = new ArrayList<>(cutoffs.length * cutoffs.length);
        for (int s = 0; s < cutoffs.length; s++) {
            for (int b = 0; b < cutoffs.length; b++) {
                final int passCount = singlePassCounts[s] + bothPassCounts[b];
                final int failCount = singleFailCounts[s] + bothFailCounts[b];
                combine.add(new CutoffResult(cutoffs[s], cutoffs[b], passCount, failCount, freqMin, freqMax));
            }
        }
        final double baseline = computeBaseline(combine);
        final int maxIndex = MathUtils.maxElementIndex(combine.stream().mapToDouble(c -> computeCutoffScore(c, baseline)).toArray());
        return combine.get(maxIndex);
    }

    private static double computeBaseline(final List<CutoffResult> list) {
        for (final CutoffResult result : list) {
            if (result.fracSingle == 0 && result.fracBoth == 0) {
                return result.countPass;
            }
        }
        throw new IllegalArgumentException("List did not contain 0-fraction entry");
    }

    final double computeCutoffScore(final CutoffResult cutoffResult, final double baseline) {
        final double a = cutoffResult.countFail / (double) (cutoffResult.countFail + cutoffResult.countPass);
        final double b = (cutoffResult.countPass / baseline) - 1;
        return -((a * a) + (b * b));
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
        final double[] homArr = homCounts.stream().mapToDouble(Double::valueOf).toArray();
        final double[] hetArr = hetCounts.stream().mapToDouble(Double::valueOf).toArray();
        final double homMedian = MEDIAN.evaluate(homArr);
        final double hetMedian = MEDIAN.evaluate(hetArr);
        final double hetMadValue = MEDIAN.evaluate(DoubleStream.of(hetArr).map(d -> Math.abs(d - hetMedian)).toArray());
        final double sdHet = 1.645 * 1.4826 * hetMadValue;
        secondPassMade = true;
        // Free training accumulation data that is no longer needed
        firstPassCounts.clear();
        hetCounts.clear();
        homCounts.clear();
        return new SplitReadGenotypeParameters(trainingCountCutoff, homMedian, sdHet);
    }
    public SplitReadGenotypeResult genotype(final SVCallRecord record,
                                                    final List<SplitReadEvidence> startEvidence,
                                                    final List<SplitReadEvidence> endEvidence,
                                                    final SplitReadGenotypeMetrics metrics,
                                                    final int medianHomIns,
                                                    final double medianHomCutoffMultiplier,
                                                    final List<String> samples) {
        final SplitReadGenotypeParameters parameters = metrics.parameters;
        final SplitReadGenotypeFrequencyCutoffs frequencyCutoffs = metrics.cutoffs;
        final Map<String, Double> startCounts = normalizeCounts(startEvidence);
        final Map<String, Double> endCounts = normalizeCounts(endEvidence);
        final int[] genotypes = new int[samples.size()];
        final double[] countSum = new double[samples.size()];
        final int[] genotypeQuals = new int[samples.size()];
        final GATKSVVCFConstants.StructuralVariantAnnotationType svtype = record.getType();
        final double medianHom = svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.INS ? medianHomIns : parameters.medianHom;
        final double bothsideMedianHet = 0.5 * medianHom;
        final double bothsideHomCutoff = medianHomCutoffMultiplier * bothsideMedianHet;
        final double bothsideSdHet = parameters.sdHet;
        final double onesideMedianHet = 0.5 * bothsideMedianHet;
        final double onesideHomCutoff = medianHomCutoffMultiplier * onesideMedianHet;
        final double onesideSdHet = 0.5 * parameters.sdHet;
        final double twoSidedThreshold = parameters.minCount / 2.;
        int nonRefCount = 0;
        final boolean[] twoSided = new boolean[samples.size()];
        final List<Double> bothsideNonRefCounts = new ArrayList<>();
        final List<Double> onesideNonRefCounts = new ArrayList<>();
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final double startCount = startCounts.getOrDefault(sample, 0.);
            final double endCount = endCounts.getOrDefault(sample, 0.);
            countSum[i] = startCount + endCount;
            twoSided[i] = startCount > twoSidedThreshold && endCount > twoSidedThreshold;
            final double medianHet;
            final double homCutoff;
            if (twoSided[i]) {
                medianHet = bothsideMedianHet;
                homCutoff = bothsideHomCutoff;
            } else {
                medianHet = onesideMedianHet;
                homCutoff = onesideHomCutoff;
            }
            if (countSum[i] < parameters.minCount) {
                genotypes[i] = 0;
            } else if (countSum[i] <= homCutoff) {
                genotypes[i] = 1;
            } else {
                genotypes[i] = Math.max((int) ((countSum[i] / medianHet) + 0.5), 2);
            }
            if (genotypes[i] != 0) {
                ++nonRefCount;
                if (twoSided[i]) {
                    bothsideNonRefCounts.add(countSum[i]);
                } else {
                    onesideNonRefCounts.add(countSum[i]);
                }
            }
        }
        final boolean hasSplitReadEvidence = record.getEvidence().contains(GATKSVVCFConstants.EvidenceTypes.SR);
        final int twosidedMinCount = Math.max((int) parameters.minCount()/ 2, 1);
        final int twoSidedPassCount = countBothSideSupport(startCounts, endCounts, twosidedMinCount);
        final int bothsideNonZeroCount = countBothSideSupport(startCounts, endCounts, 0);
        final int samplesOverOneCount = countSummedSupport(startCounts, endCounts, 1);
        boolean bothsidePass = false;
        boolean onesidePass = false;
        if (bothsideNonZeroCount > 0) {
            final double backgroundRatio = twoSidedPassCount / (double) bothsideNonZeroCount;
            bothsidePass = (backgroundRatio >= frequencyCutoffs.rare.fracBoth && bothsideNonZeroCount <= frequencyCutoffs.rare.freqMax) ||
                    (backgroundRatio >= frequencyCutoffs.common.fracBoth && bothsideNonZeroCount >= frequencyCutoffs.common.freqMin);
        }
        if (samplesOverOneCount > 0) {
            final double genotypeRatio = nonRefCount / (double) samplesOverOneCount;
            onesidePass = (genotypeRatio >= frequencyCutoffs.rare.fracSingle && nonRefCount <= frequencyCutoffs.rare.freqMax) ||
                    (genotypeRatio >= frequencyCutoffs.common.fracSingle && nonRefCount >= frequencyCutoffs.common.freqMin);
        }
        final boolean backgroundFail = !(onesidePass || bothsidePass) && hasSplitReadEvidence && nonRefCount > 0;
        final double normalization = maxQuality / Z_DISTRIBUTION.cumulativeProbability(0);
        for (int i = 0; i < samples.size(); i++) {
            final double medianHet;
            final double sdHet;
            if (twoSided[i]) {
                medianHet = bothsideMedianHet;
                sdHet = bothsideSdHet;
            } else {
                medianHet = onesideMedianHet;
                sdHet = onesideSdHet;
            }
            if (genotypes[i] == 0) {
                if (countSum[i] == 0) {
                    genotypeQuals[i] = maxQuality;
                } else {
                    final PoissonDistribution dist = new PoissonDistribution(countSum[i]);
                    genotypeQuals[i] = Math.min((int) Math.round((normalization * (1. - dist.cumulativeProbability(0)))), maxQuality);

                }
            } else if (genotypes[i] == 1) {
                genotypeQuals[i] = Math.min((int) Math.round(normalization * Z_DISTRIBUTION.cumulativeProbability((countSum[i] - medianHet) / sdHet)), maxQuality);
            } else {
                genotypeQuals[i] = Math.min((int) Math.round(normalization * Z_DISTRIBUTION.cumulativeProbability((countSum[i] - 2. * medianHet) / sdHet)), maxQuality);
            }
        }
        final int onesideVariantQuality = computeVariantQuality(onesideNonRefCounts, onesideMedianHet);
        final int bothsideVariantQuality = computeVariantQuality(bothsideNonRefCounts, bothsideMedianHet);
        int variantQuality = 0;
        if (bothsideNonZeroCount > 0) {
            final double backgroundRatio = twoSidedPassCount / (double) bothsideNonZeroCount;
            variantQuality = (int) Math.round(backgroundRatio * bothsideVariantQuality + (1. - backgroundRatio) * onesideVariantQuality);
        } else {
            variantQuality = onesideVariantQuality;
        }
        return new SplitReadGenotypeResult(genotypes, genotypeQuals, variantQuality, bothsidePass, backgroundFail);
    }

    private int computeVariantQuality(final List<Double> nonRefCounts, final double medianHet) {
        final PoissonDistribution poissonDistribution = new PoissonDistribution(medianHet);
        final double normalizationVariant = maxQuality / (-10. * Math.log10(poissonDistribution.cumulativeProbability(0)));
        if (!nonRefCounts.isEmpty()) {
            final double median = MEDIAN.evaluate(nonRefCounts.stream().mapToDouble(Double::doubleValue).toArray());
            return (int) Math.round(Math.min(-10. * normalizationVariant * Math.log10(new PoissonDistribution(median).cumulativeProbability(0)), maxQuality));
        } else {
            return 0;
        }
    }

    public SplitReadGenotypeResult genotypeTraining(final SVCallRecord record,
                                                    final List<SplitReadEvidence> startEvidence,
                                                    final List<SplitReadEvidence> endEvidence,
                                                    final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype,
                                                    final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult discordantPairGenotype,
                                                    final SplitReadGenotypeParameters parameters,
                                                    final List<String> samples) {
        // TODO parameters.minCount and this.countCutoff are redundant ?
        final Map<String, Double> startCounts = normalizeCounts(startEvidence);
        final Map<String, Double> endCounts = normalizeCounts(endEvidence);
        final int[] genotypes = new int[samples.size()];
        final double[] countSum = new double[samples.size()];
        int nonRefCount = 0;
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final double startCount = startCounts.getOrDefault(sample, 0.);
            final double endCount = endCounts.getOrDefault(sample, 0.);
            countSum[i] = startCount + endCount;
            if (countSum[i] < parameters.minCount) {
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
        final int minCount = Math.max(trainingCountCutoff / 2, 1);
        final int twoSidedPassCount = countBothSideSupport(startCounts, endCounts, minCount);
        final int bothsideNonZeroCount = countBothSideSupport(startCounts, endCounts, 0);
        final int samplesOverOneCount = countSummedSupport(startCounts, endCounts, 1);
        for (int i = 0; i < samples.size(); i++) {
            // Note: nonRefDiscordantPair checks GQ > 0, which is always true since GQ is clamped to min 1.
            // This matches v1.1 behavior (SR_genotype.opt_part1.sh checks $NF>0 on pe.geno.withquality.txt.gz,
            // where PE GQ is also clamped to min 1), making pass/fail effectively per-variant not per-sample.
            final boolean nonRefDiscordantPair = discordantPairGenotype != null && discordantPairGenotype.genotypeQuals()[i] > 0;
            final boolean nonRefDepth = depthGenotype != null && depthGenotype.copyStates()[i] != 2;
            final boolean pass = depthGenotype != null && nonRefCount > 0 && largeEnough && hasSplitReadEvidence && (nonRefDiscordantPair || nonRefDepth);
            if (bothsideNonZeroCount > 0 && twoSidedPassCount > 0) {
                final double backgroundRatio = twoSidedPassCount / (double) bothsideNonZeroCount;
                // compute ratio with number of fully supported samples
                // recover.bothsides.txt
                if (genotypes[i] > 0) {
                    final RecoverResult result = new RecoverResult(record.getId(), samples.get(i), twoSidedPassCount, backgroundRatio);
                    if (pass) {
                        // recover.both.txt
                        recoverBothPass.add(result);
                    } else {
                        // recover.both.fail.txt
                        recoverBothFail.add(result);
                    }
                }
            } else if (samplesOverOneCount > 0 && nonRefCount > 0) {  // in original code, this is not an "else" but records are deduplicated later in optimalsrcutoffs.sh
                final double genotypeRatio = nonRefCount / (double) samplesOverOneCount;
                // One-sided support
                if (genotypes[i] > 0) {
                    final RecoverResult result = new RecoverResult(record.getId(), samples.get(i), nonRefCount, genotypeRatio);
                    // recover.txt
                    if (pass) {
                        // recover.single.txt
                        recoverSinglePass.add(result);
                    } else {
                        // recover.single.fail.txt
                        recoverSingleFail.add(result);
                    }
                }
            }
        }
        return new SplitReadGenotypeResult(genotypes, null, maxQuality, null, null);
    }

    public boolean trainableRecord(final SVCallRecord record,
                                   final boolean discordantPairEligible,
                                   final SVStratificationEngine exclusionEngine) {
        if (!discordantPairEligible) {
            return false;
        }
        if (minSize != null && (record.getLength() == null || record.getLength() < minSize)) {
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

    public boolean trainableRecord(final SVCallRecord record,
                                   final DiscordantPairEvidenceGenotyper discordantPairGenotyper,
                                   final SVStratificationEngine exclusionEngine) {
        return trainableRecord(record, discordantPairGenotyper.isTrainingRecord(record), exclusionEngine);
    }

    private static final class FirstPassResult {
        private final Map<String, Double> counts;
        private final Map<String, Integer> depthGenotypes;
        public FirstPassResult(final Set<String> passingSamples, final Map<String, Double> startCounts, final Map<String, Double> endCounts, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype, final List<String> depthGenotypeSamples) {
            this.counts = new HashMap<>();
            final Set<String> keys = Sets.union(startCounts.keySet(), endCounts.keySet());
            for (final String key : keys) {
                if (passingSamples.contains(key)) {
                    final double startCount = startCounts.getOrDefault(key, 0.);
                    final double endCount = endCounts.getOrDefault(key, 0.);
                    counts.put(key, startCount + endCount);
                }
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

    public record SplitReadGenotypeResult(int[] genotypes, int[] genotypeQuals, double variantQual, Boolean bothsidePass, Boolean backgroundFail) {}
    public record SplitReadGenotypeParameters(double minCount, double medianHom, double sdHet) {}
    public record SplitReadGenotypeFrequencyCutoffs(CutoffResult rare, CutoffResult common) {}

    public static final class SplitReadTableParser {

        private static final String MIN_COUNT_COLUMN = "sr_count";
        private static final String MEDIAN_HOM_COLUMN = "median_hom";
        private static final String SD_HET_COLUMN = "sd_het";
        private static final String RARE_MIN_COLUMN = "rare_min";
        private static final String RARE_MAX_COLUMN = "rare_max";
        private static final String COMMON_MIN_COLUMN = "common_min";
        private static final String COMMON_MAX_COLUMN = "common_max";
        private static final String RARE_PASS_COLUMN = "rare_pass";
        private static final String RARE_FAIL_COLUMN = "rare_fail";
        private static final String COMMON_PASS_COLUMN = "common_pass";
        private static final String COMMON_FAIL_COLUMN = "common_fail";
        private static final String RARE_SINGLE_COLUMN = "rare_single";
        private static final String RARE_BOTH_COLUMN = "rare_both";
        private static final String COMMON_SINGLE_COLUMN = "common_single";
        private static final String COMMON_BOTH_COLUMN = "common_both";
        public static final TableColumnCollection CUTOFFS_COLUMNS = new TableColumnCollection(Arrays.asList(
                MIN_COUNT_COLUMN, MEDIAN_HOM_COLUMN, SD_HET_COLUMN, RARE_MIN_COLUMN, RARE_MAX_COLUMN, COMMON_MIN_COLUMN,
                COMMON_MAX_COLUMN, RARE_PASS_COLUMN, RARE_FAIL_COLUMN, COMMON_PASS_COLUMN, COMMON_FAIL_COLUMN,
                RARE_SINGLE_COLUMN, RARE_BOTH_COLUMN, COMMON_SINGLE_COLUMN, COMMON_BOTH_COLUMN));

        public void composeCutoffsLine(final SplitReadGenotypeMetrics stats,
                                        final DataLine dataLine) {
            dataLine.append(stats.parameters().minCount());
            dataLine.append(stats.parameters().medianHom());
            dataLine.append(stats.parameters().sdHet());
            dataLine.append(stats.cutoffs().rare().freqMin());
            dataLine.append(stats.cutoffs().rare().freqMax());
            dataLine.append(stats.cutoffs().common().freqMin());
            dataLine.append(stats.cutoffs().common().freqMax());
            dataLine.append(stats.cutoffs().rare().countPass());
            dataLine.append(stats.cutoffs().rare().countFail());
            dataLine.append(stats.cutoffs().common().countPass());
            dataLine.append(stats.cutoffs().common().countFail());
            dataLine.append(stats.cutoffs().rare().fracSingle());
            dataLine.append(stats.cutoffs().rare().fracBoth());
            dataLine.append(stats.cutoffs().common().fracSingle());
            dataLine.append(stats.cutoffs().common().fracBoth());
        }

        public Function<DataLine, SplitReadGenotypeMetrics> tableParser(TableColumnCollection columns, Function<String, RuntimeException> exceptionFactory) {
            // Check for expected columns
            for (final String column : CUTOFFS_COLUMNS.names()) {
                if (!columns.contains(column)) {
                    throw exceptionFactory.apply("Missing column " + column);
                }
            }
            // Check there are no extra columns
            if (columns.columnCount() != CUTOFFS_COLUMNS.columnCount()) {
                throw exceptionFactory.apply("Expected " + columns.columnCount() + " columns but found " + columns.columnCount());
            }
            return this::parseTableLine;
        }

        public SplitReadGenotypeMetrics parseTableLine(final DataLine dataLine) {
            final double minCount = Double.parseDouble(dataLine.get(0));
            final double medianHom = Double.parseDouble(dataLine.get(1));
            final double sdHet = Double.parseDouble(dataLine.get(2));
            final int rareMin = Integer.parseInt(dataLine.get(3));
            final int rareMax = Integer.parseInt(dataLine.get(4));
            final int commonMin = Integer.parseInt(dataLine.get(5));
            final int commonMax = Integer.parseInt(dataLine.get(6));
            final int rarePass = Integer.parseInt(dataLine.get(7));
            final int rareFail = Integer.parseInt(dataLine.get(8));
            final int commonPass = Integer.parseInt(dataLine.get(9));
            final int commonFail = Integer.parseInt(dataLine.get(10));
            final double rareSingle = Double.parseDouble(dataLine.get(11));
            final double rareBoth = Double.parseDouble(dataLine.get(12));
            final double commonSingle = Double.parseDouble(dataLine.get(13));
            final double commonBoth = Double.parseDouble(dataLine.get(14));
            return new SplitReadGenotypeMetrics(
                    new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(minCount, medianHom, sdHet),
                    new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                            new SplitReadEvidenceGenotyper.CutoffResult(rareSingle, rareBoth, rarePass, rareFail, rareMin, rareMax),
                            new SplitReadEvidenceGenotyper.CutoffResult(commonSingle, commonBoth, commonPass, commonFail, commonMin, commonMax)
                    )
            );
        }
    }

    public record SplitReadGenotypeMetrics(SplitReadGenotypeParameters parameters,
                                           SplitReadGenotypeFrequencyCutoffs cutoffs) {}
}
