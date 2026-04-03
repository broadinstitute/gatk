package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.aggregation.DiscordantPairEvidenceTester;
import org.broadinstitute.hellbender.tools.sv.aggregation.EvidenceStatUtils;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.*;
import java.util.function.Function;
import java.util.stream.DoubleStream;

public class DiscordantPairEvidenceGenotyper {

    private final Integer minSize;
    private final int trainingMinCount;
    private final double targetCoverage;
    private final int maxQuality;
    private final Map<String, Double> sampleCoverageMap;
    private final Map<String, FirstPassResult> firstPassCounts;
    private Double hetCutoff = null;
    private final List<Double> hetCounts;
    private final List<Double> homCounts;
    private boolean firstPassMade = false;
    private boolean secondPassMade = false;

    private final Map<String, SimpleInterval> variantIntervals;
    private OverlapDetector<NamedSimpleInterval> overlapDetector = null;

    private static final Median MEDIAN = new Median();
    private static final NormalDistribution Z_DISTRIBUTION = new NormalDistribution();

    private static final Set<Integer> HET_COPY_STATES = Set.of(1, 3);

    private static class NamedSimpleInterval implements Locatable {
        private final String name;
        private final SimpleInterval interval;

        public NamedSimpleInterval(final String name, final SimpleInterval interval) {
            this.name = name;
            this.interval = interval;
        }

        @Override
        public String getContig() {
            return interval.getContig();
        }

        @Override
        public int getStart() {
            return interval.getStart();
        }

        @Override
        public int getEnd() {
            return interval.getEnd();
        }

        @Override
        public boolean equals(Object o) {
            if (o == null || getClass() != o.getClass()) return false;
            NamedSimpleInterval that = (NamedSimpleInterval) o;
            return Objects.equals(name, that.name) && Objects.equals(interval, that.interval);
        }

        @Override
        public int hashCode() {
            return Objects.hash(name, interval);
        }
    }

    public DiscordantPairEvidenceGenotyper(final Map<String, Double> sampleCoverageMap, final double qualityCutoff, final Integer minSize, final double targetCoverage, final int maxQuality) {
        this.sampleCoverageMap = Utils.nonNull(sampleCoverageMap);
        this.trainingMinCount = computeCountCutoff(qualityCutoff);
        Utils.validateArg(maxQuality > 0, "Maximum quality must be greater than zero");
        this.maxQuality = maxQuality;
        this.minSize = minSize;
        this.targetCoverage = targetCoverage;
        this.firstPassCounts = new HashMap<>();
        this.hetCounts = new ArrayList<>();
        this.homCounts = new ArrayList<>();
        this.variantIntervals = new HashMap<>();
    }

    public void registerVariantForOverlapCheck(final SVCallRecord record) {
        if ((record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL || record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) && !record.isDepthOnly()) {
            final SimpleInterval interval = new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB());
            if (variantIntervals.containsKey(record.getId())) {
                throw new IllegalArgumentException("Duplicate variant ID: " + record.getId());
            }
            variantIntervals.put(record.getId(), interval);
        }
    }

    public void aggregateOverlapCheckIntervals() {
        final List<NamedSimpleInterval> namedIntervals = variantIntervals.entrySet().stream().map(e -> new NamedSimpleInterval(e.getKey(), e.getValue())).toList();
        overlapDetector = OverlapDetector.create(namedIntervals);
        variantIntervals.clear();
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

    public void addFirstPass(final SVCallRecord record, final List<DiscordantPairEvidence> evidence, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotypeResult, final List<String> depthGenotypeSamples) {
        Utils.nonNull(evidence);
        Utils.validate(!firstPassMade, "First pass already made");
        // TODO check all samples are in coverage map
        final Map<String, Double> counts = normalizeCounts(evidence, trainingMinCount);
        firstPassCounts.put(record.getId(), new FirstPassResult(counts, depthGenotypeResult, depthGenotypeSamples));
    }

    public boolean isTrainingRecord(final SVCallRecord record) {
        return firstPassCounts.containsKey(record.getId());
    }

    // 0-counts are not added to map
    private Map<String, Double> normalizeCounts(final List<DiscordantPairEvidence> evidence, final int threshold) {
        final Map<String, Integer> countsMap = DiscordantPairEvidenceTester.countEvidence(evidence);
        final Map<String, Double> counts = new HashMap<>();
        for (final String sample : countsMap.keySet()) {
            final int count = countsMap.getOrDefault(sample, 0);
            if (count > 0) {
                final double sampleCoverage = sampleCoverageMap.getOrDefault(sample, 0.);
                final double normalizedCount = EvidenceStatUtils.getNormalizedCount(count, sampleCoverage, targetCoverage);
                if (normalizedCount >= threshold) {
                    counts.put(sample, normalizedCount);
                }
            }
        }
        return counts;
    }

    public void finalizeFirstPass() {
        Utils.validate(!firstPassCounts.isEmpty(), "No discordant pair counts after first pass");
        Utils.validate(!firstPassMade, "First pass has already been made");
        final double[] hetCounts = firstPassCounts.values().stream().map(c -> c.getCounts(HET_COPY_STATES)).flatMapToDouble(DoubleStream::of).toArray();
        final double hetMedian = MEDIAN.evaluate(hetCounts);
        final double[] deviations = DoubleStream.of(hetCounts).map(d -> Math.abs(d - hetMedian)).toArray();
        final double hetMad = MEDIAN.evaluate(deviations);
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

    public DiscordantPairGenotypeParameters finalizeSecondPass() {
        Utils.validate(!homCounts.isEmpty(), "No discordant pair counts after second pass");
        Utils.validate(!hetCounts.isEmpty(), "No discordant pair counts after second pass");
        Utils.validate(!secondPassMade, "Second pass has already been made");
        final double[] homArr = homCounts.stream().mapToDouble(Double::valueOf).toArray();
        final double[] hetArr = hetCounts.stream().mapToDouble(Double::valueOf).toArray();
        final double homMedian = MEDIAN.evaluate(homArr);
        final double hetMedian = MEDIAN.evaluate(hetArr);
        final double hetMadValue = MEDIAN.evaluate(DoubleStream.of(hetArr).map(d -> Math.abs(d - hetMedian)).toArray());
        final double sdHet = 1.645 * 1.4826 * hetMadValue;
        secondPassMade = true;
        // Free training accumulation lists that are no longer needed
        hetCounts.clear();
        homCounts.clear();
        return new DiscordantPairGenotypeParameters(trainingMinCount, homMedian, sdHet);
    }

    /**
     * Clears the first-pass training data (per-variant normalized counts and depth genotypes).
     * Call this after all consumers of {@link #isTrainingRecord} have finished (i.e. after SR
     * first pass finalization) to free significant heap memory.
     */
    public void clearTrainingData() {
        firstPassCounts.clear();
    }

    public DiscordantPairGenotypeResult genotype(final SVCallRecord record, final List<DiscordantPairEvidence> evidence,
                                                 final DiscordantPairGenotypeParameters parameters, final List<String> samples) {
        final Map<String, Double> counts = normalizeCounts(evidence, 0);
        final int[] genotypes = new int[samples.size()];
        final int[] genotypeQuals = new int[samples.size()];
        final List<Double> nonRefCounts = new ArrayList<>(samples.size());
        final double normalizationGenotype = maxQuality / (-10. * Math.log10(1 - Z_DISTRIBUTION.cumulativeProbability(0.5 * parameters.medianHom() / parameters.sdHet())) - 3.0103);
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final double discordantPairCount = counts.getOrDefault(sample, 0.);
            if (discordantPairCount < parameters.minCount) {
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
                    genotypeQuals[i] = (int) Math.round(-10. * Math.log10(1. - dist.cumulativeProbability(0)) * normalizationGenotype);
                }
            } else {
                final double z0 = discordantPairCount - genotypes[i] * parameters.medianHom * 0.5;
                final double z1 = z0 + parameters.medianHom * 0.5;
                final double z2 = z0 - parameters.medianHom * 0.5;
                final double q0 = getGenotypeLikelihood(z0, parameters.sdHet);
                final double q1 = getGenotypeLikelihood(z1, parameters.sdHet);
                final double q2 = getGenotypeLikelihood(z2, parameters.sdHet);
                genotypeQuals[i] = (int) Math.round((Math.min(q1, q2) - q0) * normalizationGenotype);
            }
            genotypeQuals[i] = Math.max(Math.min(genotypeQuals[i], maxQuality), 1);
            if (genotypes[i] != 0) {
                nonRefCounts.add(discordantPairCount);
            }
        }
        final PoissonDistribution poissonDistribution = new PoissonDistribution(0.5 * parameters.medianHom());
        final double normalizationVariant = maxQuality / (-10. * Math.log10(poissonDistribution.cumulativeProbability(0)));
        final double variantQual;
        if (nonRefCounts.isEmpty()) {
            variantQual = 0;
        } else {
            final double medianCarrierCount = MEDIAN.evaluate(nonRefCounts.stream().mapToDouble(Double::valueOf).toArray());
            final PoissonDistribution variantQualDist = new PoissonDistribution(medianCarrierCount);
            variantQual = Math.min(-10. * Math.log10(variantQualDist.cumulativeProbability(0)) * normalizationVariant, maxQuality);
        }
        return new DiscordantPairGenotypeResult(genotypes, genotypeQuals, variantQual);
    }

    private static double getGenotypeLikelihood(final double z, final double sdHet) {
        return -10. * Math.log10(1. - Z_DISTRIBUTION.cumulativeProbability(Math.abs(z) / sdHet));
    }

    public boolean trainableRecord(final SVCallRecord record,
                                   final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype,
                                   final SVStratificationEngine exclusionEngine) {
        if (record.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.DEL && record.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            return false;
        }
        if (record.isDepthOnly()) {
            return false;
        }
        if (minSize != null && (record.getLength() == null || record.getLength() < minSize)) {
            return false;
        }
        if (record.isDepthOnly()) {
            return false;
        }
        if (exclusionEngine != null && !exclusionEngine.getMatches(record, 0, 1, 1).isEmpty()) {
            return false;
        }
        final SimpleInterval interval = new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB());
        if (overlapDetector != null && overlapDetector.getOverlaps(interval).size() > 1) {
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
        private final Map<String, Integer> depthGenotypes;

        public FirstPassResult(final Map<String, Double> counts, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype, final List<String> depthGenotypeSamples) {
            this.counts = counts;
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

    public record DiscordantPairGenotypeResult(int[] genotypes, int[] genotypeQuals, double variantQual) {}
    public record DiscordantPairGenotypeParameters(double minCount, double medianHom, double sdHet) {}

    public static final class DiscordantPairTableParser {

        private static final String MIN_COUNT_COLUMN = "pe_count";
        private static final String MEDIAN_HOM_COLUMN = "median_hom";
        private static final String SD_HET_COLUMN = "sd_het";
        public static final TableColumnCollection CUTOFFS_COLUMNS = new TableColumnCollection(Arrays.asList(MIN_COUNT_COLUMN, MEDIAN_HOM_COLUMN, SD_HET_COLUMN));

        public void composeCutoffsLine(final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters stats,
                                        final DataLine dataLine) {
            dataLine.append(stats.minCount());
            dataLine.append(stats.medianHom());
            dataLine.append(stats.sdHet());
        }

        public Function<DataLine, DiscordantPairGenotypeParameters> tableParser(TableColumnCollection columns, Function<String, RuntimeException> exceptionFactory) {
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

        public DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters parseTableLine(final DataLine dataLine) {
            final double minCount = Double.parseDouble(dataLine.get(0));
            final double medianHom = Double.parseDouble(dataLine.get(1));
            final double sdHet = Double.parseDouble(dataLine.get(2));
            return new DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters(minCount, medianHom, sdHet);
        }
    }
}
