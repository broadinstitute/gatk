package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import htsjdk.variant.variantcontext.VariantContext;
import it.unimi.dsi.fastutil.ints.Int2DoubleArrayMap;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.tools.walkers.mutect.MutectStats;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.NaturalLogUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SomaticClusteringModel {

    private static final int MAX_INDEL_SIZE_IN_PRIOR_MAP = 10;
    private final Map<Integer, Double> logVariantPriors = new Int2DoubleArrayMap();
    private double logVariantVsArtifactPrior;
    private final OptionalDouble callableSites;

    private static final double INITIAL_HIGH_AF_WEIGHT = 0.01;
    private static final double INITIAL_BACKGROUND_WEIGHT = 0.01;

    private double REGULARIZING_PSEUDOCOUNT = 1;

    private double logHighAFWeight = Math.log(INITIAL_HIGH_AF_WEIGHT);
    private double logBackgroundWeight = Math.log(INITIAL_BACKGROUND_WEIGHT);

    private double logSparseClustersWeight = NaturalLogUtils.log1mexp(NaturalLogUtils.logSumExp(logHighAFWeight, logBackgroundWeight));

    private static final double CONCENTRATION = 0.5;
    private static final int NUM_ITERATIONS = 5;

    private static final int SEQUENCING_ERROR_INDEX = 0;
    private static final int HIGH_AF_INDEX = 1;
    private static final int BACKGROUND_INDEX = 2;
    private static final int OFFSET = 3;

    private static final BetaDistributionShape INITIAL_HIGH_AF_BETA = new BetaDistributionShape(10, 1);
    private static final BetaDistributionShape INITIAL_BACKGROUND_BETA = BetaDistributionShape.FLAT_BETA;
    private static final AlleleFractionCluster NEW_CLUSTER = new BetaBinomialCluster(BetaDistributionShape.FLAT_BETA);

    private final RandomDataGenerator rng = Utils.getRandomDataGenerator();


    final List<Datum> data = new ArrayList<>();
    List<AlleleFractionCluster> clusters;

    List<OptionalInt> clusterAssignments = new ArrayList<>();
    final List<MutableInt> clusterCounts = new ArrayList<>();
    final MutableInt totalSparseClusterCount = new MutableInt(0);

    private boolean firstPass = true;

    public SomaticClusteringModel(final M2FiltersArgumentCollection MTFAC, final List<MutectStats> mutectStats) {
        IntStream.range(-MAX_INDEL_SIZE_IN_PRIOR_MAP, MAX_INDEL_SIZE_IN_PRIOR_MAP + 1).forEach(n -> logVariantPriors.put(n, MTFAC.getLogIndelPrior()));
        logVariantPriors.put(0, MTFAC.getLogSnvPrior());

        logVariantVsArtifactPrior = MTFAC.initialLogPriorOfVariantVersusArtifact;
        callableSites = mutectStats.stream().filter(stat -> stat.getStatistic().equals(Mutect2Engine.CALLABLE_SITES_NAME))
                .mapToDouble(MutectStats::getValue).findFirst();

        clusters = new ArrayList<>();
        clusters.add(SEQUENCING_ERROR_INDEX, new SequencingError());
        clusters.add(HIGH_AF_INDEX, new BetaBinomialCluster(INITIAL_HIGH_AF_BETA));
        clusters.add(BACKGROUND_INDEX, new BetaBinomialCluster(INITIAL_BACKGROUND_BETA));
    }

    public double getLogPriorOfSomaticVariant(final VariantContext vc, final int altIndex) {
        final int indelLength = indelLength(vc, altIndex);
        return getLogPriorOfSomaticVariant(indelLength);
    }

    private double getLogPriorOfSomaticVariant(final int indelLength) {
        if (!logVariantPriors.containsKey(indelLength)) {
            logVariantPriors.put(indelLength, logVariantPriors.values().stream().mapToDouble(d -> d).min().getAsDouble());
        }

        return logVariantPriors.get(indelLength) + (indelLength == 0 ? MathUtils.LOG_ONE_THIRD : 0);
    }

    public double getLogPriorOfVariantVersusArtifact() { return logVariantVsArtifactPrior; }

    public double probabilityOfSequencingError(final Datum datum) {
        return clusterProbabilities(datum)[SEQUENCING_ERROR_INDEX];
    }

    public void record(final int[] tumorADs, final double[] tumorLogOdds, final double artifactProbability, final double nonSomaticProbability, final VariantContext vc) {
        final int totalAD = (int) MathUtils.sum(tumorADs);
        // split into one-vs-all biallelics for clustering
        for (int i = 0; i < tumorLogOdds.length; i++) {
            data.add(new Datum(tumorLogOdds[i], artifactProbability, nonSomaticProbability, tumorADs[i+1], totalAD, indelLength(vc, i)));
        }
    }

    public void learnAndClearAccumulatedData() {
        if (firstPass) {
            clusterAssignments.addAll(Collections.nCopies(data.size(), OptionalInt.empty()));
            IntStream.range(0, clusters.size()).forEach(n -> clusterCounts.add(new MutableInt(0)));
        }

        for (int iteration = 0; iteration < NUM_ITERATIONS; iteration++) {
            for (int datumIndex = 0; datumIndex < data.size(); datumIndex++) {
                final Datum datum = popDatum(datumIndex);

                //stochastically assign to non-sequencing error (sequencing error is handled within the clustering)
                if (rng.nextUniform(0,1) < datum.getNonSequencingErrorProb()) {
                    continue;
                }

                final double[] clusterPosteriors = clusterProbabilities(datum);
                final int[] indices = IntStream.range(0, clusters.size() + 1).toArray();
                final int clusterIndex = new EnumeratedIntegerDistribution(rng.getRandomGenerator(), indices, clusterPosteriors).sample();
                assignDatum(datumIndex, clusterIndex);
            }

            pruneEmptyClusters();
            final List<List<Datum>> dataByCluster = clusters.stream().map(c -> new ArrayList<Datum>()).collect(Collectors.toList());
            new IndexRange(0, clusterAssignments.size()).forEach(n -> clusterAssignments.get(n).ifPresent(c -> dataByCluster.get(c).add(data.get(n))));
            new IndexRange(0, clusters.size()).forEach(c -> clusters.get(c).learn(dataByCluster.get(c)));
            learnWeightsAndPriors();
        }

        firstPass = false;
        data.clear();
    }

    private void pruneEmptyClusters() {
        final Map<Integer, Integer> oldToNewClusterIndices = IntStream.range(0, OFFSET).boxed()
                .collect(Collectors.toMap(n ->n, n-> n));

        int newIndex = OFFSET;
        for (int oldIndex = OFFSET; oldIndex < clusters.size(); oldIndex++) {
            if (clusterCounts.get(oldIndex).getValue() > 0) {
                oldToNewClusterIndices.put(oldIndex, newIndex);

                if (newIndex != oldIndex) {
                    clusters.set(newIndex, clusters.get(oldIndex));
                    clusterCounts.set(newIndex, clusterCounts.get(oldIndex));
                }
                newIndex++;
            }
        }
        // at this point newIndex equals the size of new clusters, and possibly some indices past this must be deleted
        clusters.subList(newIndex,clusters.size()).clear();
        clusterCounts.subList(newIndex,clusterCounts.size()).clear();

        clusterAssignments = clusterAssignments.stream()
                .map(a -> a.isPresent() ? OptionalInt.of(oldToNewClusterIndices.get(a.getAsInt())) : a)
                .collect(Collectors.toList());
    }

    // emission likelihood of given alt count given that a variant is somatic
    public double logLikelihoodGivenSomatic(final int totalCount, final int altCount) {
        final double[] logClusterLikelihoods = IntStream.range(0, clusters.size())
                .filter(c -> c != SEQUENCING_ERROR_INDEX)
                .mapToDouble(c -> {
                    final double logLikelihood = clusters.get(c).logLikelihood(totalCount, altCount);
                    if (c == HIGH_AF_INDEX) {
                        return logHighAFWeight + logLikelihood;
                    } else if (c == BACKGROUND_INDEX) {
                        return logBackgroundWeight + logLikelihood;
                    } else {   // sparse cluster
                        return logSparseClustersWeight + logCRPWeight(c) + logLikelihood;
                    }
                }).toArray();
        return NaturalLogUtils.logSumExp(logClusterLikelihoods);
    }

    private double[] clusterProbabilities(final Datum datum) {
        final double logVariantPrior = getLogPriorOfSomaticVariant(datum.getIndelLength());
        final double logNoVariantPrior = NaturalLogUtils.log1mexp(logVariantPrior);

        final double[] logClusterPosteriors = new IndexRange(0, clusters.size() + 1).mapToDouble(c -> {
            final double logLikelihood = c < clusters.size() ? clusters.get(c).logLikelihood(datum) :
                    NEW_CLUSTER.logLikelihood(datum);
            if (c == SEQUENCING_ERROR_INDEX) {
                return logNoVariantPrior + logLikelihood;
            } else if (c == HIGH_AF_INDEX) {
                return logVariantPrior + logHighAFWeight + logLikelihood;
            } else if (c == BACKGROUND_INDEX) {
                return logVariantPrior + logBackgroundWeight + logLikelihood;
            } else if (c < clusters.size()) {   // existing sparse cluster
                return logVariantPrior + logSparseClustersWeight + logCRPWeight(c)
                        + logLikelihood;
            } else {    // new sparse cluster
                return logVariantPrior + logSparseClustersWeight + logCRPWeight(c)
                        + logLikelihood;
            }
        });

        return NaturalLogUtils.normalizeLog(logClusterPosteriors, false, false);
    }

    private double logCRPWeight(final int clusterIndex) {
        Utils.validate(clusterIndex >= OFFSET, "Chinese restaurant process does not apply to error, high-AF, and backgorund clusters");
        final double numerator = clusterIndex == clusters.size() ? CONCENTRATION : clusterCounts.get(clusterIndex).getValue();
        final double denominator = totalSparseClusterCount.getValue() + CONCENTRATION;
        return Math.log(numerator / denominator);

    }

    private Datum popDatum(final int datumIndex) {
        // pop datum off its cluster and decrement the sparse cluster count if appropriate
        clusterAssignments.get(datumIndex).ifPresent(c -> {
            clusterCounts.get(c).decrement();
            if (OFFSET <= c) {
                totalSparseClusterCount.decrement();
            }
        });
        clusterAssignments.set(datumIndex, OptionalInt.empty());
        return data.get(datumIndex);
    }

    private void assignDatum(int datumIndex, int clusterIndex) {
        final Datum datum = data.get(datumIndex);

        // make new cluster
        if (clusterIndex == clusters.size()) {
            final double newClusterAlleleFraction = new BetaDistribution(rng.getRandomGenerator(), datum.getAltCount() + 1, datum.getTotalCount() - datum.getAltCount() + 1).sample();
            clusters.add(new BinomialCluster(newClusterAlleleFraction));
            clusterCounts.add(new MutableInt(0));
        }

        if (OFFSET <= clusterIndex) {
            totalSparseClusterCount.increment();
        }

        clusterAssignments.set(datumIndex, OptionalInt.of(clusterIndex));
        clusterCounts.get(clusterIndex).increment();
    }

    private void learnWeightsAndPriors() {
        final double totalVariants = clusterCounts.get(HIGH_AF_INDEX).getValue() + clusterCounts.get(BACKGROUND_INDEX).getValue()
                + totalSparseClusterCount.getValue() + REGULARIZING_PSEUDOCOUNT * 3; // high-AF, background, and sparse each get pseudocounts

        logHighAFWeight = Math.log((REGULARIZING_PSEUDOCOUNT + clusterCounts.get(HIGH_AF_INDEX).getValue()) / totalVariants);
        logBackgroundWeight = Math.log((REGULARIZING_PSEUDOCOUNT + clusterCounts.get(BACKGROUND_INDEX).getValue()) / totalVariants);
        logSparseClustersWeight = Math.log((REGULARIZING_PSEUDOCOUNT + totalSparseClusterCount.getValue()) /totalVariants);

        final Map<Integer, Long> variantCountsByIndelLength = IntStream.range(0, data.size())
                .filter( n -> clusterAssignments.get(n).orElse(SEQUENCING_ERROR_INDEX) != SEQUENCING_ERROR_INDEX)
                .map(n -> data.get(n).getIndelLength())
                .boxed()
                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting() ));


        final double technicalArtifactCount = data.stream().mapToDouble(Datum::getArtifactProb).sum();

        if (callableSites.isPresent()) {
            IntStream.range(-MAX_INDEL_SIZE_IN_PRIOR_MAP, MAX_INDEL_SIZE_IN_PRIOR_MAP + 1).forEach(n -> {
                final double empiricalRatio = variantCountsByIndelLength.getOrDefault(n, 0L) / callableSites.getAsDouble();
                logVariantPriors.put(n, Math.log(Math.max(empiricalRatio, n == 0 ? 1.0e-8 : 1.0e-9)));
            });
        }
        final long variantCount = variantCountsByIndelLength.values().stream().mapToLong(n -> n).sum();
        logVariantVsArtifactPrior = Math.log((variantCount + REGULARIZING_PSEUDOCOUNT) / (variantCount + technicalArtifactCount + REGULARIZING_PSEUDOCOUNT * 2));
    }

    public List<Pair<String, String>> clusteringMetadata() {
        final List<Pair<String, String>> result = new ArrayList<>();
        IntStream.range(-MAX_INDEL_SIZE_IN_PRIOR_MAP, MAX_INDEL_SIZE_IN_PRIOR_MAP + 1).forEach(n -> {
            final double logPrior = logVariantPriors.get(n);
            final String type = n == 0 ? "SNV" :
                    (n < 0 ? "deletion" : "insertion") + " of length " + Math.abs(n);
            result.add(ImmutablePair.of("Ln prior of " + type, Double.toString(logPrior)));
        });

        result.add(ImmutablePair.of("High-AF beta-binomial cluster",
                String.format("weight = %.4f, %s", Math.exp(logHighAFWeight), clusters.get(HIGH_AF_INDEX).toString())));
        result.add(ImmutablePair.of("Background beta-binomial cluster",
                String.format("weight = %.4f, %s", Math.exp(logBackgroundWeight), clusters.get(BACKGROUND_INDEX).toString())));

        final MutableInt clusterIndex = new MutableInt(1);
        IntStream.range(OFFSET, clusters.size()).boxed()
                .sorted(Comparator.comparingDouble(c -> -logCRPWeight(c)))
                .forEach(c -> result.add(ImmutablePair.of("Binomial cluster " + clusterIndex.toString(),
                        String.format("weight = %.4f, %s", Math.exp(logCRPWeight(c)), clusters.get(c).toString()))));
        return result;
    }

    public static int indelLength(final VariantContext vc, final int altIndex) {
        return vc.getAlternateAllele(altIndex).length() - vc.getReference().length();
    }
}
