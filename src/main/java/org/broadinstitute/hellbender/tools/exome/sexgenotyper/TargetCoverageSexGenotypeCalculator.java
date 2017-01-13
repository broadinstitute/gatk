package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.Sets;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Infers sex genotypes of samples from their read counts.
 *
 * Refer to the documentation of {@link TargetCoverageSexGenotyper} for additional information.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class TargetCoverageSexGenotypeCalculator {

    public static final Logger logger = LogManager.getLogger(TargetCoverageSexGenotypeCalculator.class);

    private final ReadCountCollection processedReadCounts;

    /**
     * If the germline ploidy of a sex genotype class is 0, it will be replaced by the following value
     */
    private final double baselineMappingErrorProbability;

    private final List<Target> autosomalTargetList;
    private final List<Target> allosomalTargetList;

    private final Set<String> sexGenotypeIdentifiers;

    /**
     * Autosomal target plodies
     */
    private final int[] autosomalTargetPloidies;
    /**
     * Map from genotype names to the list of target germline ploidies in the same order as {@code allosomalTargetList}
     */
    private final Map<String, int[]> allosomalTargetPloidies;

    /**
     * Public constructor.
     *
     * @param rawReadCounts raw read count collection
     * @param targetList list of targets
     * @param contigPloidyAnnotsList list of contig germline ploidy annotations
     * @param baselineMappingErrorProbability typical mapping error probability
     */
    public TargetCoverageSexGenotypeCalculator(@Nonnull final ReadCountCollection rawReadCounts,
                                               @Nonnull final List<Target> targetList,
                                               @Nonnull final List<ContigGermlinePloidyAnnotation> contigPloidyAnnotsList,
                                               final double baselineMappingErrorProbability) {
        this.baselineMappingErrorProbability = ParamUtils.inRange(baselineMappingErrorProbability,
                0, 1, "The baseline mapping error probability must be a number in [0, 1] interval");

        /* check if targetList is a subset of targets in the raw read counts */
        if (!Sets.difference(new HashSet<>(targetList), new HashSet<>(rawReadCounts.targets())).isEmpty()) {
            throw new UserException.BadInput("The provided target list must be a subset of targets in the provided" +
                    " read count collection.");
        }

        final ImmutablePair<ReadCountCollection, List<Target>> processedReadCountsAndTargets =
                processReadCountsAndTargets(rawReadCounts, targetList);
        processedReadCounts = processedReadCountsAndTargets.left;
        final List<Target> processedTargetList = processedReadCountsAndTargets.right;

        /* annotate targets with germline ploidy data */
        final GermlinePloidyAnnotatedTargetCollection ploidyAnnots =
                new GermlinePloidyAnnotatedTargetCollection(contigPloidyAnnotsList, processedTargetList);

        /* populate lists and maps */
        autosomalTargetList = ploidyAnnots.getAutosomalTargetList();
        allosomalTargetList = ploidyAnnots.getAllosomalTargetList();
        sexGenotypeIdentifiers = contigPloidyAnnotsList.get(0).getGenotypesSet();

        if (allosomalTargetList.isEmpty()) {
            throw new UserException.BadInput("No allosomal targets could be identified. Please ensure (1) the " +
                    "read count collections covers some allosomal targets, and (2) allosomal targets are properly " +
                    "annotated");
        }

        if (autosomalTargetList.isEmpty()) {
            throw new UserException.BadInput("No autosomal targets could be identified. Please ensure (1) the " +
                    "read count collections covers some autosomal targets, and (2) autosomal targets are properly " +
                    "annotated");
        }

        logger.info("Number of autosomal targets: " + autosomalTargetList.size());
        logger.info("Number of allosomal targets: " + allosomalTargetList.size());

        /* autosomal targets have the same ploidy for all ploidy classes (this is asserted earlier); so pick the first one */
        final String somePloidyTag = sexGenotypeIdentifiers.iterator().next();
        autosomalTargetPloidies = autosomalTargetList.stream()
                .mapToInt(t -> ploidyAnnots.getTargetGermlinePloidyByGenotype(t, somePloidyTag)).toArray();

        /* allosomal targets may have different ploidies for different ploidy classes */
        allosomalTargetPloidies = Collections.unmodifiableMap(sexGenotypeIdentifiers.stream()
                .map(sexGenotype -> ImmutablePair.of(
                        sexGenotype,
                        allosomalTargetList.stream()
                                .mapToInt(t -> ploidyAnnots.getTargetGermlinePloidyByGenotype(t, sexGenotype)).toArray()))
                .collect(Collectors.toMap(ImmutablePair::getLeft, ImmutablePair::getRight)));
    }

    /**
     * Infer sex genotypes for all samples
     *
     * @return an instance of {@link SexGenotypeDataCollection}
     */
    public SexGenotypeDataCollection inferSexGenotypes() {
        final SexGenotypeDataCollection col = new SexGenotypeDataCollection();
        IntStream.range(0, processedReadCounts.columnNames().size())
                .forEach(si -> col.add(calculateSexGenotypeData(si)));
        return col;
    }

    /**
     * Processes raw read counts and targets:
     * <dl>
     *     <dt> If more than one sample is present in the collection, filters out fully uncovered targets
     *     from read counts and removes the uncovered targets from the target list</dt>
     *
     *     <dt> Otherwise, does nothing and warns the user
     *     </dt>
     * </dl>
     *
     * @param rawReadCounts raw read count collection
     * @param targetList user provided target list
     * @return pair of processed read counts and targets
     */
    private ImmutablePair<ReadCountCollection, List<Target>> processReadCountsAndTargets(
            @Nonnull final ReadCountCollection rawReadCounts,
            @Nonnull final List<Target> targetList) {
        final ReadCountCollection finalReadCounts;
        final List<Target> finalTargetList;

        /* remove totally uncovered targets */
        if (rawReadCounts.columnNames().size() > 1) {
            finalReadCounts = ReadCountCollectionUtils.removeTotallyUncoveredTargets(rawReadCounts, logger);
            final Set<Target> targetSetFromProcessedReadCounts = new HashSet<>(finalReadCounts.targets());
            finalTargetList = targetList.stream()
                    .filter(targetSetFromProcessedReadCounts::contains)
                    .collect(Collectors.toList());
        } else {
            final long numUncoveredTargets = rawReadCounts.records().stream()
                    .filter(rec -> (int)rec.getDouble(0) == 0).count();
            final long numAllTargets = rawReadCounts.targets().size();
            logger.info("Since only one sample is given for genotyping, the user is responsible for asserting" +
                    " the aptitude of targets. Fully uncovered (irrelevant) targets can not be automatically" +
                    " identified (total targets: " + numAllTargets + ", uncovered targets: " + numUncoveredTargets + ")");
            finalReadCounts = rawReadCounts;
            finalTargetList = targetList;
        }
        return ImmutablePair.of(finalReadCounts, finalTargetList);
    }

    /**
     * Estimates read depth per target per homolog for a given sample index in the collection.
     *
     * @param sampleIndex integer index of the sample in the read count collection
     * @return read depth per target per homolog
     */
    private double getSampleReadDepthFromAutosomalTargets(final int sampleIndex) {
        final double[] readCounts = processedReadCounts.getColumnOnSpecifiedTargets(sampleIndex,
                autosomalTargetList, false);
        final double[] readCountsNormalizedByPloidy = IntStream.range(0, readCounts.length)
                .mapToDouble(i -> readCounts[i] / (double)autosomalTargetPloidies[i])
                .toArray();
        return new Median().evaluate(readCountsNormalizedByPloidy);

        /*
         * @implNote Old code:
         *
         * (assume copy ratios and multiplicative biases are 1.0 at this stage)
         *
         * final int[] mask = IntStream.range(0, readCounts.length).map(i -> 1).toArray();
         * final double[] copyRatio = IntStream.range(0, readCounts.length).mapToDouble(i -> 1.0).toArray();
         * final double[] multBias = IntStream.range(0, readCounts.length).mapToDouble(i -> 1.0).toArray();
         *
         * return CoverageModelUtils.estimateReadDepthFromPoissonModel(readCounts, multBias, autosomalTargetPloidies,
         *      copyRatio, mask);
         *
         * @implNote This is prone to over-estimation because the Gaussian approximation of the Poisson distribution,
         * which is currently used in CoverageModelUtils.estimateReadDepthFromPoissonModel breaks down
         * for outlier read counts (if n >> depth). For the time being, we use medians which is a robust
         * statistic. In the future, we must find a better approximation for the Poisson distribution which
         * is both analytically tractable and robust.
         */
    }

    /**
     * Calculates the likelihood of different sex genotypes for a given sample index
     *
     * @param sampleIndex sample index
     * @return an instance of {@link SexGenotypeData}
     */
    private SexGenotypeData calculateSexGenotypeData(final int sampleIndex) {
        final int[] allosomalReadCounts = Arrays.stream(processedReadCounts.getColumnOnSpecifiedTargets(sampleIndex,
                allosomalTargetList, true)).mapToInt(n -> (int)n).toArray();
        final double readDepth = getSampleReadDepthFromAutosomalTargets(sampleIndex);

        logger.debug("Read depth for " + processedReadCounts.columnNames().get(sampleIndex) + ": " + readDepth);
        final List<Double> logLikelihoods = new ArrayList<>();
        final List<String> sexGenotypesList = new ArrayList<>();
        final int numAllosomalTargets = allosomalTargetList.size();

        for (final String genotypeName : sexGenotypeIdentifiers) {
            /* calculate log likelihood */
            final int[] currentAllosomalTargetPloidies = allosomalTargetPloidies.get(genotypeName);
            final double[] poissonParameters = IntStream.range(0, numAllosomalTargets)
                    .mapToDouble(i -> readDepth *
                            (currentAllosomalTargetPloidies[i] > 0 ? currentAllosomalTargetPloidies[i]
                                    : baselineMappingErrorProbability))
                    .toArray();
            final double currentLogLikelihood = IntStream.range(0, numAllosomalTargets)
                    .mapToDouble(i -> {
                        final PoissonDistribution pois = new PoissonDistribution(poissonParameters[i]);
                        return pois.logProbability(allosomalReadCounts[i]);
                    }).sum();
            sexGenotypesList.add(genotypeName);
            logLikelihoods.add(currentLogLikelihood);
        }

        /* infer the most likely sex genotype */
        final Integer[] indices = new Integer[sexGenotypesList.size()];
        IntStream.range(0, sexGenotypesList.size()).forEach(i -> indices[i] = i);
        Arrays.sort(indices, (li, ri) -> Double.compare(logLikelihoods.get(ri), logLikelihoods.get(li)));

        return new SexGenotypeData(processedReadCounts.columnNames().get(sampleIndex),
                sexGenotypesList.get(indices[0]), sexGenotypesList,
                logLikelihoods.stream().mapToDouble(d -> d).toArray());
    }

    /**
     * Returns a list of all possible sex genotypes (their string identifers)
     * @return set of sex genotype string identifiers
     */
    public Set<String> getSexGenotypeIdentifiers() {
        return sexGenotypeIdentifiers;
    }

    /**
     * Returns the ploidies of autosomal targets in the same order as targets in the input read count collection
     * @return an integer array of autosomal target ploidies
     */
    public int[] getAutosomalTargetGermlinePloidies() {
        return autosomalTargetPloidies;
    }

    /**
     * Returns a map from sex genotype string identifiers to the corresponding ploidies of allosomal targets in the
     * same order as allosomal targets in the input read count collection
     * @return a map from sex genotype string identifiers to an integer array of allosomal target ploidies
     */
    public Map<String, int[]> getAllosomalTargetGermlinePloidiesMap() {
        return allosomalTargetPloidies;
    }

    /**
     * Returns a list of autosomal targets fetched from the input read count collection in the same order
     * @return list of autosomal targets
     */
    public List<Target> getAutosomalTargetList() {
        return autosomalTargetList;
    }

    /**
     * Returns a list of allosomal targets fetched from the input read count collection in the same order
     * @return list of allosomal targets
     */
    public List<Target> getAllosomalTargetList() {
        return allosomalTargetList;
    }
}
