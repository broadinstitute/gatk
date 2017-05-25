package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import avro.shaded.com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetAnnotation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
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

    private final ReadCountCollection readCounts;

    /**
     * If the germline ploidy of a sex genotype class is 0, it will be replaced by the following value
     */
    private final double baselineMappingErrorProbability;

    private final List<Target> autosomalTargetList;

    private final List<Target> allosomalTargetList;

    private final Set<String> sexGenotypeIdentifiers;

    /**
     * Autosomal target ploidies in the same order as {@link #autosomalTargetList}
     */
    private final int[] autosomalTargetPloidies;

    /**
     * Autosomal target masks in the same order as {@link #autosomalTargetList}.
     * Targets with mask = 0 will be ignored.
     */
    private final int[] autosomalTargetMasks;

    /**
     * Mean ploidy of autosomal targets
     */
    private final double meanAutosomalPloidy;

    /**
     * Autosomal target multiplicative biases in the same order as {@link #autosomalTargetList}
     */
    private final double[] autosomalTargetBiases;

    /**
     * Map from genotype string identifiers to the germline ploidies of allosomal targets.
     * The array entries are in the same order as {@link #allosomalTargetList}
     */
    private final Map<String, int[]> allosomalTargetPloidies;

    /**
     * Autosomal target multiplicative biases in the same order as {@link #allosomalTargetList}
     */
    private final double[] allosomalTargetBiases;

    /**
     * Allosomal target masks in the same order as {@link #allosomalTargetList}.
     * Targets with mask = 0 will be ignored.
     */
    private final int[] allosomalTargetMasks;

    /**
     * If bait count per target is lower than the following value, that target will be excluded
     */
    public static final double MIN_ALLOWED_BAIT_COUNT_PER_TARGET = 1.0;

    /**
     * Public constructor.
     *
     * @param readCounts raw read count collection
     * @param genotypingTargets list of targets to consider for genotyping; target annotations, if present, will be used
     *                          for estimating biases
     * @param contigPloidyAnnotsList list of contig germline ploidy annotations
     * @param genotypingIntervals (nullable) list of intervals to consider for genotyping. If it is null, all targets
     *                            will be considered; otherwise, only targets overlapping with these intervals will
     *                            be used
     * @param baselineMappingErrorProbability typical mapping error probability
     */
    public TargetCoverageSexGenotypeCalculator(@Nonnull final ReadCountCollection readCounts,
                                               @Nonnull final List<Target> genotypingTargets,
                                               @Nonnull final List<ContigGermlinePloidyAnnotation> contigPloidyAnnotsList,
                                               @Nullable final List<SimpleInterval> genotypingIntervals,
                                               final double baselineMappingErrorProbability) {
        this.readCounts = Utils.nonNull(readCounts);
        this.baselineMappingErrorProbability = ParamUtils.inRange(baselineMappingErrorProbability,
                0, 1, "The baseline mapping error probability must be a number in [0, 1] interval");
        Utils.nonNull(genotypingTargets);
        Utils.nonNull(contigPloidyAnnotsList);

        /* keep genotyping targets that exist in the read count collection */
        final Set<Target> readCountTargets = new HashSet<>(readCounts.targets());
        final List<Target> mutualTargets = genotypingTargets.stream()
                .filter(readCountTargets::contains)
                .collect(Collectors.toList());

        /* annotate targets with germline ploidy data */
        final GermlinePloidyAnnotatedTargetCollection ploidyAnnots =
                new GermlinePloidyAnnotatedTargetCollection(contigPloidyAnnotsList, mutualTargets);

        /* populate lists and maps */
        autosomalTargetList = ploidyAnnots.getAutosomalTargetList();
        allosomalTargetList = ploidyAnnots.getAllosomalTargetList();
        sexGenotypeIdentifiers = contigPloidyAnnotsList.get(0).getGenotypesSet();

        /* autosomal targets have the same ploidy for all ploidy classes (this is asserted earlier); so pick the first one */
        autosomalTargetPloidies = getPloidyArray(sexGenotypeIdentifiers.iterator().next(), autosomalTargetList, ploidyAnnots);

        /* calculate mean autosomal ploidy */
        meanAutosomalPloidy = getAutosomalMeanPloidy(autosomalTargetList, autosomalTargetPloidies);

        /* allosomal targets may have different ploidies for different ploidy classes */
        allosomalTargetPloidies = getAllosomalPloidyMap(allosomalTargetList, sexGenotypeIdentifiers, ploidyAnnots);

        /* initialize the multiplicative biases and masks */
        autosomalTargetBiases = new double[autosomalTargetList.size()];
        Arrays.fill(autosomalTargetBiases, 1.0);
        autosomalTargetMasks = new int[autosomalTargetList.size()];
        Arrays.fill(autosomalTargetMasks, 1);
        allosomalTargetBiases = new double[allosomalTargetList.size()];
        Arrays.fill(allosomalTargetBiases, 1.0);
        allosomalTargetMasks = new int[allosomalTargetList.size()];
        Arrays.fill(allosomalTargetMasks, 1);

        /* take into account bait count bias */
        if (processBaitCountAnnotations(autosomalTargetList, autosomalTargetBiases, autosomalTargetMasks, "autosomal")) {
            processBaitCountAnnotations(allosomalTargetList, allosomalTargetBiases, allosomalTargetMasks, "allosomal");
        }

        /* mask totally uncovered targets */
        maskTotallyUncoveredTargets(autosomalTargetList, autosomalTargetMasks, readCounts, "autosomal");
        maskTotallyUncoveredTargets(allosomalTargetList, allosomalTargetMasks, readCounts, "allosomal");

        /* mask targets that do not overlap with genotyping intervals */
        maskTargetsNotOverlappingWithGenotypingIntervals(autosomalTargetList, autosomalTargetMasks, genotypingIntervals, "autosomal");
        maskTargetsNotOverlappingWithGenotypingIntervals(allosomalTargetList, allosomalTargetMasks, genotypingIntervals, "allosomal");

        /* assert that the generated mask and bias is valid */
        assertBiasAndMaskValidity(autosomalTargetBiases, autosomalTargetMasks);
        assertBiasAndMaskValidity(allosomalTargetBiases, allosomalTargetMasks);

        /* assert that we have enough unmasked autosomal and allosomal targets */
        assertSomeUnmaskedTargetsExist(autosomalTargetMasks, allosomalTargetMasks);

        /* normalize mean biases to 1.0 */
        normalizeTargetBiases(autosomalTargetBiases, autosomalTargetMasks, allosomalTargetBiases, allosomalTargetMasks);

        logger.info("Final number of autosomal targets to be used (contig, count): " +
                getTargetCountPerContigStatString(autosomalTargetList, autosomalTargetMasks));
        logger.info("Final number of allosomal targets to be used (contig, count): " +
                getTargetCountPerContigStatString(allosomalTargetList, allosomalTargetMasks));
    }

    private static void normalizeTargetBiases(double[] autosomalTargetBiases, int[] autosomalTargetMasks,
                                              double[] allosomalTargetBiases, int[] allosomalTargetMasks) {
        Utils.validateArg(autosomalTargetBiases.length > 0, "Autosomal bias array must be non-empty");
        Utils.validateArg(allosomalTargetBiases.length > 0, "Allosomal bias array must be non-empty");
        Utils.validate(autosomalTargetBiases.length == autosomalTargetMasks.length, "Autosomal target bias and mask" +
                " arrays must have the same length");
        Utils.validate(allosomalTargetBiases.length == allosomalTargetMasks.length, "Allosomal target bias and mask" +
                " arrays must have the same length");
        double biasSum = 0;
        double activeTargetCount = 0;
        for (int ti = 0; ti < autosomalTargetBiases.length; ti++) {
            biasSum += autosomalTargetBiases[ti] * autosomalTargetMasks[ti];
            activeTargetCount += autosomalTargetMasks[ti];
        }
        for (int ti = 0; ti < allosomalTargetBiases.length; ti++) {
            biasSum += allosomalTargetBiases[ti] * allosomalTargetMasks[ti];
            activeTargetCount += allosomalTargetMasks[ti];
        }
        if ((int)activeTargetCount == 0) {
            throw new GATKException.ShouldNeverReachHereException("All targets are masked");
        }
        final double meanBias = biasSum / activeTargetCount;
        for (int ti = 0; ti < autosomalTargetBiases.length; ti++) {
            autosomalTargetBiases[ti] /= meanBias;
        }
        for (int ti = 0; ti < allosomalTargetBiases.length; ti++) {
            allosomalTargetBiases[ti] /= meanBias;
        }
    }

    private static void assertSomeUnmaskedTargetsExist(final int[] autosomalTargetMasks, final int[] allosomalTargetMasks) {
        final int numUnmaskedAutosomalTargets = Arrays.stream(autosomalTargetMasks).sum();
        if (numUnmaskedAutosomalTargets == 0) {
            throw new UserException.BadInput("No autosomal targets could be identified. Please ensure (1) the " +
                    "read count collections covers some autosomal targets, and (2) autosomal targets are properly " +
                    "annotated in the provided germline contig ploidy annotations table");
        }

        final int numUnmaskedAllosomalTargets = Arrays.stream(allosomalTargetMasks).sum();
        if (numUnmaskedAllosomalTargets == 0) {
            throw new UserException.BadInput("No allosomal targets could be identified. Please ensure (1) the " +
                    "read count collections covers some allosomal targets, and (2) autosomal targets are properly " +
                    "annotated in the provided germline contig ploidy annotations table");
        }
    }

    private static int[] getPloidyArray(final String sexGenotypeIdentifier,
                                        final List<Target> autosomalTargetList,
                                        final GermlinePloidyAnnotatedTargetCollection ploidyAnnots) {
        return autosomalTargetList.stream()
                .mapToInt(t -> ploidyAnnots.getTargetGermlinePloidyByGenotype(t, sexGenotypeIdentifier)).toArray();
    }

    private static Map<String, int[]> getAllosomalPloidyMap(final List<Target> allosomalTargetList,
                                                            final Set<String> sexGenotypeIdentifiers,
                                                            final GermlinePloidyAnnotatedTargetCollection ploidyAnnots) {
        return Collections.unmodifiableMap(sexGenotypeIdentifiers.stream()
                .map(sexGenotype -> ImmutablePair.of(sexGenotype, getPloidyArray(sexGenotype, allosomalTargetList, ploidyAnnots)))
                .collect(Collectors.toMap(ImmutablePair::getLeft, ImmutablePair::getRight)));
    }

    private static double getAutosomalMeanPloidy(final List<Target> autosomalTargetList, final int[] autosomalTargetPloidies) {
        double autosomalPloidySum = 0;
        for (int ti = 0; ti < autosomalTargetList.size(); ti++) {
            autosomalPloidySum += autosomalTargetPloidies[ti];
        }
        return autosomalPloidySum / autosomalTargetList.size();
    }

    private static void assertBiasAndMaskValidity(final double[] bias, final int[] mask) {
        if (Arrays.stream(mask).anyMatch(entry -> entry != 0 && entry != 1)) {
            throw new GATKException.ShouldNeverReachHereException("The target mask array contains values other than 0 and 1");
        }
        if (Arrays.stream(bias).anyMatch(entry -> entry < 0)) {
            throw new GATKException.ShouldNeverReachHereException("The target bias array contains negative values");
        }
    }

    /**
     * Processes BAIT_COUNT annotations (if present) and updates the bias and mask array accordingly
     * @return true if any corrections applied; false otherwise
     */
    private static boolean processBaitCountAnnotations(final List<Target> targetList, final double[] biases, final int[] masks,
                                                       final String targetClass) {
        /* check if targets have BAIT_COUNT annotations */
        if (targetList.stream().map(Target::getAnnotations)
                .allMatch(annotations -> annotations != null && annotations.hasAnnotation(TargetAnnotation.BAIT_COUNT))) {
            logger.info(String.format("The input target list contains BAIT_COUNT annotations for all %s targets; using" +
                    " bait count as a multiplicative bias factor", targetClass));
            final List<Target> targetsWithSuspiciousBaitAnnotations = targetList.stream()
                    .filter(target -> target.getAnnotations().getDouble(TargetAnnotation.BAIT_COUNT) < MIN_ALLOWED_BAIT_COUNT_PER_TARGET)
                    .collect(Collectors.toList());
            if (!targetsWithSuspiciousBaitAnnotations.isEmpty()) {
                logger.warn(String.format("Some %s targets have suspicious bait counts (less than 1.0) and they were removed;" +
                        " it is highly recommended to reassess the target/bait tables; number of suspicious targets: %d",
                        targetClass, targetsWithSuspiciousBaitAnnotations.size()));
                logger.debug("Suspicious targets: " + targetsWithSuspiciousBaitAnnotations.stream()
                        .map(Target::getName).collect(Collectors.joining(", ")));
            }
            IntStream.range(0, targetList.size())
                    .forEach(ti -> {
                        final double baitCount = targetList.get(ti).getAnnotations().getDouble(TargetAnnotation.BAIT_COUNT);
                        biases[ti] *= baitCount;
                        masks[ti] = (baitCount < MIN_ALLOWED_BAIT_COUNT_PER_TARGET) ? 0 : masks[ti];
                    });
            return true;
        } else {
            logger.warn(String.format("(ignore if the input data is obtained from WGS) some (or all) of the %s targets miss BAIT_COUNT" +
                    " annotation. If the input data is WES, it is highly recommended to provide a target list with" +
                    " BAIT_COUNT annotations. The BAIT_COUNT annotation can be added to a target list using the" +
                    " AnnotateTarget tool and the bait interval list", targetClass));
            return false;
        }
    }

    private static String getTargetCountPerContigStatString(final List<Target> targets, final int[] mask) {
        return IntStream.range(0, targets.size())
                .filter(ti -> mask[ti] == 1)
                .mapToObj(targets::get)
                .collect(Collectors.groupingBy(Target::getContig)).entrySet().stream()
                .map(entry -> String.format("(%s, %s)", entry.getKey(), String.valueOf(entry.getValue().size())))
                .collect(Collectors.joining(", "));
    }

    /**
     * Infer sex genotypes for all samples
     *
     * @return an instance of {@link SexGenotypeDataCollection}
     */
    public SexGenotypeDataCollection inferSexGenotypes() {
        final SexGenotypeDataCollection col = new SexGenotypeDataCollection();
        IntStream.range(0, readCounts.columnNames().size())
                .forEach(si -> {
                    final SexGenotypeData sampleSexGenotype = calculateSexGenotypeData(si);
                    logger.info(sampleSexGenotype);
                    col.add(sampleSexGenotype);
                });
        return col;
    }

    private static void maskTotallyUncoveredTargets(@Nonnull final List<Target> targetList,
                                                    @Nonnull final int[] masks,
                                                    @Nonnull final ReadCountCollection readCounts,
                                                    @Nonnull final String targetClass) {
        if (readCounts.columnNames().size() == 1) {
            logger.warn(String.format("Since only one sample is given for genotyping, the user is responsible for asserting" +
                    " the quality of capture/mapping. In particular, fully uncovered %s targets (e.g. due to bad" +
                    " baits and mapping issues) could not be automatically identified and masked. The presence of" +
                    " such targets may lead to unreliable genotyping results.", targetClass));
            return;
        }
        Utils.validateArg(new HashSet<>(readCounts.targets()).containsAll(targetList), "Some targets are missing in the" +
                " given read count collection");
        Utils.validate(targetList.size() == masks.length, "The target mask and target list must have the same length");
        final int[] uncoveredTargets = IntStream.range(0, targetList.size())
                .filter(ti -> Arrays.stream(readCounts.getRow(targetList.get(ti))).allMatch(count -> (int)count == 0))
                .toArray();
        if (uncoveredTargets.length > 0) {
            logger.info(String.format("Some %s targets were uncovered in the given read count collection for all samples and" +
                    " were excluded (count: %d)", targetClass, uncoveredTargets.length));
            Arrays.stream(uncoveredTargets).forEach(ti -> masks[ti] = 0);
            logger.debug("Excluded targets: " + Arrays.stream(uncoveredTargets).mapToObj(targetList::get)
                .map(Target::getName).collect(Collectors.joining(", ")));
        } else {
            logger.info("No fully uncovered " + targetClass + " targets detected");
        }
    }

    private static void maskTargetsNotOverlappingWithGenotypingIntervals(@Nonnull final List<Target> targetList,
                                                                         @Nonnull final int[] masks,
                                                                         @Nullable final List<SimpleInterval> genotypingIntervals,
                                                                         @Nonnull final String targetClass) {
        if (genotypingIntervals == null) {
            return;
        }
        final int[] outlierTargetIndices = IntStream.range(0, targetList.size())
                .filter(ti -> genotypingIntervals.stream().noneMatch(interval -> interval.overlaps(targetList.get(ti))))
                .toArray();
        if (outlierTargetIndices.length > 0) {
            logger.info(String.format("Some %s targets did not overlap with genotyping intervals and were excluded (count: %d)",
                    targetClass, outlierTargetIndices.length));
            Arrays.stream(outlierTargetIndices).forEach(ti -> masks[ti] = 0);
            logger.debug("Excluded targets: " + Arrays.stream(outlierTargetIndices).mapToObj(targetList::get)
                    .map(Target::getName).collect(Collectors.joining(", ")));
        } else {
            logger.info(String.format("All %s targets overlap with genotyping intervals; keeping all", targetClass));
        }
    }

    /**
     * Estimates read depth per target per homolog for a given sample index in the collection.
     *
     * @param sampleIndex integer index of the sample in the read count collection
     * @return read depth per target per homolog
     */
    private double getSampleReadDepthFromAutosomalTargets(final int sampleIndex) {
        final double[] sampleReadCounts = readCounts.getColumnOnSpecifiedTargets(sampleIndex,
                autosomalTargetList, false);
        final double ploidyError = baselineMappingErrorProbability * meanAutosomalPloidy;
        final double[] readDepthEstimates = IntStream.range(0, sampleReadCounts.length)
                .filter(ti -> autosomalTargetMasks[ti] == 1)
                .mapToDouble(ti -> sampleReadCounts[ti] / ((1.0 - baselineMappingErrorProbability) *
                        autosomalTargetPloidies[ti] * autosomalTargetBiases[ti] + ploidyError))
                .toArray();
        return new Median().evaluate(readDepthEstimates);
    }

    /**
     * Calculates the likelihood of different sex genotypes for a given sample index
     *
     * @param sampleIndex sample index
     * @return an instance of {@link SexGenotypeData}
     */
    private SexGenotypeData calculateSexGenotypeData(final int sampleIndex) {
        final int[] allosomalReadCounts = Arrays.stream(readCounts.getColumnOnSpecifiedTargets(sampleIndex,
                allosomalTargetList, true)).mapToInt(n -> (int)n).toArray();
        final double readDepth = getSampleReadDepthFromAutosomalTargets(sampleIndex);

        logger.debug("Read depth for " + readCounts.columnNames().get(sampleIndex) + ": " + readDepth);
        final List<Double> logLikelihoods = new ArrayList<>();
        final List<String> sexGenotypesList = new ArrayList<>();
        final double ploidyError = baselineMappingErrorProbability * meanAutosomalPloidy;
        for (final String genotypeName : sexGenotypeIdentifiers) {
            /* calculate log likelihood */
            final int[] currentAllosomalTargetPloidies = allosomalTargetPloidies.get(genotypeName);
            final double logLikelihood = IntStream.range(0, allosomalTargetList.size())
                    .filter(ti -> allosomalTargetMasks[ti] == 1)
                    .mapToDouble(ti -> {
                        final double lambda = readDepth * ((1 - baselineMappingErrorProbability) * currentAllosomalTargetPloidies[ti]
                                * allosomalTargetBiases[ti] + ploidyError);
                        return new PoissonDistribution(lambda).logProbability(allosomalReadCounts[ti]);
                    }).sum();
            sexGenotypesList.add(genotypeName);
            logLikelihoods.add(logLikelihood);
        }

        /* infer the most likely sex genotype */
        final Integer[] indices = new Integer[sexGenotypesList.size()];
        IntStream.range(0, sexGenotypesList.size()).forEach(i -> indices[i] = i);
        Arrays.sort(indices, (li, ri) -> Double.compare(logLikelihoods.get(ri), logLikelihoods.get(li)));

        return new SexGenotypeData(readCounts.columnNames().get(sampleIndex),
                sexGenotypesList.get(indices[0]), sexGenotypesList,
                logLikelihoods.stream().mapToDouble(d -> d).toArray());
    }

    @VisibleForTesting
    Set<String> getSexGenotypeIdentifiers() {
        return sexGenotypeIdentifiers;
    }

    @VisibleForTesting
    int[] getAutosomalTargetGermlinePloidies() {
        return autosomalTargetPloidies;
    }

    @VisibleForTesting
    Map<String, int[]> getAllosomalTargetGermlinePloidiesMap() {
        return allosomalTargetPloidies;
    }

    @VisibleForTesting
    List<Target> getAutosomalTargetList() {
        return autosomalTargetList;
    }

    @VisibleForTesting
    List<Target> getAllosomalTargetList() {
        return allosomalTargetList;
    }

    @VisibleForTesting
    int[] getAutosomalTargetMasks() {
        return autosomalTargetMasks;
    }

    @VisibleForTesting
    int[] getAllosomalTargetMasks() {
        return allosomalTargetMasks;
    }

    @VisibleForTesting
    double[] getAutosomalTargetBiases() {
        return  autosomalTargetBiases;
    }

    @VisibleForTesting
    double[] getAllosomalTargetBiases() {
        return allosomalTargetBiases;
    }
}
