package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.coveragemodel.CoverageModelUtils;
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
     * Autosomal target lengths (number of covered bases)
     */
    private final int[] autosomalTargetLengths;

    /**
     * Map from genotype names to the list of target germline ploidies in the same order as {@code allosomalTargetList}
     */
    private final Map<String, int[]> allosomalTargetPloidies;

    /**
     * Allosomal target lengths
     */
    private final int[] allosomalTargetLengths;

    /**
     * Public constructor.
     *
     * @param rawReadCounts raw read count collection
     * @param contigPloidyAnnotsList list of contig germline ploidy annotations
     * @param baselineMappingErrorProbability typical mapping error probabilty
     */
    public TargetCoverageSexGenotypeCalculator(@Nonnull final ReadCountCollection rawReadCounts,
                                               @Nonnull final List<ContigGermlinePloidyAnnotation> contigPloidyAnnotsList,
                                               final double baselineMappingErrorProbability) {
        this.baselineMappingErrorProbability = ParamUtils.inRange(baselineMappingErrorProbability,
                0, 1, "The baseline mapping error probability must be a number in [0, 1] interval");
        processedReadCounts = processReadCounts(rawReadCounts);

        /* annotate targets with germline ploidy data */
        final GermlinePloidyAnnotatedTargetCollection ploidyAnnots =
                new GermlinePloidyAnnotatedTargetCollection(contigPloidyAnnotsList, processedReadCounts.targets());

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

        /* autosomal targets have the same ploidy for all ploidy classes (this is asserted earlier); so pick the first one */
        final String somePloidyTag = sexGenotypeIdentifiers.iterator().next();
        autosomalTargetPloidies = autosomalTargetList.stream()
                .mapToInt(t -> ploidyAnnots.getTargetGermlinePloidyByGenotype(t, somePloidyTag)).toArray();
        autosomalTargetLengths = autosomalTargetList.stream()
                .mapToInt(t -> t.getEnd() - t.getStart() + 1).toArray();

        /* allosomal targets may have different ploidies for different ploidy classes */
        allosomalTargetPloidies = Collections.unmodifiableMap(sexGenotypeIdentifiers.stream()
                .map(ploidyTag -> ImmutablePair.of(
                        ploidyTag,
                        allosomalTargetList.stream()
                                .mapToInt(t -> ploidyAnnots.getTargetGermlinePloidyByGenotype(t, ploidyTag)).toArray()))
                .collect(Collectors.toMap(ImmutablePair::getLeft, ImmutablePair::getRight)));
        allosomalTargetLengths = allosomalTargetList.stream()
                .mapToInt(t -> t.getEnd() - t.getStart() + 1).toArray();
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
     * Processes raw read counts. If more than one sample is present in the collection,
     * <dl>
     *     <dt> Filters out fully uncovered targets </dt>
     * </dl>
     * @param rawReadCounts raw read count collection
     * @return processed read count collection
     */
    private ReadCountCollection processReadCounts(@Nonnull final ReadCountCollection rawReadCounts) {
        /* remove totally uncovered targets */
        if (rawReadCounts.columnNames().size() > 1) {
            return ReadCountCollectionUtils.removeTotallyUncoveredTargets(rawReadCounts, logger);
        } else {
            final long numUncoveredTargets = rawReadCounts.records().stream()
                    .filter(rec -> (int)rec.getDouble(0) == 0).count();
            final long numAllTargets = rawReadCounts.targets().size();
            logger.info("Since only one sample is given for genotyping, the user is responsible for asserting" +
                    " the aptitude of targets. Fully uncovered (irrelevant) targets can not be automatically" +
                    " identified (total targets: " + numAllTargets + ", uncovered targets: " + numUncoveredTargets + ")");
            return rawReadCounts;
        }
    }

    /**
     * Estimates read depth per base per homolog for a given sample index in the collection.
     *
     * @param sampleIndex integer index of the sample in the read count collection
     * @return read depth per base per homolog
     */
    private double getSampleReadDepthDensityFromAutosomalTargets(final int sampleIndex) {
        final double[] readCounts = processedReadCounts.getColumnOnSpecifiedTargets(sampleIndex, autosomalTargetList, false);
        final int[] mask = IntStream.range(0, readCounts.length).map(i -> 1).toArray();
        /* assume copy ratios and multiplicative biases are 1.0 at this stage */
        final double[] copyRatio = IntStream.range(0, readCounts.length).mapToDouble(i -> 1.0).toArray();
        final double[] multBias = IntStream.range(0, readCounts.length).mapToDouble(i -> 1.0).toArray();

        return CoverageModelUtils.estimateReadDepthDensityFromPoissonModel(readCounts, autosomalTargetLengths,
                multBias, autosomalTargetPloidies, copyRatio, mask);
    }

    /**
     * Calculates the likelihood of different sex genotypes for a given sample index
     *
     * @param sampleIndex sample index
     * @return an instance of {@link SexGenotypeData}
     */
    private SexGenotypeData calculateSexGenotypeData(final int sampleIndex) {
        final int[] allosomalReadCounts = Arrays.stream(processedReadCounts.getColumnOnSpecifiedTargets(sampleIndex, allosomalTargetList,
                false)).mapToInt(n -> (int)n).toArray();
        final double rho = getSampleReadDepthDensityFromAutosomalTargets(sampleIndex);

        final List<Double> logLikelihoods = new ArrayList<>();
        final List<String> sexGenotypesList = new ArrayList<>();
        final int numAllosomalTargets = allosomalTargetLengths.length;

        for (final String genotypeName : sexGenotypeIdentifiers) {
            /* calculate log likelihood */
            final int[] currentAllosomalTargetPloidies = allosomalTargetPloidies.get(genotypeName);
            final double[] poissonParameters = IntStream.range(0, numAllosomalTargets)
                    .mapToDouble(i -> rho * allosomalTargetLengths[i] *
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
