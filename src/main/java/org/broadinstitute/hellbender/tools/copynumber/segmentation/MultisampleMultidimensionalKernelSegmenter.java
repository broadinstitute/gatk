package org.broadinstitute.hellbender.tools.copynumber.segmentation;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Segments copy-ratio data and/or alternate-allele-fraction data from one or more samples using kernel segmentation.
 * Copy-ratio intervals and/or allele-fraction sites must be identical in all samples.  Segments do not span chromosomes.
 * If both types of data are provided, only the first allele-fraction site in each copy-ratio interval is used, and
 * the alternate-allele fraction in copy-ratio intervals that do not contain any sites is imputed to be balanced at 0.5.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class MultisampleMultidimensionalKernelSegmenter {
    private static final Logger logger = LogManager.getLogger(MultisampleMultidimensionalKernelSegmenter.class);

    private enum Mode {
        COPY_RATIO_ONLY, ALLELE_FRACTION_ONLY, COPY_RATIO_AND_ALLELE_FRACTION
    }

    private static final int MIN_NUM_POINTS_REQUIRED_PER_CHROMOSOME = 10;

    //assume alternate-allele fraction is 0.5 for missing data
    private static final SimpleInterval DUMMY_INTERVAL = new SimpleInterval("DUMMY", 1, 1);
    private static final AllelicCount BALANCED_ALLELIC_COUNT = new AllelicCount(DUMMY_INTERVAL, 1, 1);

    //Gaussian kernel for a specified variance; if variance is zero, use a linear kernel
    private static final Function<Double, BiFunction<Double, Double, Double>> KERNEL =
            standardDeviation -> standardDeviation == 0.
                    ? (x, y) -> x * y
                    : (x, y) -> new NormalDistribution(null, x, standardDeviation).density(y);

    private static final class MultidimensionalPoint implements Locatable {
        private final SimpleInterval interval;
        private final double[] log2CopyRatios;
        private final double[] alternateAlleleFractions;

        MultidimensionalPoint(final SimpleInterval interval,
                              final double[] log2CopyRatios,
                              final double[] alternateAlleleFractions) {
            this.interval = interval;
            this.log2CopyRatios = log2CopyRatios;
            this.alternateAlleleFractions = alternateAlleleFractions;
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
    }

    private final Mode mode;
    private final int numSamples;
    private final int numPointsCopyRatio;
    private final int numPointsAlleleFraction;
    private final LocatableMetadata metadata;
    private final Map<String, List<MultidimensionalPoint>> multidimensionalPointsPerChromosome;

    /**
     * @param denoisedCopyRatiosPerSample   non-empty; all copy-ratio intervals identical across samples;
     *                                      number of samples and order identical to that in {@code allelicCountsPerSample};
     *                                      pass list of empty {@link CopyRatioCollection}s for {@code ALLELE_FRACTION_ONLY} mode
     * @param allelicCountsPerSample        non-empty; all allele-fraction sites identical across samples;
     *                                      number of samples and order identical to that in {@code denoisedCopyRatiosPerSample};
     *                                      pass list of empty {@link AllelicCountCollection}s for {@code COPY_RATIO_ONLY} mode
     */
    public MultisampleMultidimensionalKernelSegmenter(final List<CopyRatioCollection> denoisedCopyRatiosPerSample,
                                                      final List<AllelicCountCollection> allelicCountsPerSample) {
        validateInputs(denoisedCopyRatiosPerSample, allelicCountsPerSample);
        numSamples = denoisedCopyRatiosPerSample.size();
        final CopyRatioCollection denoisedCopyRatiosFirstSample = denoisedCopyRatiosPerSample.get(0);
        final AllelicCountCollection allelicCountsFirstSample = allelicCountsPerSample.get(0);
        metadata = denoisedCopyRatiosFirstSample.getMetadata();
        numPointsCopyRatio = denoisedCopyRatiosFirstSample.size();
        numPointsAlleleFraction = allelicCountsFirstSample.size();

        if (numPointsAlleleFraction == 0) {
            mode = Mode.COPY_RATIO_ONLY;
            multidimensionalPointsPerChromosome = IntStream.range(0, numPointsCopyRatio).boxed()
                    .map(i -> new MultidimensionalPoint(
                            denoisedCopyRatiosFirstSample.getRecords().get(i).getInterval(),
                            denoisedCopyRatiosPerSample.stream()
                                    .mapToDouble(cr -> cr.getRecords().get(i).getLog2CopyRatioValue())
                                    .toArray(),
                            null))
                    .collect(Collectors.groupingBy(
                            MultidimensionalPoint::getContig,
                            LinkedHashMap::new,
                            Collectors.toList()));
        } else if (numPointsCopyRatio == 0) {
            mode = Mode.ALLELE_FRACTION_ONLY;
            multidimensionalPointsPerChromosome = IntStream.range(0, numPointsAlleleFraction).boxed()
                    .map(i -> new MultidimensionalPoint(
                            allelicCountsFirstSample.getRecords().get(i).getInterval(),
                            null,
                            allelicCountsPerSample.stream()
                                    .mapToDouble(ac -> ac.getRecords().get(i).getAlternateAlleleFraction())
                                    .toArray()))
                    .collect(Collectors.groupingBy(
                            MultidimensionalPoint::getContig,
                            LinkedHashMap::new,
                            Collectors.toList()));
        } else {
            mode = Mode.COPY_RATIO_AND_ALLELE_FRACTION;
            final OverlapDetector<AllelicCount> allelicCountOverlapDetector = allelicCountsFirstSample.getOverlapDetector();
            final Comparator<Locatable> comparator = denoisedCopyRatiosFirstSample.getComparator();
            final Map<SimpleInterval, Integer> allelicSiteToIndexMap = IntStream.range(0, numPointsAlleleFraction).boxed()
                    .collect(Collectors.toMap(
                            i -> allelicCountsFirstSample.getRecords().get(i).getInterval(),
                            Function.identity(),
                            (u, v) -> {
                                throw new GATKException.ShouldNeverReachHereException("Cannot have duplicate sites.");
                            },   //sites should already be distinct
                            LinkedHashMap::new));
            final Map<Integer, Integer> intervalIndexToSiteIndexMap = IntStream.range(0, numPointsCopyRatio).boxed()
                    .collect(Collectors.toMap(
                            Function.identity(),
                            i -> allelicCountOverlapDetector.getOverlaps(denoisedCopyRatiosFirstSample.getRecords().get(i)).stream()
                                    .map(AllelicCount::getInterval)
                                    .min(comparator::compare)
                                    .map(allelicSiteToIndexMap::get)
                                    .orElse(-1),
                            (u, v) -> {
                                throw new GATKException.ShouldNeverReachHereException("Cannot have duplicate indices.");
                            },
                            LinkedHashMap::new));
            final int numAllelicCountsToUse = (int) intervalIndexToSiteIndexMap.values().stream()
                    .filter(i -> i != -1)
                    .count();
            logger.info(String.format("Using first allelic-count site in each copy-ratio interval (%d / %d) for multidimensional segmentation...",
                    numAllelicCountsToUse, numPointsAlleleFraction));
            multidimensionalPointsPerChromosome = IntStream.range(0, numPointsCopyRatio).boxed()
                    .map(i -> new MultidimensionalPoint(
                            denoisedCopyRatiosFirstSample.getRecords().get(i).getInterval(),
                            denoisedCopyRatiosPerSample.stream()
                                    .mapToDouble(denoisedCopyRatios -> denoisedCopyRatios.getRecords().get(i).getLog2CopyRatioValue())
                                    .toArray(),
                            allelicCountsPerSample.stream()
                                    .map(allelicCounts -> intervalIndexToSiteIndexMap.get(i) != -1
                                            ? allelicCounts.getRecords().get(intervalIndexToSiteIndexMap.get(i))
                                            : BALANCED_ALLELIC_COUNT)
                                    .mapToDouble(AllelicCount::getAlternateAlleleFraction)
                                    .toArray()))
                    .collect(Collectors.groupingBy(
                            MultidimensionalPoint::getContig,
                            LinkedHashMap::new,
                            Collectors.toList()));
        }
    }

    private static void validateInputs(final List<CopyRatioCollection> denoisedCopyRatiosPerSample,
                                       final List<AllelicCountCollection> allelicCountsPerSample) {
        Utils.nonEmpty(denoisedCopyRatiosPerSample);
        Utils.nonEmpty(allelicCountsPerSample);
        Utils.validateArg(denoisedCopyRatiosPerSample.size() == allelicCountsPerSample.size(),
                "Number of copy-ratio and allelic-count collections must be equal.");

        Utils.validateArg(IntStream.range(0, denoisedCopyRatiosPerSample.size())
                        .allMatch(i -> denoisedCopyRatiosPerSample.get(i).getMetadata().equals(allelicCountsPerSample.get(i).getMetadata())),
                "Metadata do not match across copy-ratio and allelic-count collections for the samples.  " +
                        "Check that the sample orders for the corresponding inputs are identical.");

        CopyNumberArgumentValidationUtils.getValidatedSequenceDictionary(
                Stream.of(denoisedCopyRatiosPerSample, allelicCountsPerSample)
                        .flatMap(Collection::stream)
                        .toArray(AbstractLocatableCollection[]::new));

        Utils.validateArg((int) denoisedCopyRatiosPerSample.stream()
                        .map(CopyRatioCollection::getIntervals)
                        .distinct()
                        .count() == 1,
                "Copy-ratio intervals must be identical across all samples.");

        Utils.validateArg((int) allelicCountsPerSample.stream()
                        .map(AllelicCountCollection::getIntervals)
                        .distinct()
                        .count() == 1,
                "Allelic-count sites must be identical across all samples.");
    }

    /**
     * Segments the internally held {@link CopyRatioCollection} and {@link AllelicCountCollection}
     * using a separate {@link KernelSegmenter} for each chromosome.
     * @param kernelVarianceCopyRatio       variance of the Gaussian kernel used for copy-ratio data;
     *                                      if zero, a linear kernel is used instead
     * @param kernelVarianceAlleleFraction  variance of the Gaussian kernel used for allele-fraction data;
     *                                      if zero, a linear kernel is used instead
     * @param kernelScalingAlleleFraction   relative scaling S of the kernel K_AF for allele-fraction data
     *                                      to the kernel K_CR for copy-ratio data;
     *                                      the total kernel is K_CR + S * K_AF
     */
    public SimpleIntervalCollection findSegmentation(final int maxNumSegmentsPerChromosome,
                                                     final double kernelVarianceCopyRatio,
                                                     final double kernelVarianceAlleleFraction,
                                                     final double kernelScalingAlleleFraction,
                                                     final int kernelApproximationDimension,
                                                     final List<Integer> windowSizes,
                                                     final double numChangepointsPenaltyLinearFactor,
                                                     final double numChangepointsPenaltyLogLinearFactor) {
        ParamUtils.isPositive(maxNumSegmentsPerChromosome, "Maximum number of segments must be positive.");
        ParamUtils.isPositiveOrZero(kernelVarianceCopyRatio, "Variance of copy-ratio Gaussian kernel must be non-negative (if zero, a linear kernel will be used).");
        ParamUtils.isPositiveOrZero(kernelVarianceAlleleFraction, "Variance of allele-fraction Gaussian kernel must be non-negative (if zero, a linear kernel will be used).");
        ParamUtils.isPositiveOrZero(kernelScalingAlleleFraction, "Scaling of allele-fraction Gaussian kernel must be non-negative.");
        ParamUtils.isPositive(kernelApproximationDimension, "Dimension of kernel approximation must be positive.");
        Utils.validateArg(windowSizes.stream().allMatch(ws -> ws > 0), "Window sizes must all be positive.");
        Utils.validateArg(new HashSet<>(windowSizes).size() == windowSizes.size(), "Window sizes must all be unique.");
        ParamUtils.isPositiveOrZero(numChangepointsPenaltyLinearFactor,
                "Linear factor for the penalty on the number of changepoints per chromosome must be non-negative.");
        ParamUtils.isPositiveOrZero(numChangepointsPenaltyLogLinearFactor,
                "Log-linear factor for the penalty on the number of changepoints per chromosome must be non-negative.");

        final BiFunction<MultidimensionalPoint, MultidimensionalPoint, Double> kernel = constructKernel(
                kernelVarianceCopyRatio, kernelVarianceAlleleFraction, kernelScalingAlleleFraction);

        final int maxNumChangepointsPerChromosome = maxNumSegmentsPerChromosome - 1;

        logger.info(String.format("Finding changepoints in (%d, %d) data points and %d chromosomes across %d sample(s)...",
                numPointsCopyRatio, numPointsAlleleFraction, multidimensionalPointsPerChromosome.size(), numSamples));

        //loop over chromosomes, find changepoints, and create segments
        final List<SimpleInterval> segments = new ArrayList<>();
        for (final String chromosome : multidimensionalPointsPerChromosome.keySet()) {
            final List<MultidimensionalPoint> multidimensionalPointsInChromosome = multidimensionalPointsPerChromosome.get(chromosome);
            final int numMultidimensionalPointsInChromosome = multidimensionalPointsInChromosome.size();
            logger.info(String.format("Finding changepoints in %d data points in chromosome %s...",
                    numMultidimensionalPointsInChromosome, chromosome));

            if (numMultidimensionalPointsInChromosome < MIN_NUM_POINTS_REQUIRED_PER_CHROMOSOME) {
                logger.warn(String.format("Number of points in chromosome %s (%d) is less than that required (%d), skipping segmentation...",
                        chromosome, numMultidimensionalPointsInChromosome, MIN_NUM_POINTS_REQUIRED_PER_CHROMOSOME));
                final int start = multidimensionalPointsInChromosome.get(0).getStart();
                final int end = multidimensionalPointsInChromosome.get(numMultidimensionalPointsInChromosome - 1).getEnd();
                segments.add(new SimpleInterval(chromosome, start, end));
                continue;
            }

            final List<Integer> changepoints = new ArrayList<>(new KernelSegmenter<>(multidimensionalPointsInChromosome)
                .findChangepoints(maxNumChangepointsPerChromosome, kernel, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, KernelSegmenter.ChangepointSortOrder.INDEX));

            if (!changepoints.contains(numMultidimensionalPointsInChromosome)) {
                changepoints.add(numMultidimensionalPointsInChromosome - 1);
            }
            int previousChangepoint = -1;
            for (final int changepoint : changepoints) {
                final int start = multidimensionalPointsPerChromosome.get(chromosome).get(previousChangepoint + 1).getStart();
                final int end = multidimensionalPointsPerChromosome.get(chromosome).get(changepoint).getEnd();
                segments.add(new SimpleInterval(chromosome, start, end));
                previousChangepoint = changepoint;
            }
        }
        logger.info(String.format("Found %d segments in %d chromosomes across %d sample(s).", segments.size(), multidimensionalPointsPerChromosome.size(), numSamples));
        return new SimpleIntervalCollection(metadata, segments);
    }

    private BiFunction<MultidimensionalPoint, MultidimensionalPoint, Double> constructKernel(final double kernelVarianceCopyRatio,
                                                                                             final double kernelVarianceAlleleFraction,
                                                                                             final double kernelScalingAlleleFraction) {
        final double standardDeviationCopyRatio = Math.sqrt(kernelVarianceCopyRatio);
        final double standardDeviationAlleleFraction = Math.sqrt(kernelVarianceAlleleFraction);
        switch (mode) {
            case COPY_RATIO_ONLY:
                return (p1, p2) -> {
                    double sum = 0.;
                    for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                        sum += KERNEL.apply(standardDeviationCopyRatio).apply(p1.log2CopyRatios[sampleIndex], p2.log2CopyRatios[sampleIndex]);
                    }
                    return sum;
                };
            case ALLELE_FRACTION_ONLY:
                return (p1, p2) -> {
                    double sum = 0.;
                    for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                        sum += KERNEL.apply(standardDeviationAlleleFraction).apply(p1.alternateAlleleFractions[sampleIndex], p2.alternateAlleleFractions[sampleIndex]);
                    }
                    return sum;
                };
            case COPY_RATIO_AND_ALLELE_FRACTION:
                return (p1, p2) -> {
                    double sum = 0.;
                    for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                        sum += KERNEL.apply(standardDeviationCopyRatio).apply(p1.log2CopyRatios[sampleIndex], p2.log2CopyRatios[sampleIndex]) +
                                kernelScalingAlleleFraction * KERNEL.apply(standardDeviationAlleleFraction).apply(p1.alternateAlleleFractions[sampleIndex], p2.alternateAlleleFractions[sampleIndex]);
                    }
                    return sum;
                };
            default:
                throw new GATKException.ShouldNeverReachHereException("Encountered unknown Mode.");
        }
    }
}
