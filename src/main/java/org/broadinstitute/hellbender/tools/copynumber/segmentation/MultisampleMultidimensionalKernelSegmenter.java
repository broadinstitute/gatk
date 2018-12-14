package org.broadinstitute.hellbender.tools.copynumber.segmentation;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.MultidimensionalSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.MultidimensionalSegment;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Segments copy-ratio and alternate-allele-fraction data using kernel segmentation.  Segments do not span chromosomes.
 * Only the first allele-fraction site in each copy-ratio interval is used.  The alternate-allele fraction in
 * copy-ratio intervals that do not contain any sites is imputed to be balanced at 0.5.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Marton Kanasz-Nagy &lt;mkanaszn@broadinstitute.org&gt;
 */

public final class MultisampleMultidimensionalKernelSegmenter {
    private static final Logger logger = LogManager.getLogger(MultisampleMultidimensionalKernelSegmenter.class);

    private static final int MIN_NUM_POINTS_REQUIRED_PER_CHROMOSOME = 10;

    //assume alternate-allele fraction is 0.5 for missing data
    private static final SimpleInterval DUMMY_INTERVAL = new SimpleInterval("DUMMY", 1, 1);
    private static final AllelicCount BALANCED_ALLELIC_COUNT = new AllelicCount(DUMMY_INTERVAL, 1, 1);

    //Gaussian kernel for a specified variance; if variance is zero, use a linear kernel
    private static final Function<Double, BiFunction<Double, Double, Double>> KERNEL =
            standardDeviation -> standardDeviation == 0.
                    ? (x, y) -> x * y
                    : (x, y) -> new NormalDistribution(null, x, standardDeviation).density(y);

    private static final class MultisampleMultidimensionalPoint implements Locatable {
        private final SimpleInterval interval;
        private final List<Double> log2CopyRatiosPerSample;
        private final List<Double> alternateAlleleFractionsPerSample;
        private int nSamples;

        MultisampleMultidimensionalPoint(final SimpleInterval interval,
                                         List<Double> log2CopyRatiosPerSample,
                                         List<Double> alternateAlleleFractionsPerSample) {
            Utils.nonNull(log2CopyRatiosPerSample);
            Utils.nonNull(alternateAlleleFractionsPerSample);
            Utils.validateArg(log2CopyRatiosPerSample.size() == alternateAlleleFractionsPerSample.size(),
                    "Number of copy ratio and allele fraciont samples needs to be the same.");
            this.interval = interval;
            this.log2CopyRatiosPerSample = log2CopyRatiosPerSample;
            this.alternateAlleleFractionsPerSample = alternateAlleleFractionsPerSample;
            this.nSamples = log2CopyRatiosPerSample.size();
        }

        public final void add(SimpleInterval interval, double log2CopyRatio, double alternateAlleleFraction) {
            Utils.validateArg(this.interval == interval, "Intervals need to match.");
            this.log2CopyRatiosPerSample.add(log2CopyRatio);
            this.alternateAlleleFractionsPerSample.add(alternateAlleleFraction);
            this.nSamples += 1;
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

        public int getNumberOfSamples() {
            return this.nSamples;
        }

        public List<Double> getLog2CopyRatiosPerSample() {
            return this.log2CopyRatiosPerSample;
        }

        public List<Double> getAlternateAlleleFractionsPerSample() {
            return this.alternateAlleleFractionsPerSample;
        }
    }

    private final List<CopyRatioCollection> denoisedCopyRatiosPerSample;
    private final List<OverlapDetector<CopyRatio>> copyRatioMidpointOverlapDetectorPerSample;
    private final List<AllelicCountCollection> allelicCountsPerSample;
    private final List<OverlapDetector<AllelicCount>> allelicCountOverlapDetectorPerSample;
    private final List<Comparator<Locatable>> comparatorsPerSample;
    private final Map<String, List<MultisampleMultidimensionalPoint>> multisampleMultidimensionalPointsPerChromosome;
    private final int numberOfSamples;


    public MultisampleMultidimensionalKernelSegmenter(final List<CopyRatioCollection> denoisedCopyRatiosPerSample,
                                                      final List<AllelicCountCollection> allelicCountsPerSample) {
        Utils.nonNull(denoisedCopyRatiosPerSample);
        Utils.nonNull(allelicCountsPerSample);
        Utils.validateArg(denoisedCopyRatiosPerSample.size() == allelicCountsPerSample.size()
                        && denoisedCopyRatiosPerSample.size() > 0,
                "Number of copy ratio and allele fraciont samples needs to be the same and non-zero.");
        for (int i=0; i<denoisedCopyRatiosPerSample.size(); i++) {
            Utils.nonNull(denoisedCopyRatiosPerSample.get(i));
            Utils.nonNull(allelicCountsPerSample.get(i));
            Utils.validateArg(denoisedCopyRatiosPerSample.get(i).getMetadata().equals(allelicCountsPerSample.get(i).getMetadata()),
                    "Metadata do not match.");
        }
        this.numberOfSamples = denoisedCopyRatiosPerSample.size();
        this.denoisedCopyRatiosPerSample = denoisedCopyRatiosPerSample;
        this.copyRatioMidpointOverlapDetectorPerSample = denoisedCopyRatiosPerSample
                .stream()
                .map(cr -> cr.getMidpointOverlapDetector()).collect(Collectors.toList());
        this.allelicCountsPerSample = allelicCountsPerSample;
        this.allelicCountOverlapDetectorPerSample = allelicCountsPerSample
                .stream()
                .map(ac -> ac.getOverlapDetector()).collect(Collectors.toList());
        final List<Integer> numAllelicCountsToUsePerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<this.numberOfSamples; i_sample++) {
            numAllelicCountsToUsePerSample.add((int) denoisedCopyRatiosPerSample.get(i_sample).getRecords().stream()
                    .filter(allelicCountOverlapDetectorPerSample.get(i_sample)::overlapsAny)
                    .count());
            logger.info(String.format("Using first allelic-count site in each copy-ratio interval (%d / %d) for multidimensional segmentation...",
                    numAllelicCountsToUsePerSample.get(i_sample), allelicCountsPerSample.get(i_sample).size()));
        }
        this.comparatorsPerSample = denoisedCopyRatiosPerSample.stream().map(cr -> cr.getComparator()).collect(Collectors.toList());
        List<MultisampleMultidimensionalPoint> multisampleMultidimensionalPointList = denoisedCopyRatiosPerSample
                .get(0)
                .getRecords()
                .stream()
                .map(cr -> new MultisampleMultidimensionalKernelSegmenter.MultisampleMultidimensionalPoint(
                        cr.getInterval(),
                        new ArrayList<Double>(),
                        new ArrayList<Double>())).collect(Collectors.toList());

        for (int i_sample=0; i_sample<denoisedCopyRatiosPerSample.size(); i_sample++) {
            for (int i_segment=0; i_segment<denoisedCopyRatiosPerSample.get(i_sample).getRecords().size(); i_segment++) {
                CopyRatio cr = denoisedCopyRatiosPerSample.get(i_sample).getRecords().get(i_segment);
                multisampleMultidimensionalPointList.get(i_segment).add(
                        cr.getInterval(),
                        cr.getLog2CopyRatioValue(),
                        allelicCountOverlapDetectorPerSample.get(i_sample).getOverlaps(cr).stream()
                                .min(comparatorsPerSample.get(i_sample)::compare)
                                .orElse(BALANCED_ALLELIC_COUNT).getAlternateAlleleFraction()
                );
            }
        }

        multisampleMultidimensionalPointsPerChromosome = multisampleMultidimensionalPointList.stream()
                .collect(Collectors.groupingBy(
                        MultisampleMultidimensionalKernelSegmenter.MultisampleMultidimensionalPoint::getContig,
                        LinkedHashMap::new,
                        Collectors.toList()));
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
    public List<MultidimensionalSegmentCollection> findSegmentation(final int maxNumChangepointsPerChromosome,
                                                                    final double kernelVarianceCopyRatio,
                                                                    final double kernelVarianceAlleleFraction,
                                                                    final double kernelScalingAlleleFraction,
                                                                    final int kernelApproximationDimension,
                                                                    final List<Integer> windowSizes,
                                                                    final double numChangepointsPenaltyLinearFactor,
                                                                    final double numChangepointsPenaltyLogLinearFactor) {
        ParamUtils.isPositiveOrZero(maxNumChangepointsPerChromosome, "Maximum number of changepoints must be non-negative.");
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

        final BiFunction<MultisampleMultidimensionalPoint, MultisampleMultidimensionalPoint, Double> kernel = constructKernel(
                kernelVarianceCopyRatio, kernelVarianceAlleleFraction, kernelScalingAlleleFraction);

        for (int i_sample=0; i_sample<this.numberOfSamples; i_sample++) {
            logger.info(String.format("Finding changepoints in (%d, %d) data points and %d chromosomes...",
                    denoisedCopyRatiosPerSample.get(i_sample).size(),
                    allelicCountsPerSample.get(i_sample).size(),
                    multisampleMultidimensionalPointsPerChromosome.size()));
        }

        //loop over chromosomes, find changepoints, and create allele-fraction segments
        final List<ArrayList<MultidimensionalSegment>> segmentsPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<this.numberOfSamples; i_sample++) {
            segmentsPerSample.add(new ArrayList<>());
        }
        for (final String chromosome : this.multisampleMultidimensionalPointsPerChromosome.keySet()) {
            final List<MultisampleMultidimensionalPoint> multisampleMultidimensionalPointsInChromosome = multisampleMultidimensionalPointsPerChromosome.get(chromosome);
            final int numMultisampleMultidimensionalPointsInChromosome = multisampleMultidimensionalPointsInChromosome.size();
            logger.info(String.format("Finding changepoints in %d data points in chromosome %s...",
                    numMultisampleMultidimensionalPointsInChromosome, chromosome));

            if (numMultisampleMultidimensionalPointsInChromosome < MIN_NUM_POINTS_REQUIRED_PER_CHROMOSOME) {
                logger.warn(String.format("Number of points in chromosome %s (%d) is less than that required (%d), skipping segmentation...",
                        chromosome, numMultisampleMultidimensionalPointsInChromosome, MIN_NUM_POINTS_REQUIRED_PER_CHROMOSOME));
                final int start = multisampleMultidimensionalPointsInChromosome.get(0).getStart();
                final int end = multisampleMultidimensionalPointsInChromosome.get(numMultisampleMultidimensionalPointsInChromosome - 1).getEnd();
                for (int i_sample=0; i_sample<this.numberOfSamples; i_sample++) {
                    segmentsPerSample.get(i_sample).add(new MultidimensionalSegment(
                            new SimpleInterval(chromosome, start, end),
                            this.comparatorsPerSample.get(i_sample),
                            this.copyRatioMidpointOverlapDetectorPerSample.get(i_sample),
                            this.allelicCountOverlapDetectorPerSample.get(i_sample)));
                }
                continue;
            }

            final List<Integer> changepoints = new ArrayList<>(new KernelSegmenter<>(multisampleMultidimensionalPointsInChromosome)
                    .findChangepoints(maxNumChangepointsPerChromosome, kernel, kernelApproximationDimension,
                            windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, KernelSegmenter.ChangepointSortOrder.INDEX));

            if (!changepoints.contains(numMultisampleMultidimensionalPointsInChromosome)) {
                changepoints.add(numMultisampleMultidimensionalPointsInChromosome - 1);
            }
            int previousChangepoint = -1;
            for (final int changepoint : changepoints) {
                final int start = this.multisampleMultidimensionalPointsPerChromosome.get(chromosome).get(previousChangepoint + 1).getStart();
                final int end = this.multisampleMultidimensionalPointsPerChromosome.get(chromosome).get(changepoint).getEnd();
                for (int i_sample=0; i_sample<this.numberOfSamples; i_sample++) {
                    segmentsPerSample.get(i_sample).add(new MultidimensionalSegment(
                            new SimpleInterval(chromosome, start, end),
                            this.comparatorsPerSample.get(i_sample),
                            this.copyRatioMidpointOverlapDetectorPerSample.get(i_sample),
                            this.allelicCountOverlapDetectorPerSample.get(i_sample)));
                }
                previousChangepoint = changepoint;
            }
        }
        logger.info(String.format("Found %d segments in %d chromosomes.", segmentsPerSample.get(0).size(), multisampleMultidimensionalPointsPerChromosome.keySet().size()));

        List<MultidimensionalSegmentCollection> segmentationsPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<this.numberOfSamples; i_sample++) {
            segmentationsPerSample.add(new MultidimensionalSegmentCollection(
                    allelicCountsPerSample.get(i_sample).getMetadata(),
                    segmentsPerSample.get(i_sample)
            ));
        }
        return segmentationsPerSample;
    }

    private BiFunction<MultisampleMultidimensionalPoint, MultisampleMultidimensionalPoint, Double> constructKernel(final double kernelVarianceCopyRatio,
                                                                                                                   final double kernelVarianceAlleleFraction,
                                                                                                                   final double kernelScalingAlleleFraction) {
        final double standardDeviationCopyRatio = Math.sqrt(kernelVarianceCopyRatio);
        final double standardDeviationAlleleFraction = Math.sqrt(kernelVarianceAlleleFraction);
        return (p1, p2) -> {
            double kernelSum = 0.;
            for (int i_sample=0; i_sample<this.numberOfSamples; i_sample++) {
                kernelSum += KERNEL.apply(standardDeviationCopyRatio).apply(p1.log2CopyRatiosPerSample.get(i_sample), p2.log2CopyRatiosPerSample.get(i_sample)) +
                        kernelScalingAlleleFraction * KERNEL.apply(standardDeviationAlleleFraction).apply(p1.alternateAlleleFractionsPerSample.get(i_sample), p2.alternateAlleleFractionsPerSample.get(i_sample));
            }
            return kernelSum;
        };
    }
}
