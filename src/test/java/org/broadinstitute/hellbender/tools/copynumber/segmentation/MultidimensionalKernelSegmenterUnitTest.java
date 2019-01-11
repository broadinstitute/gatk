package org.broadinstitute.hellbender.tools.copynumber.segmentation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.MultidimensionalSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.MultidimensionalSegment;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenterUnitTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class MultidimensionalKernelSegmenterUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 1;   //reset seed before each simulated test case

    /**
     * Generates same Gaussian and test data as {@link KernelSegmenterUnitTest#dataKernelSegmenter()}
     * and alternate-allele-fraction-like data (similar to zero-mean multimodal test data
     * in {@link KernelSegmenterUnitTest#dataKernelSegmenter()}),
     * but introduces further segments by placing data on different chromosomes.
     * This is just a combination of the test data from
     * {@link CopyRatioKernelSegmenterUnitTest} and {@link AlleleFractionKernelSegmenterUnitTest}.
     */
    @DataProvider(name = "dataMultidimensionalKernelSegmenter")
    public Object[][] dataMultidimensionalKernelSegmenter() {
        final Random rng = new Random(RANDOM_SEED);
        rng.setSeed(RANDOM_SEED);

        //generate numIntervals copy ratios
        final int numIntervals = 1000;
        final List<Double> dataGaussian = IntStream.range(0, numIntervals).boxed()
                .map(i -> Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());             //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899

        //generate numAllSites alternate-allele fractions
        final int numAllSites = 1000;
        final double noiseLevel = 0.001;
        final double homFraction = 0.1;     //low hom fraction minimizes uncertainty in the changepoints coming from runs of adjacent homs near the changepoints
        final List<Double> allMinorAlleleFractions = Arrays.asList(0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25);
        final List<Double> allAlternateAlleleFractions = IntStream.range(0, numAllSites).boxed()
                .map(i -> rng.nextFloat() < homFraction
                        ? rng.nextBoolean()
                            ? 0. + noiseLevel * Math.abs(rng.nextGaussian())                                //hom ref
                            : 1. - noiseLevel * Math.abs(rng.nextGaussian())                                //hom alt
                        : rng.nextBoolean()
                            ? Math.max(allMinorAlleleFractions.get(i / 1000) + noiseLevel * rng.nextGaussian(), 0.)          //het alt minor
                            : Math.min(1. - allMinorAlleleFractions.get(i / 1000) + noiseLevel * rng.nextGaussian(), 1.))    //het ref minor
                .collect(Collectors.toList());             //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899

        //generate intervals for copy-ratio data
        final List<SimpleInterval> intervals = IntStream.range(0, numIntervals).boxed()
                .map(i -> new SimpleInterval(
                        Integer.toString(i / 250 + 1),  //start a new chromosome every 250 points, which adds additional changepoints
                        (i % 250) * 10 + 1,
                        (i % 250) * 10 + 10))         //intervals for copy-ratio data points have length = 10
                .collect(Collectors.toList());

        //generate sites for allele-fraction data
        final List<SimpleInterval> allSites = IntStream.range(0, numAllSites).boxed()
                .map(i -> new SimpleInterval(
                        Integer.toString(i / 250 + 1),   //start a new chromosome every 250 points, which adds additional changepoints
                        (i % 250) * 10 + 1,
                        (i % 250) * 10 + 1))           //one site per copy-ratio interval, sites have length = 1
                .collect(Collectors.toList());

        //drop half of the sites at random to give some copy-ratio intervals with no allele-fraction sites (to test imputation of allele fraction at 0.5)
        final List<Boolean> isNotDropped = IntStream.range(0, numAllSites).boxed()
                .map(i -> rng.nextBoolean())
                .collect(Collectors.toList());
        final List<Double> alternateAlleleFractions = IntStream.range(0, numAllSites).boxed()
                .filter(isNotDropped::get)
                .map(allAlternateAlleleFractions::get)
                .collect(Collectors.toList());
        final List<SimpleInterval> sites = IntStream.range(0, numAllSites).boxed()
                .filter(isNotDropped::get)
                .map(allSites::get)
                .collect(Collectors.toList());

        final SampleLocatableMetadata metadata = new SimpleSampleLocatableMetadata(
                "test-sample",
                new SAMSequenceDictionary(intervals.stream()
                        .map(SimpleInterval::getContig)
                        .distinct()
                        .map(c -> new SAMSequenceRecord(c, 10000))
                        .collect(Collectors.toList())));

        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(
                metadata,
                IntStream.range(0, intervals.size()).boxed()
                        .map(i -> new CopyRatio(intervals.get(i), dataGaussian.get(i)))
                        .collect(Collectors.toList()));

        final int globalDepth = 100;
        final List<AllelicCount> allelicCountsList = IntStream.range(0, sites.size()).boxed()
                .map(i -> new AllelicCount(
                        sites.get(i),
                        (int) ((1 - alternateAlleleFractions.get(i)) * globalDepth),
                        (int) (alternateAlleleFractions.get(i) * globalDepth)))
                .collect(Collectors.toList());
        final AllelicCountCollection allelicCounts = new AllelicCountCollection(metadata, allelicCountsList);

        final Comparator<Locatable> comparator = denoisedCopyRatios.getComparator();
        final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector = denoisedCopyRatios.getMidpointOverlapDetector();
        final OverlapDetector<AllelicCount> allelicCountOverlapDetector = allelicCounts.getOverlapDetector();
        final MultidimensionalSegmentCollection segmentsExpected =
                new MultidimensionalSegmentCollection(
                        metadata,
                        Arrays.asList(
                                new MultidimensionalSegment(new SimpleInterval("1", 1, 1000), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector),
                                new MultidimensionalSegment(new SimpleInterval("1", 1001, 2000), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector),
                                new MultidimensionalSegment(new SimpleInterval("1", 2001, 2500), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector),
                                new MultidimensionalSegment(new SimpleInterval("2", 1, 500), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector),
                                new MultidimensionalSegment(new SimpleInterval("2", 501, 1500), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector),
                                new MultidimensionalSegment(new SimpleInterval("2", 1501, 2500), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector),
                                new MultidimensionalSegment(new SimpleInterval("3", 1, 1000), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector),
                                new MultidimensionalSegment(new SimpleInterval("3", 1001, 2000), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector),
                                new MultidimensionalSegment(new SimpleInterval("3", 2001, 2500), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector),
                                new MultidimensionalSegment(new SimpleInterval("4", 1, 500), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector),
                                new MultidimensionalSegment(new SimpleInterval("4", 501, 1500), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector),
                                new MultidimensionalSegment(new SimpleInterval("4", 1501, 2500), comparator, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector)));

        return new Object[][]{
                {denoisedCopyRatios, allelicCounts, segmentsExpected}
        };
    }

    @Test(dataProvider = "dataMultidimensionalKernelSegmenter")
    public void testMultidimensionalKernelSegmenter(final CopyRatioCollection denoisedCopyRatios,
                                                    final AllelicCountCollection allelicCounts,
                                                    final MultidimensionalSegmentCollection segmentsExpected) {
        final int maxNumChangepointsPerChromosome = 25;
        final double kernelVarianceCopyRatio = 0.;
        final double kernelVarianceAlleleFraction = 0.025;
        final double kernelScalingAlleleFraction = 1.;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 2.;
        final double numChangepointsPenaltyLogLinearFactor = 2.;

        final MultidimensionalSegmentCollection segments = new MultidimensionalKernelSegmenter(denoisedCopyRatios, allelicCounts)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelVarianceAlleleFraction,
                        kernelScalingAlleleFraction, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);
        Assert.assertEquals(segments, segmentsExpected);
    }
}