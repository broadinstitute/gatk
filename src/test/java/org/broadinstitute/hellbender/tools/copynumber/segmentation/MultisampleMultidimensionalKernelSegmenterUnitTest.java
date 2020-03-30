package org.broadinstitute.hellbender.tools.copynumber.segmentation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenterUnitTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class MultisampleMultidimensionalKernelSegmenterUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 1;   //reset seed before each simulated test case

    /**
     * Generates same Gaussian and test data as {@link KernelSegmenterUnitTest#dataKernelSegmenter()}
     * and alternate-allele-fraction-like data (similar to zero-mean multimodal test data
     * in {@link KernelSegmenterUnitTest#dataKernelSegmenter()}),
     * but introduces further segments by placing data on different chromosomes.
     * This is just a combination of the test data from
     * {@link CopyRatioKernelSegmenterUnitTest} and {@link AlleleFractionKernelSegmenterUnitTest},
     * but for multiple samples and at lower signal-to-noise ratio.
     */
    @DataProvider(name = "dataMultisampleMultidimensionalKernelSegmenter")
    public Object[][] dataMultisampleMultidimensionalKernelSegmenter() {
        final Random rng = new Random(RANDOM_SEED);

        final int numSamples = 20;
        final List<CopyRatioCollection> denoisedCopyRatiosPerSample = new ArrayList<>(numSamples);
        final List<AllelicCountCollection> allelicCountsPerSample = new ArrayList<>(numSamples);

        final int numIntervals = 1000;
        final int numAllSites = 1000;
        final double noiseLevelCopyRatio = 1.0;
        final double noiseLevelAlleleFraction = 0.25;
        final double homFraction = 0.025;   //low hom fraction minimizes uncertainty in the changepoints coming from runs of adjacent homs near the changepoints

        //generate intervals for copy-ratio data
        final List<SimpleInterval> intervals = IntStream.range(0, numIntervals).boxed()
                .map(i -> new SimpleInterval(
                        Integer.toString(i / 250 + 1),  //start a new chromosome every 250 points, which adds additional changepoints
                        (i % 250) * 10 + 1,
                        (i % 250) * 10 + 10))           //intervals for copy-ratio data points have length = 10
                .collect(Collectors.toList());

        //generate sites for allele-fraction data
        final List<SimpleInterval> allSites = IntStream.range(0, numAllSites).boxed()
                .map(i -> new SimpleInterval(
                        Integer.toString(i / 250 + 1),  //start a new chromosome every 250 points, which adds additional changepoints
                        (i % 250) * 10 + 1,
                        (i % 250) * 10 + 1))            //one site per copy-ratio interval, sites have length = 1
                .collect(Collectors.toList());

        //in all samples, drop half of the sites at random to give some copy-ratio intervals with no allele-fraction sites (to test imputation of allele fraction at 0.5)
        final List<Boolean> isNotDropped = IntStream.range(0, numAllSites).boxed()
                .map(i -> rng.nextBoolean())
                .collect(Collectors.toList());

        final SAMSequenceDictionary sequenceDictionary = new SAMSequenceDictionary(intervals.stream()
                .map(SimpleInterval::getContig)
                .distinct()
                .map(c -> new SAMSequenceRecord(c, 10000))
                .collect(Collectors.toList()));

        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            //generate numIntervals copy ratios
            final List<Double> dataGaussian = IntStream.range(0, numIntervals).boxed()
                    .map(i -> Math.abs(i / 100 - 5) + noiseLevelCopyRatio * rng.nextGaussian())
                    .collect(Collectors.toList());             //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899

            //generate numAllSites alternate-allele fractions
            final List<Double> allMinorAlleleFractions = Arrays.asList(0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25);
            final List<Double> allAlternateAlleleFractions = IntStream.range(0, numAllSites).boxed()
                    .map(i -> rng.nextFloat() < homFraction
                            ? rng.nextBoolean()
                                ? 0. + noiseLevelAlleleFraction * Math.abs(rng.nextGaussian())                                                //hom ref
                                : 1. - noiseLevelAlleleFraction * Math.abs(rng.nextGaussian())                                                //hom alt
                            : rng.nextBoolean()
                                ? Math.max(allMinorAlleleFractions.get(i / 1000) + noiseLevelAlleleFraction * rng.nextGaussian(), 0.)          //het alt minor
                                : Math.min(1. - allMinorAlleleFractions.get(i / 1000) + noiseLevelAlleleFraction * rng.nextGaussian(), 1.))    //het ref minor
                    .map(f -> Math.min(Math.max(f, 0.), 1.))
                    .collect(Collectors.toList());             //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899

            final List<Double> alternateAlleleFractions = IntStream.range(0, numAllSites).boxed()
                    .filter(isNotDropped::get)
                    .map(allAlternateAlleleFractions::get)
                    .collect(Collectors.toList());
            final List<SimpleInterval> sites = IntStream.range(0, numAllSites).boxed()
                    .filter(isNotDropped::get)
                    .map(allSites::get)
                    .collect(Collectors.toList());

            final SampleLocatableMetadata metadata = new SimpleSampleLocatableMetadata(
                    String.format("test-sample-%d", sampleIndex),
                    sequenceDictionary);

            final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(
                    metadata,
                    IntStream.range(0, intervals.size()).boxed()
                            .map(i -> new CopyRatio(intervals.get(i), dataGaussian.get(i)))
                            .collect(Collectors.toList()));
            denoisedCopyRatiosPerSample.add(denoisedCopyRatios);

            final int globalDepth = 100;
            final List<AllelicCount> allelicCountsList = IntStream.range(0, sites.size()).boxed()
                    .map(i -> new AllelicCount(
                            sites.get(i),
                            (int) ((1 - alternateAlleleFractions.get(i)) * globalDepth),
                            (int) (alternateAlleleFractions.get(i) * globalDepth))
            )
                    .collect(Collectors.toList());
            final AllelicCountCollection allelicCounts = new AllelicCountCollection(metadata, allelicCountsList);
            allelicCountsPerSample.add(allelicCounts);
        }

        final SimpleIntervalCollection segmentsExpected =
                new SimpleIntervalCollection(
                        denoisedCopyRatiosPerSample.get(0).getMetadata(),
                        Arrays.asList(
                                new SimpleInterval("1", 1, 1000),
                                new SimpleInterval("1", 1001, 2000),
                                new SimpleInterval("1", 2001, 2500),
                                new SimpleInterval("2", 1, 500),
                                new SimpleInterval("2", 501, 1500),
                                new SimpleInterval("2", 1501, 2500),
                                new SimpleInterval("3", 1, 1000),
                                new SimpleInterval("3", 1001, 2000),
                                new SimpleInterval("3", 2001, 2500),
                                new SimpleInterval("4", 1, 500),
                                new SimpleInterval("4", 501, 1500),
                                new SimpleInterval("4", 1501, 2500)));

        //we do not expect to find all breakpoints in subsets of the generated data (i.e., subsets with fewer samples)
        return new Object[][]{
                {denoisedCopyRatiosPerSample.subList(0, 1), allelicCountsPerSample.subList(0, 1), segmentsExpected, false},
                {denoisedCopyRatiosPerSample.subList(0, 5), allelicCountsPerSample.subList(0, 5), segmentsExpected, false},
                {denoisedCopyRatiosPerSample, allelicCountsPerSample, segmentsExpected, true}
        };
    }

    @Test(dataProvider = "dataMultisampleMultidimensionalKernelSegmenter")
    public void testMultisampleMultidimensionalKernelSegmenter(final List<CopyRatioCollection> denoisedCopyRatiosPerSample,
                                                               final List<AllelicCountCollection> allelicCountsPerSample,
                                                               final SimpleIntervalCollection segmentsExpected,
                                                               final boolean isPassing) {
        final int maxNumChangepointsPerChromosome = 25;
        final double kernelVarianceCopyRatio = 0.;
        final double kernelVarianceAlleleFraction = 0.05;
        final double kernelScalingAlleleFraction = 1.;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 10.;
        final double numChangepointsPenaltyLogLinearFactor = 10.;

        final SimpleIntervalCollection segments = new MultisampleMultidimensionalKernelSegmenter(denoisedCopyRatiosPerSample, allelicCountsPerSample)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelVarianceAlleleFraction,
                        kernelScalingAlleleFraction, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);

        Assert.assertEquals(segments.equals(segmentsExpected), isPassing);
    }
}