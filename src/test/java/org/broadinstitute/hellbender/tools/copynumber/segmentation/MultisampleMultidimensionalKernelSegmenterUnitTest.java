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

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Marton Kanasz-Nagy &lt;mkanaszn@broadinstitute.org&gt;
 */

public final class MultisampleMultidimensionalKernelSegmenterUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 1;   //reset seed before each simulated test case

    /**
     * Generates a normal and a tumor sample data. The tumor sample is the
     * same Gaussian and test data as {@link KernelSegmenterUnitTest#dataKernelSegmenter()}
     * and alternate-allele-fraction-like data (similar to zero-mean multimodal test data
     * in {@link KernelSegmenterUnitTest#dataKernelSegmenter()}),
     * but introduces further segments by placing data on different chromosomes.
     * This is just a combination of the test data from
     * {@link CopyRatioKernelSegmenterUnitTest} and {@link AlleleFractionKernelSegmenterUnitTest}.
     * The data are identical to the ones in
     * {@link MultidimensionalKernelSegmenterUnitTest#dataMultidimensionalKernelSegmenter()}.
     */

    @DataProvider(name = "dataMultisampleMultidimensionalKernelSegmenterNormalTumor")
    public Object[][] dataMultisampleMultidimensionalKernelSegmenterNormalTumor() {
        final Random rng = new Random(RANDOM_SEED);
        rng.setSeed(RANDOM_SEED);

        //generate numIntervals copy ratios
        final int numIntervals = 1000;
        final List<Double> dataGaussianTumor = IntStream.range(0, numIntervals).boxed()
                .map(i -> Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());             //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899
        final List<Double> dataGaussianNormal = IntStream.range(0, numIntervals).boxed()
                .map(i -> 1 + 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());

        //generate numAllSites alternate-allele fractions
        final int numAllSites = 1000;
        final double noiseLevel = 0.001;
        final double homFraction = 0.1;     //low hom fraction minimizes uncertainty in the changepoints coming from runs of adjacent homs near the changepoints
        final List<Double> allMinorAlleleFractionsTumor = Arrays.asList(0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25);
        final List<Double> allMinorAlleleFractionsNormal = Arrays.asList(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);
        final List<Double> allAlternateAlleleFractionsTumor = IntStream.range(0, numAllSites).boxed()
                .map(i -> rng.nextFloat() < homFraction
                        ? rng.nextBoolean()
                        ? 0. + noiseLevel * Math.abs(rng.nextGaussian())                                //hom ref
                        : 1. - noiseLevel * Math.abs(rng.nextGaussian())                                //hom alt
                        : rng.nextBoolean()
                        ? Math.max(allMinorAlleleFractionsTumor.get(i / 1000) + noiseLevel * rng.nextGaussian(), 0.)          //het alt minor
                        : Math.min(1. - allMinorAlleleFractionsTumor.get(i / 1000) + noiseLevel * rng.nextGaussian(), 1.))    //het ref minor
                .collect(Collectors.toList());             //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899

        final List<Double> allAlternateAlleleFractionsNormal = IntStream.range(0, numAllSites).boxed()
                .map(i -> rng.nextFloat() < homFraction
                        ? rng.nextBoolean()
                        ? 0. + noiseLevel * Math.abs(rng.nextGaussian())                                //hom ref
                        : 1. - noiseLevel * Math.abs(rng.nextGaussian())                                //hom alt
                        : rng.nextBoolean()
                        ? Math.max(allMinorAlleleFractionsNormal.get(i / 1000) + noiseLevel * rng.nextGaussian(), 0.)          //het alt minor
                        : Math.min(1. - allMinorAlleleFractionsNormal.get(i / 1000) + noiseLevel * rng.nextGaussian(), 1.))    //het ref minor
                .collect(Collectors.toList());

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
        final List<Double> alternateAlleleFractionsTumor = IntStream.range(0, numAllSites).boxed()
                .filter(isNotDropped::get)
                .map(allAlternateAlleleFractionsTumor::get)
                .collect(Collectors.toList());
        final List<Double> alternateAlleleFractionsNormal = IntStream.range(0, numAllSites).boxed()
                .filter(isNotDropped::get)
                .map(allAlternateAlleleFractionsNormal::get)
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

        final CopyRatioCollection denoisedCopyRatiosTumor = new CopyRatioCollection(
                metadata,
                IntStream.range(0, intervals.size()).boxed()
                        .map(i -> new CopyRatio(intervals.get(i), dataGaussianTumor.get(i)))
                        .collect(Collectors.toList()));
        final CopyRatioCollection denoisedCopyRatiosNormal = new CopyRatioCollection(
                metadata,
                IntStream.range(0, intervals.size()).boxed()
                        .map(i -> new CopyRatio(intervals.get(i), dataGaussianNormal.get(i)))
                        .collect(Collectors.toList()));

        final int globalDepth = 100;
        final List<AllelicCount> allelicCountsListTumor = IntStream.range(0, sites.size()).boxed()
                .map(i -> new AllelicCount(
                        sites.get(i),
                        (int) ((1 - alternateAlleleFractionsTumor.get(i)) * globalDepth),
                        (int) (alternateAlleleFractionsTumor.get(i) * globalDepth)))
                .collect(Collectors.toList());
        final AllelicCountCollection allelicCountsTumor = new AllelicCountCollection(metadata, allelicCountsListTumor);

        final List<AllelicCount> allelicCountsListNormal = IntStream.range(0, sites.size()).boxed()
                .map(i -> new AllelicCount(
                        sites.get(i),
                        (int) ((1 - alternateAlleleFractionsNormal.get(i)) * globalDepth),
                        (int) (alternateAlleleFractionsNormal.get(i) * globalDepth)))
                .collect(Collectors.toList());
        final AllelicCountCollection allelicCountsNormal = new AllelicCountCollection(metadata, allelicCountsListNormal);


        final Comparator<Locatable> comparatorTumor = denoisedCopyRatiosTumor.getComparator();
        final Comparator<Locatable> comparatorNormal = denoisedCopyRatiosNormal.getComparator();
        final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetectorTumor = denoisedCopyRatiosTumor.getMidpointOverlapDetector();
        final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetectorNormal = denoisedCopyRatiosNormal.getMidpointOverlapDetector();
        final OverlapDetector<AllelicCount> allelicCountOverlapDetectorTumor = allelicCountsTumor.getOverlapDetector();
        final OverlapDetector<AllelicCount> allelicCountOverlapDetectorNormal = allelicCountsNormal.getOverlapDetector();
        final MultidimensionalSegmentCollection segmentsExpectedNormal =
                new MultidimensionalSegmentCollection(
                        metadata,
                        Arrays.asList(
                                new MultidimensionalSegment(new SimpleInterval("1", 1, 2500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("2", 1, 2500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("3", 1, 2500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("4", 1, 2500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal)));
        final MultidimensionalSegmentCollection segmentsExpectedTumor =
                new MultidimensionalSegmentCollection(
                        metadata,
                        Arrays.asList(
                                new MultidimensionalSegment(new SimpleInterval("1", 1, 1000), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor),
                                new MultidimensionalSegment(new SimpleInterval("1", 1001, 2000), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor),
                                new MultidimensionalSegment(new SimpleInterval("1", 2001, 2500), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor),
                                new MultidimensionalSegment(new SimpleInterval("2", 1, 500), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor),
                                new MultidimensionalSegment(new SimpleInterval("2", 501, 1500), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor),
                                new MultidimensionalSegment(new SimpleInterval("2", 1501, 2500), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor),
                                new MultidimensionalSegment(new SimpleInterval("3", 1, 1000), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor),
                                new MultidimensionalSegment(new SimpleInterval("3", 1001, 2000), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor),
                                new MultidimensionalSegment(new SimpleInterval("3", 2001, 2500), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor),
                                new MultidimensionalSegment(new SimpleInterval("4", 1, 500), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor),
                                new MultidimensionalSegment(new SimpleInterval("4", 501, 1500), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor),
                                new MultidimensionalSegment(new SimpleInterval("4", 1501, 2500), comparatorTumor, copyRatioMidpointOverlapDetectorTumor, allelicCountOverlapDetectorTumor)));
        final List<MultidimensionalSegmentCollection> segmentsExpectedNormalTumor = new ArrayList<>();
        segmentsExpectedNormalTumor.add(
                new MultidimensionalSegmentCollection(
                        metadata,
                        Arrays.asList(
                                new MultidimensionalSegment(new SimpleInterval("1", 1, 1000), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("1", 1001, 2000), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("1", 2001, 2500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("2", 1, 500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("2", 501, 1500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("2", 1501, 2500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("3", 1, 1000), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("3", 1001, 2000), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("3", 2001, 2500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("4", 1, 500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("4", 501, 1500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal),
                                new MultidimensionalSegment(new SimpleInterval("4", 1501, 2500), comparatorNormal, copyRatioMidpointOverlapDetectorNormal, allelicCountOverlapDetectorNormal))));
        segmentsExpectedNormalTumor.add(segmentsExpectedTumor);

        return new Object[][]{
                {denoisedCopyRatiosNormal, allelicCountsNormal, segmentsExpectedNormal,
                        denoisedCopyRatiosTumor, allelicCountsTumor, segmentsExpectedTumor,
                        segmentsExpectedNormalTumor}
        };
    }


    /**
     * Generates the same tumor data as
     * {@link MultisampleMultidimensionalKernelSegmenterUnitTest#dataMultisampleMultidimensionalKernelSegmenterNormalTumor()}
     * but in multiple copies that differ only in their noise.
     */
    @DataProvider(name = "dataMultisampleMultidimensionalKernelSegmenterIdenticalSamples")
    public Object[][] dataMultisampleMultidimensionalKernelSegmenterIdenticalSamples() {
        int numSamples = 5;

        final Random rng = new Random(RANDOM_SEED);
        rng.setSeed(RANDOM_SEED);

        //generate numIntervals copy ratios
        final int numIntervals = 1000;
        final List<List<Double>> dataGaussianPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            dataGaussianPerSample.add(IntStream.range(0, numIntervals).boxed()
                    .map(i -> Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian())
                    .collect(Collectors.toList())          //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899
            );
        }

        //generate numAllSites alternate-allele fractions
        final int numAllSites = 1000;
        final double noiseLevel = 0.001;
        final double homFraction = 0.1;     //low hom fraction minimizes uncertainty in the changepoints coming from runs of adjacent homs near the changepoints
        final List<Double> allMinorAlleleFractions = Arrays.asList(0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25);
        final List<List<Double>> allAlternateAlleleFractionsPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            allAlternateAlleleFractionsPerSample.add(IntStream.range(0, numAllSites).boxed()
                    .map(i -> rng.nextFloat() < homFraction
                            ? rng.nextBoolean()
                            ? 0. + noiseLevel * Math.abs(rng.nextGaussian())                                //hom ref
                            : 1. - noiseLevel * Math.abs(rng.nextGaussian())                                //hom alt
                            : rng.nextBoolean()
                            ? Math.max(allMinorAlleleFractions.get(i / 1000) + noiseLevel * rng.nextGaussian(), 0.)          //het alt minor
                            : Math.min(1. - allMinorAlleleFractions.get(i / 1000) + noiseLevel * rng.nextGaussian(), 1.))    //het ref minor
                    .collect(Collectors.toList())             //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899
            );
        }

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
        final List<List<Double>> alternateAlleleFractionsPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            alternateAlleleFractionsPerSample.add(
                    IntStream.range(0, numAllSites).boxed()
                            .filter(isNotDropped::get)
                            .map(allAlternateAlleleFractionsPerSample.get(i_sample)::get)
                            .collect(Collectors.toList())
            );
        }
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

        final List<CopyRatioCollection> denoisedCopyRatiosPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            List<Double> randomGaussianData = dataGaussianPerSample.get(i_sample);
            denoisedCopyRatiosPerSample.add(new CopyRatioCollection(
                    metadata,
                    IntStream.range(0, intervals.size()).boxed()
                            .map(i -> new CopyRatio(intervals.get(i), randomGaussianData.get(i)))
                            .collect(Collectors.toList()))
            );
        }

        final int globalDepth = 100;
        final List<List<AllelicCount>> allelicCountsListPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            List<Double> aafSample = alternateAlleleFractionsPerSample.get(i_sample);
            allelicCountsListPerSample.add(
                    IntStream.range(0, sites.size()).boxed()
                            .map(i -> new AllelicCount(
                                    sites.get(i),
                                    (int) ((1 - aafSample.get(i)) * globalDepth),
                                    (int) (aafSample.get(i) * globalDepth)))
                            .collect(Collectors.toList())
            );
        }
        final List<AllelicCountCollection> allelicCountsPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            allelicCountsPerSample.add(new AllelicCountCollection(metadata, allelicCountsListPerSample.get(i_sample)));
        }

        final List<Comparator<Locatable>> comparatorsPerSample = denoisedCopyRatiosPerSample
                .stream().map(cr -> cr.getComparator()).collect(Collectors.toList());
        final List<OverlapDetector<CopyRatio>> copyRatioMidpointOverlapDetectorsPerSample = denoisedCopyRatiosPerSample
                .stream().map(cr -> cr.getMidpointOverlapDetector()).collect(Collectors.toList());
        final List<OverlapDetector<AllelicCount>> allelicCountOverlapDetectorsPerSample = allelicCountsPerSample
                .stream().map(ac->ac.getOverlapDetector()).collect(Collectors.toList());
        final List<MultidimensionalSegmentCollection> segmentsExpectedPerSample =
                IntStream.range(0, numSamples).boxed()
                        .map(i_sample ->
                                new MultidimensionalSegmentCollection(
                                        metadata,
                                        Arrays.asList(
                                                new MultidimensionalSegment(new SimpleInterval("1", 1, 1000), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("1", 1001, 2000), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("1", 2001, 2500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("2", 1, 500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("2", 501, 1500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("2", 1501, 2500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("3", 1, 1000), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("3", 1001, 2000), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("3", 2001, 2500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("4", 1, 500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("4", 501, 1500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("4", 1501, 2500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample))))
                        ).collect(Collectors.toList());
        return new Object[][]{
                {denoisedCopyRatiosPerSample, allelicCountsPerSample, segmentsExpectedPerSample}
        };
    }

    /**
     * Generates the same tumor data as
     * {@link MultisampleMultidimensionalKernelSegmenterUnitTest#dataMultisampleMultidimensionalKernelSegmenterNormalTumor()}
     * but in multiple copies that differ only in their noise.
     */
    @DataProvider(name = "dataMultisampleMultidimensionalKernelSegmenterTumorMixture")
    public Object[][] dataMultisampleMultidimensionalKernelSegmenterTumorMixture() {
        int numSamples = 3;

        final Random rng = new Random(RANDOM_SEED);
        rng.setSeed(RANDOM_SEED);

        //generate numIntervals copy ratios
        final int numIntervals = 1000;
        final List<List<Double>> dataGaussianPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            double random_cr_level = 5 * rng.nextDouble();
            dataGaussianPerSample.add(IntStream.range(0, numIntervals).boxed()
                    .map(i -> Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian() + random_cr_level)
                    .collect(Collectors.toList())          //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899
            );
        }

        //generate numAllSites alternate-allele fractions
        final int numAllSites = 1000;
        final double noiseLevel = 0.001;
        final double homFraction = 0.1;     //low hom fraction minimizes uncertainty in the changepoints coming from runs of adjacent homs near the changepoints
        final List<List<Double>> allMinorAlleleFractions = new ArrayList<>();
        allMinorAlleleFractions.add(Arrays.asList(0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25));
        allMinorAlleleFractions.add(Arrays.asList(0.45, 0.05, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.45, 0.05, 0.25));
        allMinorAlleleFractions.add(Arrays.asList(0.45, 0.05, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.45, 0.05, 0.25));
        final List<List<Double>> allAlternateAlleleFractionsPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            List<Double> allMinorAlleleFractionsOfSample = allMinorAlleleFractions.get(i_sample);
            allAlternateAlleleFractionsPerSample.add(IntStream.range(0, numAllSites).boxed()
                    .map(i -> rng.nextFloat() < homFraction
                            ? rng.nextBoolean()
                            ? 0. + noiseLevel * Math.abs(rng.nextGaussian())                                //hom ref
                            : 1. - noiseLevel * Math.abs(rng.nextGaussian())                                //hom alt
                            : rng.nextBoolean()
                            ? Math.max(allMinorAlleleFractionsOfSample.get(i / 1000) + noiseLevel * rng.nextGaussian(), 0.)          //het alt minor
                            : Math.min(1. - allMinorAlleleFractionsOfSample.get(i / 1000) + noiseLevel * rng.nextGaussian(), 1.))    //het ref minor
                    .collect(Collectors.toList())             //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899
            );
        }

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

        //drop some fraction of the sites at random for each sample to give some copy-ratio intervals with
        //no allele-fraction sites (to test imputation of allele fraction at 0.5)
        final List<List<Boolean>> isNotDropped = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            isNotDropped.add(
                    IntStream.range(0, numAllSites).boxed()
                            .map(i -> rng.nextDouble() > 0.5/numSamples ? true : false)
                            .collect(Collectors.toList())
            );
        }
        final List<List<Double>> alternateAlleleFractionsPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            alternateAlleleFractionsPerSample.add(
                    IntStream.range(0, numAllSites).boxed()
                            .filter(isNotDropped.get(i_sample)::get)
                            .map(allAlternateAlleleFractionsPerSample.get(i_sample)::get)
                            .collect(Collectors.toList())
            );
        }
        final List<List<SimpleInterval>> sites = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            sites.add(
                    IntStream.range(0, numAllSites).boxed()
                            .filter(isNotDropped.get(i_sample)::get)
                            .map(allSites::get)
                            .collect(Collectors.toList())
            );
        }

        final SampleLocatableMetadata metadata = new SimpleSampleLocatableMetadata(
                "test-sample",
                new SAMSequenceDictionary(intervals.stream()
                        .map(SimpleInterval::getContig)
                        .distinct()
                        .map(c -> new SAMSequenceRecord(c, 10000))
                        .collect(Collectors.toList())));

        final List<CopyRatioCollection> denoisedCopyRatiosPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            List<Double> randomGaussianData = dataGaussianPerSample.get(i_sample);
            denoisedCopyRatiosPerSample.add(new CopyRatioCollection(
                    metadata,
                    IntStream.range(0, intervals.size()).boxed()
                            .map(i -> new CopyRatio(intervals.get(i), randomGaussianData.get(i)))
                            .collect(Collectors.toList()))
            );
        }

        final int globalDepth = 100;
        final List<List<AllelicCount>> allelicCountsListPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            List<Double> aafSample = alternateAlleleFractionsPerSample.get(i_sample);
            List<SimpleInterval> sitesOfSample = sites.get(i_sample);
            allelicCountsListPerSample.add(
                    IntStream.range(0, sitesOfSample.size()).boxed()
                            .map(i -> new AllelicCount(
                                    sitesOfSample.get(i),
                                    (int) ((1 - aafSample.get(i)) * globalDepth),
                                    (int) (aafSample.get(i) * globalDepth)))
                            .collect(Collectors.toList())
            );
        }
        final List<AllelicCountCollection> allelicCountsPerSample = new ArrayList<>();
        for (int i_sample=0; i_sample<numSamples; i_sample++) {
            allelicCountsPerSample.add(new AllelicCountCollection(metadata, allelicCountsListPerSample.get(i_sample)));
        }

        final List<Comparator<Locatable>> comparatorsPerSample = denoisedCopyRatiosPerSample
                .stream().map(cr -> cr.getComparator()).collect(Collectors.toList());
        final List<OverlapDetector<CopyRatio>> copyRatioMidpointOverlapDetectorsPerSample = denoisedCopyRatiosPerSample
                .stream().map(cr -> cr.getMidpointOverlapDetector()).collect(Collectors.toList());
        final List<OverlapDetector<AllelicCount>> allelicCountOverlapDetectorsPerSample = allelicCountsPerSample
                .stream().map(ac->ac.getOverlapDetector()).collect(Collectors.toList());
        final List<MultidimensionalSegmentCollection> segmentsExpectedPerSample =
                IntStream.range(0, numSamples).boxed()
                        .map(i_sample ->
                                new MultidimensionalSegmentCollection(
                                        metadata,
                                        Arrays.asList(
                                                new MultidimensionalSegment(new SimpleInterval("1", 1, 1000), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("1", 1001, 2000), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("1", 2001, 2500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("2", 1, 500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("2", 501, 1500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("2", 1501, 2500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("3", 1, 1000), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("3", 1001, 2000), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("3", 2001, 2500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("4", 1, 500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("4", 501, 1500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample)),
                                                new MultidimensionalSegment(new SimpleInterval("4", 1501, 2500), comparatorsPerSample.get(i_sample), copyRatioMidpointOverlapDetectorsPerSample.get(i_sample), allelicCountOverlapDetectorsPerSample.get(i_sample))))
                        ).collect(Collectors.toList());
        return new Object[][]{
                {denoisedCopyRatiosPerSample, allelicCountsPerSample, segmentsExpectedPerSample}
        };
    }

    /**
     * Tests if multiple sampels but with different levels of noise generate the same segmentation as we would get with
     * a single sample.
     */
    @Test(dataProvider = "dataMultisampleMultidimensionalKernelSegmenterIdenticalSamples")
    public void testMultisampleMultidimensionalKernelSegmenterIdenticalSamples(
            final List<CopyRatioCollection> denoisedCopyRatiosPerSample,
            final List<AllelicCountCollection> allelicCountsPerSample,
            final List<MultidimensionalSegmentCollection> segmentsExpectedPerSample) {

        final int maxNumChangepointsPerChromosome = 25;
        final double kernelVarianceCopyRatio = 0.;
        final double kernelVarianceAlleleFraction = 0.025;
        final double kernelScalingAlleleFraction = 1.;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 2.;
        final double numChangepointsPenaltyLogLinearFactor = 2.;
        final List<CopyRatioCollection> denoisedCopyRatiosSingleSample = new ArrayList<>(Arrays.asList(denoisedCopyRatiosPerSample.get(0)));
        final List<AllelicCountCollection> allelicCountsSingleSample = new ArrayList<>(Arrays.asList(allelicCountsPerSample.get(0)));

        final List<MultidimensionalSegmentCollection> segmentsSingleSample = new MultisampleMultidimensionalKernelSegmenter(
                denoisedCopyRatiosSingleSample,
                allelicCountsSingleSample)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelVarianceAlleleFraction,
                        kernelScalingAlleleFraction, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);
        final List<MultidimensionalSegmentCollection> segmentsPerSample = new MultisampleMultidimensionalKernelSegmenter(
                denoisedCopyRatiosPerSample,
                allelicCountsPerSample)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelVarianceAlleleFraction,
                        kernelScalingAlleleFraction, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);

        Assert.assertEquals(segmentsSingleSample.get(0), segmentsExpectedPerSample.get(0));
        for (int i_sample=0; i_sample<segmentsExpectedPerSample.size(); i_sample++) {
            Assert.assertEquals(segmentsPerSample, segmentsExpectedPerSample);
        }
    }

    @Test(dataProvider = "dataMultisampleMultidimensionalKernelSegmenterNormalTumor")
    public void testMultisampleMultidimensionalKernelSegmenterNormalTumor(
            final CopyRatioCollection denoisedCopyRatiosNormal,
            final AllelicCountCollection allelicCountsNormal,
            final MultidimensionalSegmentCollection segmentsExpectedNormal,
            final CopyRatioCollection denoisedCopyRatiosTumor,
            final AllelicCountCollection allelicCountsTumor,
            final MultidimensionalSegmentCollection segmentsExpectedTumor,
            final List<MultidimensionalSegmentCollection> segmentsExpectedNormalTumor) {

        final int maxNumChangepointsPerChromosome = 25;
        final double kernelVarianceCopyRatio = 0.;
        final double kernelVarianceAlleleFraction = 0.025;
        final double kernelScalingAlleleFraction = 1.;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 2.;
        final double numChangepointsPenaltyLogLinearFactor = 2.;
        final List<CopyRatioCollection> denoisedCopyRatiosPerSampleNormal = new ArrayList<>(Arrays.asList(denoisedCopyRatiosNormal));
        final List<AllelicCountCollection> allelicCountsPerSampleNormal = new ArrayList<>(Arrays.asList(allelicCountsNormal));
        final List<CopyRatioCollection> denoisedCopyRatiosPerSampleTumor = new ArrayList<>(Arrays.asList(denoisedCopyRatiosTumor));
        final List<AllelicCountCollection> allelicCountsPerSampleTumor = new ArrayList<>(Arrays.asList(allelicCountsTumor));
        List<CopyRatioCollection> denoisedCopyRatiosPerSampleNormalTumor = new ArrayList<>();
        denoisedCopyRatiosPerSampleNormalTumor.add(denoisedCopyRatiosNormal);
        denoisedCopyRatiosPerSampleNormalTumor.add(denoisedCopyRatiosTumor);
        List<AllelicCountCollection> allelicCountsPerSampleNormalTumor = new ArrayList<>();
        allelicCountsPerSampleNormalTumor.add(allelicCountsNormal);
        allelicCountsPerSampleNormalTumor.add(allelicCountsTumor);

        final List<MultidimensionalSegmentCollection> segmentsTumor = new MultisampleMultidimensionalKernelSegmenter(
                denoisedCopyRatiosPerSampleTumor,
                allelicCountsPerSampleTumor)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelVarianceAlleleFraction,
                        kernelScalingAlleleFraction, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);

        final List<MultidimensionalSegmentCollection> segmentsNormal = new MultisampleMultidimensionalKernelSegmenter(
                denoisedCopyRatiosPerSampleNormal,
                allelicCountsPerSampleNormal)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelVarianceAlleleFraction,
                        kernelScalingAlleleFraction, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);

        final List<MultidimensionalSegmentCollection> segmentsNormalTumor = new MultisampleMultidimensionalKernelSegmenter(
                denoisedCopyRatiosPerSampleNormalTumor,
                allelicCountsPerSampleNormalTumor)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelVarianceAlleleFraction,
                        kernelScalingAlleleFraction, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);

        Assert.assertEquals(segmentsNormal.get(0), segmentsExpectedNormal);
        Assert.assertEquals(segmentsTumor.get(0), segmentsExpectedTumor);
        Assert.assertEquals(segmentsNormalTumor, segmentsExpectedNormalTumor);
    }

    @Test(dataProvider = "dataMultisampleMultidimensionalKernelSegmenterTumorMixture")
    public void testMultisampleMultidimensionalKernelSegmenterTumorMixture(
            final List<CopyRatioCollection> denoisedCopyRatiosPerSample,
            final List<AllelicCountCollection> allelicCountsPerSample,
            final List<MultidimensionalSegmentCollection> segmentsExpectedPerSample) {

        final int maxNumChangepointsPerChromosome = 25;
        final double kernelVarianceCopyRatio = 0.;
        final double kernelVarianceAlleleFraction = 0.025;
        final double kernelScalingAlleleFraction = 1.;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 2.;
        final double numChangepointsPenaltyLogLinearFactor = 2.;
        final List<CopyRatioCollection> denoisedCopyRatiosSingleSample = new ArrayList<>(Arrays.asList(denoisedCopyRatiosPerSample.get(0)));
        final List<AllelicCountCollection> allelicCountsSingleSample = new ArrayList<>(Arrays.asList(allelicCountsPerSample.get(0)));

        final List<MultidimensionalSegmentCollection> segmentsSingleSample = new MultisampleMultidimensionalKernelSegmenter(
                denoisedCopyRatiosSingleSample,
                allelicCountsSingleSample)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelVarianceAlleleFraction,
                        kernelScalingAlleleFraction, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);
        final List<MultidimensionalSegmentCollection> segmentsPerSample = new MultisampleMultidimensionalKernelSegmenter(
                denoisedCopyRatiosPerSample,
                allelicCountsPerSample)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelVarianceAlleleFraction,
                        kernelScalingAlleleFraction, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);

        Assert.assertEquals(segmentsSingleSample.get(0), segmentsExpectedPerSample.get(0));
        for (int i_sample=0; i_sample<segmentsExpectedPerSample.size(); i_sample++) {
            Assert.assertEquals(segmentsPerSample, segmentsExpectedPerSample);
        }
    }
}

