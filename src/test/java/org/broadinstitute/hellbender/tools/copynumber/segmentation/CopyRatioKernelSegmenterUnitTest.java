package org.broadinstitute.hellbender.tools.copynumber.segmentation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenterUnitTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioKernelSegmenterUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 1;   //reset seed before each simulated test case

    /**
     * Generates same Gaussian test data as {@link KernelSegmenterUnitTest#dataKernelSegmenter()},
     * but introduces further segments by placing data on different chromosomes.
     */
    @DataProvider(name = "dataCopyRatioKernelSegmenter")
    public Object[][] dataCopyRatioKernelSegmenter() {
        final int numPoints = 1000;

        final Random rng = new Random(RANDOM_SEED);
        rng.setSeed(RANDOM_SEED);
        final List<Double> dataGaussian = IntStream.range(0, numPoints).boxed()
                .map(i -> Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());             //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899

        final List<SimpleInterval> intervals = IntStream.range(0, numPoints).boxed()
                .map(i -> new SimpleInterval(
                        Integer.toString(i / 250 + 1),  //start a new chromosome every 250 points, which adds additional changepoints
                        (i % 250) * 10 + 1,
                        (i % 250) * 10 + 10))         //intervals for copy-ratio data points have length = 10
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

        final CopyRatioSegmentCollection segmentsExpected =
                new CopyRatioSegmentCollection(
                        metadata,
                        Arrays.asList(
                                new CopyRatioSegment(new SimpleInterval("1", 1, 1000), denoisedCopyRatios.getRecords().subList(0, 100)),
                                new CopyRatioSegment(new SimpleInterval("1", 1001, 2000), denoisedCopyRatios.getRecords().subList(100, 200)),
                                new CopyRatioSegment(new SimpleInterval("1", 2001, 2500), denoisedCopyRatios.getRecords().subList(200, 250)),
                                new CopyRatioSegment(new SimpleInterval("2", 1, 500), denoisedCopyRatios.getRecords().subList(250, 300)),
                                new CopyRatioSegment(new SimpleInterval("2", 501, 1500), denoisedCopyRatios.getRecords().subList(300, 400)),
                                new CopyRatioSegment(new SimpleInterval("2", 1501, 2500), denoisedCopyRatios.getRecords().subList(400, 500)),
                                new CopyRatioSegment(new SimpleInterval("3", 1, 1000), denoisedCopyRatios.getRecords().subList(500, 600)),
                                new CopyRatioSegment(new SimpleInterval("3", 1001, 2000), denoisedCopyRatios.getRecords().subList(600, 700)),
                                new CopyRatioSegment(new SimpleInterval("3", 2001, 2500), denoisedCopyRatios.getRecords().subList(700, 750)),
                                new CopyRatioSegment(new SimpleInterval("4", 1, 500), denoisedCopyRatios.getRecords().subList(750, 800)),
                                new CopyRatioSegment(new SimpleInterval("4", 501, 1500), denoisedCopyRatios.getRecords().subList(800, 900)),
                                new CopyRatioSegment(new SimpleInterval("4", 1501, 2500), denoisedCopyRatios.getRecords().subList(900, 1000))));

        return new Object[][]{
                {denoisedCopyRatios, segmentsExpected}
        };
    }

    @Test(dataProvider = "dataCopyRatioKernelSegmenter")
    public void testCopyRatioKernelSegmenter(final CopyRatioCollection denoisedCopyRatios,
                                             final CopyRatioSegmentCollection segmentsExpected) {
        final int maxNumChangepointsPerChromosome = 25;
        final double kernelVariance = 0.;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 2.;
        final double numChangepointsPenaltyLogLinearFactor = 2.;

        final CopyRatioSegmentCollection segments = new CopyRatioKernelSegmenter(denoisedCopyRatios)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVariance, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);
        Assert.assertEquals(segments, segmentsExpected);
    }
}