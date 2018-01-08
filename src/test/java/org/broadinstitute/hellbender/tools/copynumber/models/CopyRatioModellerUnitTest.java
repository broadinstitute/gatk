package org.broadinstitute.hellbender.tools.copynumber.models;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.ParameterDecileCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.ModeledSegment;
import org.broadinstitute.hellbender.utils.mcmc.Decile;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tests the MCMC inference performed by {@link CopyRatioModeller}.  Only recovery of posterior centers is tested.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioModellerUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 13;

    // note: the following tolerance could actually be made much smaller if we used more segments and/or
    // more intervals -- most of the error is the sampling error of a finite simulated data set, not numerical error of MCMC
    private static final double ABSOLUTE_TOLERANCE = 0.01;

    @Test
    public void testMCMC() {
        final double variance = 0.01;
        final double outlierProbability = 0.05;
        final int numSegments = 100;
        final double averageIntervalsPerSegment = 100.;
        final int numSamples = 150;
        final int numBurnIn = 50;
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

        final SampleLocatableMetadata metadata = new SimpleSampleLocatableMetadata(
                "test-sample",
                new SAMSequenceDictionary(IntStream.range(0, numSegments)
                        .mapToObj(i -> new SAMSequenceRecord("chr" + i + 1, 10000))
                        .collect(Collectors.toList())));
        final CopyRatioSimulatedData simulatedData = new CopyRatioSimulatedData(
                metadata, variance, outlierProbability, numSegments, averageIntervalsPerSegment, rng);

        final CopyRatioModeller modeller = new CopyRatioModeller(simulatedData.getData().getCopyRatios(), simulatedData.getData().getSegments());
        modeller.fitMCMC(numSamples, numBurnIn);

        assertCopyRatioPosteriorCenters(modeller, simulatedData);
    }

    static void assertCopyRatioPosteriorCenters(final CopyRatioModeller modeller,
                                                final CopyRatioSimulatedData simulatedData) {
        final CopyRatioState trueState = simulatedData.getTrueState();
        final int numSegments = simulatedData.getData().getNumSegments();

        //check centers from samples
        final List<Double> varianceSamples = modeller.getVarianceSamples();
        final List<Double> outlierProbabilitySamples = modeller.getOutlierProbabilitySamples();
        final List<CopyRatioState.SegmentMeans> segmentMeansSamples = modeller.getSegmentMeansSamples();

        Assert.assertEquals(numSegments, segmentMeansSamples.get(0).size());
        final List<List<Double>> segmentMeansSamplesBySegment = IntStream.range(0, numSegments)
                .mapToObj(i -> segmentMeansSamples.stream().map(s -> s.get(i)).collect(Collectors.toList()))
                .collect(Collectors.toList());

        final double varianceResult = varianceSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final double outlierProbabilityResult = outlierProbabilitySamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final List<Double> segmentMeansResult = segmentMeansSamplesBySegment
                .stream().map(list -> list.stream().mapToDouble(x -> x).average().getAsDouble())
                .collect(Collectors.toList());

        final double totalSegmentError = IntStream.range(0, numSegments)
                .mapToDouble(s -> Math.abs(segmentMeansResult.get(s) - trueState.segmentMean(s)))
                .sum();

        Assert.assertEquals(varianceResult, trueState.variance(), ABSOLUTE_TOLERANCE);
        Assert.assertEquals(outlierProbabilityResult, trueState.outlierProbability(), ABSOLUTE_TOLERANCE);
        Assert.assertEquals(totalSegmentError / numSegments, 0.0, ABSOLUTE_TOLERANCE);

        //check centers from summaries
        final ParameterDecileCollection<CopyRatioParameter> globalParameterDeciles = modeller.getGlobalParameterDeciles();
        final DecileCollection varianceDeciles = globalParameterDeciles.getDeciles(CopyRatioParameter.VARIANCE);
        final double variancePosteriorCenter = varianceDeciles.get(Decile.DECILE_50);
        Assert.assertEquals(variancePosteriorCenter, trueState.variance(), ABSOLUTE_TOLERANCE);

        final DecileCollection outlierProbabilityDeciles = globalParameterDeciles.getDeciles(CopyRatioParameter.OUTLIER_PROBABILITY);
        final double outlierProbabilityPosteriorCenter = outlierProbabilityDeciles.get(Decile.DECILE_50);
        Assert.assertEquals(outlierProbabilityPosteriorCenter, trueState.outlierProbability(), ABSOLUTE_TOLERANCE);

        final List<ModeledSegment.SimplePosteriorSummary> segmentMeansPosteriorSummaries = modeller.getSegmentMeansPosteriorSummaries();
        Assert.assertEquals(numSegments, segmentMeansPosteriorSummaries.size());
        final List<Double> segmentMeansPosteriorCenters = segmentMeansPosteriorSummaries.stream().map(ModeledSegment.SimplePosteriorSummary::getDecile50).collect(Collectors.toList());
        double totalPosteriorCentersSegmentError = 0.0;
        for (int segment = 0; segment < numSegments; segment++) {
            totalPosteriorCentersSegmentError += Math.abs(segmentMeansPosteriorCenters.get(segment) - trueState.segmentMean(segment));
        }
        Assert.assertEquals(totalPosteriorCentersSegmentError / numSegments, 0.0, ABSOLUTE_TOLERANCE);
    }
}