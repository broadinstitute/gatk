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
 * Tests the MCMC inference performed by {@link AlleleFractionModeller}.  Only recovery of posterior centers is tested.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionModellerUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 13;

    // note: the following tolerance could actually be made much smaller if we used more segments and/or
    // more hets -- most of the error is the sampling error of a finite simulated data set, not numerical error of MCMC
    private static final double ABSOLUTE_TOLERANCE = 0.01;

    @Test
    public void testMCMC() {
        final double meanBias = 1.2;
        final double biasVariance = 0.04;
        final double outlierProbability = 0.02;
        final AlleleFractionGlobalParameters globalParameters = new AlleleFractionGlobalParameters(meanBias, biasVariance, outlierProbability);
        final double minorAlleleFractionPriorAlpha = 1.;
        final AlleleFractionPrior prior = new AlleleFractionPrior(minorAlleleFractionPriorAlpha);
        final int numSegments = 50;
        final double averageHetsPerSegment = 50.;
        final double averageDepth = 50.;
        final int numSamples = 150;
        final int numBurnIn = 50;
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

        final SampleLocatableMetadata metadata = new SimpleSampleLocatableMetadata(
                "test-sample",
                new SAMSequenceDictionary(IntStream.range(0, numSegments)
                        .mapToObj(i -> new SAMSequenceRecord("chr" + i + 1, 10000))
                        .collect(Collectors.toList())));
        final AlleleFractionSimulatedData simulatedData = new AlleleFractionSimulatedData(
                metadata, globalParameters, numSegments, averageHetsPerSegment, averageDepth, rng);

        final AlleleFractionModeller modeller = new AlleleFractionModeller(simulatedData.getData().getAllelicCounts(), simulatedData.getData().getSegments(), prior);
        modeller.fitMCMC(numSamples, numBurnIn);

        assertAlleleFractionPosteriorCenters(modeller, simulatedData);
    }

    static void assertAlleleFractionPosteriorCenters(final AlleleFractionModeller modeller,
                                                     final AlleleFractionSimulatedData simulatedData) {
        final AlleleFractionState trueState = simulatedData.getTrueState();
        final int numSegments = simulatedData.getData().getNumSegments();

        //check centers from samples
        final List<Double> meanBiasSamples = modeller.getMeanBiasSamples();
        final List<Double> biasVarianceSamples = modeller.getBiasVarianceSamples();
        final List<Double> outlierProbabilitySamples = modeller.getOutlierProbabilitySamples();
        final List<AlleleFractionState.MinorFractions> minorFractionsSamples = modeller.getMinorFractionsSamples();

        Assert.assertEquals(numSegments, minorFractionsSamples.get(0).size());
        final List<List<Double>> minorFractionsSamplesBySegment = IntStream.range(0, numSegments)
                .mapToObj(i -> minorFractionsSamples.stream().map(s -> s.get(i)).collect(Collectors.toList()))
                .collect(Collectors.toList());

        final double meanBiasResult = meanBiasSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final double biasVarianceResult = biasVarianceSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final double outlierProbabilityResult = outlierProbabilitySamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final List<Double> minorFractionsResult = minorFractionsSamplesBySegment
                .stream().map(list -> list.stream().mapToDouble(x -> x).average().getAsDouble())
                .collect(Collectors.toList());

        final double totalSegmentError = IntStream.range(0, numSegments)
                .mapToDouble(s -> Math.abs(minorFractionsResult.get(s) - trueState.segmentMinorFraction(s)))
                .sum();

        Assert.assertEquals(meanBiasResult, trueState.meanBias(), ABSOLUTE_TOLERANCE);
        Assert.assertEquals(biasVarianceResult, trueState.biasVariance(), ABSOLUTE_TOLERANCE);
        Assert.assertEquals(outlierProbabilityResult, trueState.outlierProbability(), ABSOLUTE_TOLERANCE);
        Assert.assertEquals(totalSegmentError / numSegments, 0.0, ABSOLUTE_TOLERANCE);

        //check centers from summaries
        final ParameterDecileCollection<AlleleFractionParameter> globalParameterDeciles = modeller.getGlobalParameterDeciles();
        final DecileCollection meanBiasDeciles = globalParameterDeciles.getDeciles(AlleleFractionParameter.MEAN_BIAS);
        final double meanBiasPosteriorCenter = meanBiasDeciles.get(Decile.DECILE_50);
        Assert.assertEquals(meanBiasPosteriorCenter, trueState.meanBias(), ABSOLUTE_TOLERANCE);

        final DecileCollection biasVarianceDeciles = globalParameterDeciles.getDeciles(AlleleFractionParameter.BIAS_VARIANCE);
        final double biasVariancePosteriorCenter = biasVarianceDeciles.get(Decile.DECILE_50);
        Assert.assertEquals(biasVariancePosteriorCenter, trueState.biasVariance(), ABSOLUTE_TOLERANCE);

        final DecileCollection outlierProbabilityDeciles = globalParameterDeciles.getDeciles(AlleleFractionParameter.OUTLIER_PROBABILITY);
        final double outlierProbabilityPosteriorCenter = outlierProbabilityDeciles.get(Decile.DECILE_50);
        Assert.assertEquals(outlierProbabilityPosteriorCenter, trueState.outlierProbability(), ABSOLUTE_TOLERANCE);

        final List<ModeledSegment.SimplePosteriorSummary> minorFractionsPosteriorSummaries = modeller.getMinorAlleleFractionsPosteriorSummaries();
        Assert.assertEquals(numSegments, minorFractionsPosteriorSummaries.size());
        final List<Double> minorFractionsPosteriorCenters = minorFractionsPosteriorSummaries.stream().map(ModeledSegment.SimplePosteriorSummary::getDecile50).collect(Collectors.toList());
        double totalPosteriorCentersSegmentError = 0.0;
        for (int segment = 0; segment < numSegments; segment++) {
            totalPosteriorCentersSegmentError += Math.abs(minorFractionsPosteriorCenters.get(segment) - trueState.segmentMinorFraction(segment));
        }
        Assert.assertEquals(totalPosteriorCentersSegmentError / numSegments, 0.0, ABSOLUTE_TOLERANCE);
    }
}