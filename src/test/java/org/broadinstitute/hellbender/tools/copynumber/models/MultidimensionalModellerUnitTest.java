package org.broadinstitute.hellbender.tools.copynumber.models;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.MultidimensionalSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.MultidimensionalSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tests the MCMC inference performed by {@link MultidimensionalModeller}.  Only recovery of posterior centers is tested.
 * Merging of adjacent similar segments is also tested.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class MultidimensionalModellerUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 13;

    @Test
    public void testMCMC() {
        final int numSegments = 25;
        final int numSamples = 150;
        final int numBurnIn = 50;
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

        //copy-ratio model parameters
        final double varianceCR = 0.01;
        final double outlierProbabilityCR = 0.05;
        final double averageIntervalsPerSegment = 100.;

        //allele-fraction model parameters
        final double meanBiasAF = 1.2;
        final double biasVarianceAF = 0.04;
        final double outlierProbabilityAF = 0.02;
        final AlleleFractionGlobalParameters globalParametersAF = new AlleleFractionGlobalParameters(meanBiasAF, biasVarianceAF, outlierProbabilityAF);
        final double minorAlleleFractionPriorAlpha = 1.;
        final AlleleFractionPrior priorAF = new AlleleFractionPrior(minorAlleleFractionPriorAlpha);
        final double averageHetsPerSegment = 50.;
        final double averageDepthAF = 50.;

        //similar-segment merging parameters
        final int maxNumSmoothingIterations = 10;
        final int numSmoothingIterationsPerFit = 0;
        final double smoothingCredibleIntervalThresholdCopyRatio = 2.;
        final double smoothingCredibleIntervalThresholdAlleleFraction = 2.;

        //recall that both CR and AF data points are at loci 1, 2, 3, etc. and that each segment is on a different contig
        final SampleLocatableMetadata metadata = new SimpleSampleLocatableMetadata(
                "test-sample",
                new SAMSequenceDictionary(IntStream.range(0, numSegments)
                        .mapToObj(i -> new SAMSequenceRecord("chr" + i + 1, 10000))
                        .collect(Collectors.toList())));
        final CopyRatioSimulatedData simulatedDataCR = new CopyRatioSimulatedData(
                metadata, varianceCR, outlierProbabilityCR, numSegments, averageIntervalsPerSegment, rng);
        final AlleleFractionSimulatedData simulatedDataAF = new AlleleFractionSimulatedData(
                metadata, globalParametersAF, numSegments, averageHetsPerSegment, averageDepthAF, rng);

        //we introduce extra segments, which we will later merge to test similar-segment merging
        final MultidimensionalSegmentCollection oversegmentedSegments = new MultidimensionalSegmentCollection(
                metadata,
                constructOversegmentedSegments(simulatedDataCR, simulatedDataAF));

        final MultidimensionalModeller modeller = new MultidimensionalModeller(
                oversegmentedSegments,
                simulatedDataCR.getCopyRatios(),
                simulatedDataAF.getAllelicCounts(), priorAF,
                numSamples, numBurnIn, numSamples, numBurnIn);
        modeller.smoothSegments(maxNumSmoothingIterations, numSmoothingIterationsPerFit, smoothingCredibleIntervalThresholdCopyRatio, smoothingCredibleIntervalThresholdAlleleFraction);

        CopyRatioModellerUnitTest.assertCopyRatioPosteriorCenters(modeller.getCopyRatioModeller(), simulatedDataCR);
        AlleleFractionModellerUnitTest.assertAlleleFractionPosteriorCenters(modeller.getAlleleFractionModeller(), simulatedDataAF);
    }

    private List<MultidimensionalSegment> constructOversegmentedSegments(final CopyRatioSimulatedData simulatedDataCR,
                                                                         final AlleleFractionSimulatedData simulatedDataAF) {
        final int numSegments = simulatedDataCR.getData().getNumSegments();
        final List<String> contigs = simulatedDataCR.getData().getSegments().getRecords().stream()
                .map(SimpleInterval::getContig)
                .distinct()
                .collect(Collectors.toList());
        final List<MultidimensionalSegment> segments = new ArrayList<>(2 * numSegments);    //we split every real segment into two
        for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
            final String contig = contigs.get(segmentIndex);
            final List<CopyRatio> copyRatiosInSegment =
                    simulatedDataCR.getData().getIndexedCopyRatiosInSegment(segmentIndex).stream()
                            .map(icr -> (CopyRatio) icr)
                            .collect(Collectors.toList());
            final List<AllelicCount> allelicCountsInSegment =
                    simulatedDataAF.getData().getIndexedAllelicCountsInSegment(segmentIndex).stream()
                            .map(iac -> (AllelicCount) iac)
                            .collect(Collectors.toList());

            //take half of whichever data source has fewer points, take the same number from the other data source, and make a segment
            final int numPointsCR = copyRatiosInSegment.size();
            final int numPointsAF = allelicCountsInSegment.size();
            final int numPointsMinHalf = Math.min(numPointsCR, numPointsAF) / 2;
            segments.add(new MultidimensionalSegment(
                    new SimpleInterval(contig, 1, numPointsMinHalf),
                    copyRatiosInSegment.subList(0, numPointsMinHalf + 1),
                    allelicCountsInSegment.subList(0, numPointsMinHalf + 1)));
            //add the remaining points to another segment
            final int numPointsMax = Math.max(numPointsCR, numPointsAF);
            segments.add(new MultidimensionalSegment(
                    new SimpleInterval(contig, numPointsMinHalf + 1, numPointsMax),
                    copyRatiosInSegment.subList(numPointsMinHalf + 1, numPointsCR),
                    allelicCountsInSegment.subList(numPointsMinHalf + 1, numPointsAF)));
        }
        return segments;
    }
}