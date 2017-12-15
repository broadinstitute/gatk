package org.broadinstitute.hellbender.tools.copynumber.models;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tests the initialization performed by {@link AlleleFractionInitializer}.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionInitializerUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 13;
    private static final double ABSOLUTE_TOLERANCE = 0.01;

    @Test
    public void testInitialization() {
        final double meanBias = 1.2;
        final double biasVariance = 0.04;
        final double outlierProbability = 0.02;
        final AlleleFractionGlobalParameters globalParameters = new AlleleFractionGlobalParameters(meanBias, biasVariance, outlierProbability);
        final int numSegments = 100;
        final double averageHetsPerSegment = 50.;
        final double averageDepth = 50.;
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

        final SampleLocatableMetadata metadata = new SimpleSampleLocatableMetadata(
                "test-sample",
                new SAMSequenceDictionary(IntStream.range(0, numSegments)
                        .mapToObj(i -> new SAMSequenceRecord("chr" + i + 1, 1000))
                        .collect(Collectors.toList())));
        final AlleleFractionSimulatedData simulatedData = new AlleleFractionSimulatedData(
                metadata, globalParameters, numSegments, averageHetsPerSegment, averageDepth, rng);

        final AlleleFractionSegmentedData data = simulatedData.getData();
        final AlleleFractionState initializedState = new AlleleFractionInitializer(data).getInitializedState();

        Assert.assertEquals(initializedState.meanBias(), meanBias, ABSOLUTE_TOLERANCE);
        Assert.assertEquals(initializedState.biasVariance(), biasVariance, ABSOLUTE_TOLERANCE);
        Assert.assertEquals(initializedState.outlierProbability(), outlierProbability, ABSOLUTE_TOLERANCE);

        final double averageMinorFractionError = IntStream.range(0, numSegments)
                .mapToDouble(s -> Math.abs(initializedState.segmentMinorFraction(s) - simulatedData.getTrueState().segmentMinorFraction(s)))
                .average().getAsDouble();
        Assert.assertEquals(averageMinorFractionError, 0, ABSOLUTE_TOLERANCE);
    }
}