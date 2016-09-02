package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionGlobalParameters;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by davidben on 5/16/16.
 */
public final class AlleleFractionSegmenterUnitTest {
    @Test
    public void testSegmentation() {
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(563));

        final List<Double> trueWeights = Arrays.asList(0.2, 0.5, 0.3);
        final List<Double> trueMinorAlleleFractions = Arrays.asList(0.12, 0.32, 0.5);
        final double trueMemoryLength = 1e5;

        final AlleleFractionGlobalParameters trueParams = new AlleleFractionGlobalParameters(1.0, 0.01, 0.01);

        final AlleleFractionHiddenMarkovModel trueModel = new AlleleFractionHiddenMarkovModel(trueMinorAlleleFractions, trueWeights,
                trueMemoryLength, AllelicPanelOfNormals.EMPTY_PON, trueParams);

        // randomly set positions
        final int chainLength = 10000;
        final List<SimpleInterval> positions = CopyRatioSegmenterUnitTest.randomPositions("chr1", chainLength, rng, trueMemoryLength/4);
        final List<Integer> trueStates = trueModel.generateHiddenStateChain(positions);
        final List<Double> truthMinorFractions = trueStates.stream().map(trueModel::getMinorAlleleFraction).collect(Collectors.toList());
        final AllelicCountCollection counts = generateCounts(truthMinorFractions, positions, rng, trueParams);

        final AlleleFractionSegmenter segmenter = new AlleleFractionSegmenter(10, counts, AllelicPanelOfNormals.EMPTY_PON);
        final List<ModeledSegment> segments = segmenter.getModeledSegments();
        final double[] segmentMinorFractions = segments.stream()
                .flatMap(s -> Collections.nCopies((int) s.getTargetCount(), s.getSegmentMean()).stream())
                .mapToDouble(x->x).toArray();

        final double averageMinorFractionError = IntStream.range(0, truthMinorFractions.size())
                .mapToDouble(n -> Math.abs(segmentMinorFractions[n] - truthMinorFractions.get(n)))
                .average().getAsDouble();

        Assert.assertEquals(averageMinorFractionError, 0, 0.01);
    }

    @Test
    public void testChromosomesOnDifferentSegments() {
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(563));

        final double[] trueMinorAlleleFractions = new double[] {0.12, 0.32, 0.5};
        final double trueMemoryLength = 1e5;
        final AlleleFractionGlobalParameters trueParams = new AlleleFractionGlobalParameters(1.0, 0.01, 0.01);

        // randomly set positions
        final int chainLength = 100;
        final List<SimpleInterval> positions = CopyRatioSegmenterUnitTest.randomPositions("chr1", chainLength, rng, trueMemoryLength/4);
        positions.addAll(CopyRatioSegmenterUnitTest.randomPositions("chr2", chainLength, rng, trueMemoryLength/4));
        positions.addAll(CopyRatioSegmenterUnitTest.randomPositions("chr3", chainLength, rng, trueMemoryLength/4));

        final int trueState = 2;    //fix everything to the same state 2
        final List<Double> minorAlleleFractionSequence = Collections.nCopies(positions.size(), trueMinorAlleleFractions[trueState]);

        final AllelicCountCollection counts = generateCounts(minorAlleleFractionSequence, positions, rng, trueParams);

        final AlleleFractionSegmenter segmenter = new AlleleFractionSegmenter(10, counts, AllelicPanelOfNormals.EMPTY_PON);
        final List<ModeledSegment> segments = segmenter.getModeledSegments();

        //check that each chromosome has at least one segment
        final int numDifferentContigsInSegments = (int) segments.stream().map(ModeledSegment::getContig).distinct().count();
        Assert.assertEquals(numDifferentContigsInSegments, 3);
    }

    //visible for testing joint segmentation
    protected static AllelicCountCollection generateCounts(final List<Double> minorAlleleFractionSequence,
                                                           final List<SimpleInterval> positions,
                                                           final RandomGenerator rng,
                                                           final AlleleFractionGlobalParameters trueParams) {
        //translate to ApacheCommons' parametrization of the gamma distribution
        final GammaDistribution biasGenerator = getGammaDistribution(trueParams, rng);
        final double outlierProbability = trueParams.getOutlierProbability();

        final AllelicCountCollection counts = new AllelicCountCollection();

        for (int n = 0; n < minorAlleleFractionSequence.size(); n++) {
            counts.add(generateAllelicCount(minorAlleleFractionSequence.get(n), positions.get(n), rng, biasGenerator, outlierProbability));
        }
        return counts;
    }

    protected static GammaDistribution getGammaDistribution(final AlleleFractionGlobalParameters trueParams, final RandomGenerator rng) {
        final double gammaShape = trueParams.getMeanBias() * trueParams.getMeanBias() / trueParams.getBiasVariance();
        final double gammaScale = trueParams.getBiasVariance() / trueParams.getMeanBias();
        return new GammaDistribution(rng, gammaShape, gammaScale);
    }

    protected static AllelicCount generateAllelicCount(final double minorFraction, final SimpleInterval position,
                                                       final RandomGenerator rng, final GammaDistribution biasGenerator, final double outlierProbability) {
        final int numReads = 100;
        final double bias = biasGenerator.sample();

        //flip a coin to decide alt minor (alt fraction = minor fraction) or ref minor (alt fraction = 1 - minor fraction)
        final double altFraction =  rng.nextDouble() < 0.5 ? minorFraction : 1 - minorFraction;

        //the probability of an alt read is the alt fraction modified by the bias or, in the case of an outlier, random
        final double pAlt = rng.nextDouble() < outlierProbability ? rng.nextDouble()
                : altFraction / (altFraction + (1 - altFraction) * bias);

        final int numAltReads = new BinomialDistribution(rng, numReads, pAlt).sample();
        final int numRefReads = numReads - numAltReads;
        return new AllelicCount(position, numAltReads, numRefReads);
    }
}