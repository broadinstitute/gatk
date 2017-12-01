package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionGlobalParameters;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by davidben on 10/6/16.
 */
public final class JointAFCRSegmenterUnitTest {

    @Test(enabled=false)
    public void testSegmentation() {
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(563));
        final double hetProportion = 0.25; // probability that a datum is a het i.e. #hets / (#hets + #targets)
        final List<Double> trueMinorAlleleFractions = Arrays.asList(0.12, 0.32, 0.5);
        final List<Double>  trueLog2CopyRatios = Arrays.asList(-2.0, 0.0, 1.7);
        final List<AFCRHiddenState> trueJointStates = IntStream.range(0, trueLog2CopyRatios.size())
                .mapToObj(n -> new AFCRHiddenState(trueMinorAlleleFractions.get(n), trueLog2CopyRatios.get(n)))
                .collect(Collectors.toList());
        final double trueMemoryLength = 1e5;
        final double trueCauchyWidth = 0.2;
        final double trueOutlierProbability = 0.02;
        final int initialNumCRStates = 20;
        final int initialNumAFStates = 20;
        final AlleleFractionGlobalParameters trueAFParams = new AlleleFractionGlobalParameters(1.0, 0.01, 0.01);
        final JointAFCRHMM trueJointModel = new JointAFCRHMM(trueJointStates, trueMemoryLength,
                trueCauchyWidth, trueOutlierProbability);

        // generate joint truth
        final int chainLength = 10000;
        final List<SimpleInterval> positions = CopyRatioSegmenterUnitTest.randomPositions("chr1", chainLength, rng, trueMemoryLength/4);
        final List<Integer> trueHiddenStates = trueJointModel.generateHiddenStateChain(positions);
        final List<AFCRHiddenState> trueAFCRSequence = trueHiddenStates.stream().map(trueJointModel::getHiddenStateValue).collect(Collectors.toList());
        final List<Double> trueLog2CopyRatioSequence = trueAFCRSequence.stream().map(AFCRHiddenState::getLog2CopyRatio).collect(Collectors.toList());
        final List<Double> trueMinorFractionSequence = trueAFCRSequence.stream().map(AFCRHiddenState::getMinorAlleleFraction).collect(Collectors.toList());

        // generate separate af and cr data
        final GammaDistribution biasGenerator = AlleleFractionSegmenterUnitTest.getGammaDistribution(trueAFParams, rng);
        final double outlierProbability = trueAFParams.getOutlierProbability();
        final AllelicCountCollection afData = new AllelicCountCollection();
        final List<Double> crData = new ArrayList<>();
        final List<Target> crTargets = new ArrayList<>();
        for (int n = 0; n < positions.size(); n++) {
            final SimpleInterval position = positions.get(n);
            final AFCRHiddenState jointState = trueAFCRSequence.get(n);
            final double minorFraction = jointState.getMinorAlleleFraction();
            final double log2CopyRatio = jointState.getLog2CopyRatio();

            if (rng.nextDouble() < hetProportion) {  // het datum
                afData.add(AlleleFractionSegmenterUnitTest.generateAllelicCount(minorFraction, position, rng, biasGenerator, outlierProbability));
            } else {    //target datum
                crTargets.add(new Target(position));
                crData.add(CopyRatioSegmenterUnitTest.generateData(trueCauchyWidth, log2CopyRatio, rng));
            }
        }

        final ReadCountCollection rcc = new ReadCountCollection(crTargets, Arrays.asList("SAMPLE"), new Array2DRowRealMatrix(crData.stream().mapToDouble(x->x).toArray()));

        final JointAFCRSegmenter segmenter = JointAFCRSegmenter.createJointSegmenter(initialNumCRStates, rcc, initialNumAFStates, afData);

        final TargetCollection<SimpleInterval> tc = new HashedListTargetCollection<>(positions);
        final List<Pair<SimpleInterval, AFCRHiddenState>> segmentation = segmenter.findSegments();
        final List<ACNVModeledSegment> jointSegments = segmentation.stream()
                .map(pair -> {
                    final SimpleInterval position = pair.getLeft();
                    final AFCRHiddenState jointState = pair.getRight();
                    final PosteriorSummary crSummary = PerformJointSegmentation.errorlessPosterior(jointState.getLog2CopyRatio());
                    final PosteriorSummary afSummary = PerformJointSegmentation.errorlessPosterior(jointState.getMinorAlleleFraction());
                    return new ACNVModeledSegment(position, crSummary, afSummary);
                })
                .collect(Collectors.toList());

        final List<Double> segmentMinorFractions = jointSegments.stream()
                .flatMap(s -> Collections.nCopies(tc.targetCount(s.getInterval()), s.getMinorAlleleFractionPosteriorSummary().getCenter()).stream())
                .collect(Collectors.toList());
        final List<Double> segmentCopyRatios = jointSegments.stream()
                .flatMap(s -> Collections.nCopies(tc.targetCount(s.getInterval()), s.getSegmentMeanPosteriorSummary().getCenter()).stream())
                .collect(Collectors.toList());

        final double averageMinorFractionError = IntStream.range(0, segmentMinorFractions.size())
                .mapToDouble(i -> Math.abs(segmentMinorFractions.get(i) - trueMinorFractionSequence.get(i)))
                .average().getAsDouble();
        final double averageCopyRatioError = IntStream.range(0, trueLog2CopyRatioSequence.size())
                .mapToDouble(i -> Math.abs(segmentCopyRatios.get(i) - trueLog2CopyRatioSequence.get(i)))
                .average().getAsDouble();

        System.out.println(segmentCopyRatios);
        System.out.println(trueLog2CopyRatioSequence);

        //TODO: fix test and use stricter value of delta
        Assert.assertEquals(averageMinorFractionError, 0, 0.1);
        Assert.assertEquals(averageCopyRatioError, 0, 0.1);
    }
}