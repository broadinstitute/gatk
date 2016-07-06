package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * Created by davidben on 6/6/16.
 */
public final class CopyRatioSegmenterUnitTest {
    @Test
    public void testSegmentation() {
        final double[] trueWeights = new double[] {0.2, 0.5, 0.3};
        final double[] trueLog2CopyRatios = new double[] {-2.0, 0.0, 1.7};
        final double trueMemoryLength = 1e5;

        final double trueStandardDeviation = 0.2;

        final CopyRatioHiddenMarkovModel trueModel = new CopyRatioHiddenMarkovModel(trueLog2CopyRatios, trueWeights,
                trueMemoryLength, trueStandardDeviation);

        // randomly set positions
        final Random random = new Random(271);
        final int chainLength = 10000;
        final List<SimpleInterval> positions = randomPositions("chr1", chainLength, random, trueMemoryLength/4);
        final List<Integer> trueStates = trueModel.generateHiddenStateChain(positions);

        final List<Double> data = new ArrayList<>();
        for (int n = 0; n < positions.size(); n++) {
            final double copyRatio = trueLog2CopyRatios[trueStates.get(n)];
            final double observed = copyRatio + random.nextGaussian() * trueStandardDeviation;
            data.add(observed);
        }

        final CopyRatioSegmenter segmenter = new CopyRatioSegmenter(10, positions, data);
        final List<ModeledSegment> segments = segmenter.findSegments();
        final double[] segmentCopyRatios = segments.stream()
                .flatMap(s -> Collections.nCopies((int) s.getTargetCount(), s.getSegmentMean()).stream())
                .mapToDouble(x -> x).toArray();
        final double[] truthCopyRatios = trueStates.stream().mapToDouble(trueModel::getCopyRatio).toArray();

        final double averageCopyRatioError = IntStream.range(0, truthCopyRatios.length)
                .mapToDouble(n -> Math.abs(segmentCopyRatios[n] - truthCopyRatios[n]))
                .average().getAsDouble();

        Assert.assertEquals(averageCopyRatioError, 0, 0.025);
    }

    @Test
    public void testChromosomesOnDifferentSegments() {
        final double[] trueLog2CopyRatios = new double[] {-2.0, 0.0, 1.7};
        final double trueMemoryLength = 1e5;

        final double trueStandardDeviation = 0.2;

        // randomly set positions
        final Random random = new Random(271);
        final int chainLength = 100;
        final List<SimpleInterval> positions = randomPositions("chr1", chainLength, random, trueMemoryLength/4);
        positions.addAll(randomPositions("chr2", chainLength, random, trueMemoryLength/4));
        positions.addAll(randomPositions("chr3", chainLength, random, trueMemoryLength/4));

        final int trueState = 2;    //fix everything to the same state 2

        final List<Double> data = new ArrayList<>();
        for (int n = 0; n < positions.size(); n++) {
            final double copyRatio = trueLog2CopyRatios[trueState];
            final double observed = copyRatio + random.nextGaussian() * trueStandardDeviation;
            data.add(observed);
        }

        final CopyRatioSegmenter segmenter = new CopyRatioSegmenter(10, positions, data);
        final List<ModeledSegment> segments = segmenter.findSegments();

        //check that each chromosome has at least one segment
        final int numDifferentContigsInSegments = (int) segments.stream().map(ModeledSegment::getContig).distinct().count();
        Assert.assertEquals(numDifferentContigsInSegments, 3);
    }

    public static List<SimpleInterval> randomPositions(final String contig, final int chainLength, final Random random, final double separationScale) {
        final List<SimpleInterval> positions = new ArrayList<>();
        int position = 1;
        for (int n = 0; n < chainLength; n++) {
            position += random.nextInt((int) separationScale);
            final SimpleInterval interval = new SimpleInterval(contig, position, position);
            positions.add(interval);
        }
        return positions;
    }
}