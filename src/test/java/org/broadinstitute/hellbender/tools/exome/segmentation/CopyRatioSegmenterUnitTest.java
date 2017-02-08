package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by davidben on 6/6/16.
 */
public final class CopyRatioSegmenterUnitTest {
    @Test
    public void testSegmentation() {
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(563));

        final List<Double> trueLog2CopyRatios = Arrays.asList(-2.0, 0.0, 1.4);
        final double trueMemoryLength = 1e5;
        final double trueStandardDeviation = 0.2;

        final CopyRatioHMM trueModel = new CopyRatioHMM(trueLog2CopyRatios, trueMemoryLength, trueStandardDeviation);

        final int chainLength = 10000;
        final List<SimpleInterval> positions = randomPositions("chr1", chainLength, rng, trueMemoryLength/4);
        final List<Integer> trueStates = trueModel.generateHiddenStateChain(positions);
        final List<Double> trueLog2CopyRatioSequence = trueStates.stream().map(trueLog2CopyRatios::get).collect(Collectors.toList());

        final List<Double> data = trueLog2CopyRatioSequence.stream()
                .map(cr -> generateData(trueStandardDeviation, cr, rng)).collect(Collectors.toList());

        final List<Target> targets = positions.stream().map(Target::new).collect(Collectors.toList());
        final ReadCountCollection rcc = new ReadCountCollection(targets, Arrays.asList("SAMPLE"), new Array2DRowRealMatrix(data.stream().mapToDouble(x->x).toArray()));
        final CopyRatioSegmenter segmenter = new CopyRatioSegmenter(10, rcc);
        final List<ModeledSegment> segments = segmenter.getModeledSegments();

        final List<Double> segmentCopyRatios = segments.stream()
                .flatMap(s -> Collections.nCopies((int) s.getTargetCount(), s.getSegmentMeanInLog2CRSpace()).stream())
                .collect(Collectors.toList());

        System.out.println(trueLog2CopyRatioSequence);
        System.out.println(segmentCopyRatios);

        final double averageCopyRatioError = IntStream.range(0, trueLog2CopyRatioSequence.size())
                .mapToDouble(i -> Math.abs(segmentCopyRatios.get(i) - trueLog2CopyRatioSequence.get(i)))
                .average().getAsDouble();

        //TODO: fix test and use stricter value of delta
        Assert.assertEquals(averageCopyRatioError, 0, 0.1);
    }

    protected static double generateData(final double trueStandardDeviation, final Double cr, final RandomGenerator rng) {
        return cr + rng.nextGaussian() * trueStandardDeviation;
    }

    @Test
    public void testChromosomesOnDifferentSegments() {
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(563));
        final double[] trueLog2CopyRatios = new double[] {-2.0, 0.0, 1.7};
        final double trueMemoryLength = 1e5;

        final double trueStandardDeviation = 0.2;

        // randomly set positions
        final int chainLength = 100;
        final List<SimpleInterval> positions = randomPositions("chr1", chainLength, rng, trueMemoryLength/4);
        positions.addAll(randomPositions("chr2", chainLength, rng, trueMemoryLength/4));
        positions.addAll(randomPositions("chr3", chainLength, rng, trueMemoryLength/4));

        final int trueState = 2;    //fix everything to the same state 2

        final List<Double> data = new ArrayList<>();
        for (int n = 0; n < positions.size(); n++) {
            final double copyRatio = trueLog2CopyRatios[trueState];
            final double observed = generateData(trueStandardDeviation, copyRatio, rng);
            data.add(observed);
        }

        final List<Target> targets = positions.stream().map(Target::new).collect(Collectors.toList());
        final ReadCountCollection rcc = new ReadCountCollection(targets, Arrays.asList("SAMPLE"), new Array2DRowRealMatrix(data.stream().mapToDouble(x->x).toArray()));
        final CopyRatioSegmenter segmenter = new CopyRatioSegmenter(10, rcc);
        final List<ModeledSegment> segments = segmenter.getModeledSegments();

        //check that each chromosome has at least one segment
        final int numDifferentContigsInSegments = (int) segments.stream().map(ModeledSegment::getContig).distinct().count();
        Assert.assertEquals(numDifferentContigsInSegments, 3);
    }

    public static List<SimpleInterval> randomPositions(final String contig, final int chainLength, final RandomGenerator rng, final double separationScale) {
        final List<SimpleInterval> positions = new ArrayList<>();
        int position = 1;
        for (int n = 0; n < chainLength; n++) {
            position += rng.nextInt((int) separationScale);
            final SimpleInterval interval = new SimpleInterval(contig, position, position);
            positions.add(interval);
        }
        return positions;
    }
}