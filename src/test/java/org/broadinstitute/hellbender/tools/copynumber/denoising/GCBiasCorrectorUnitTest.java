package org.broadinstitute.hellbender.tools.copynumber.denoising;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Random;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class GCBiasCorrectorUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 13;
    private static final int MEAN_READ_DEPTH = 100;
    private static final double NON_GC_BIAS_NOISE_LEVEL = 0.01;

    //a reasonable default GC bias curve
    private static Function<Double, Double> QUADRATIC_GC_BIAS_CURVE = gc -> 0.5 + 2 * gc * (1 - gc);

    private static Pair<RealMatrix, double[]> simulateData(final int numSamples,
                                                           final int numIntervals) {
        final Random random = new Random(RANDOM_SEED);
        final double[] intervalGCContent = IntStream.range(0, numIntervals)
                .mapToDouble(n -> 0.5 + 0.2 * random.nextGaussian())
                .map(x -> Math.min(x, 0.95)).map(x -> Math.max(x, 0.05))
                .toArray();
        final double[] intervalGCBias = Arrays.stream(intervalGCContent)
                .map(QUADRATIC_GC_BIAS_CURVE::apply)
                .toArray();

        //model GC bias along with a small random amount of uniform non-GC-bias noise;
        //remaining noise after GC-bias correction should be only arise from the latter
        final RealMatrix readCounts = new Array2DRowRealMatrix(numSamples, numIntervals);
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int sample, final int interval, final double value) {
                return (int) (MEAN_READ_DEPTH * intervalGCBias[interval] * (1. + NON_GC_BIAS_NOISE_LEVEL * random.nextDouble()));
            }
        });
        return new ImmutablePair<>(readCounts, intervalGCContent);
    }

    @Test
    public void testGCCorrection() {
        final int numSamples = 5;
        final int numIntervals = 10000;

        final Pair<RealMatrix, double[]> data = simulateData(numSamples, numIntervals);
        final RealMatrix readCounts = data.getLeft();
        final double[] intervalGCContent = data.getRight();

        //
        final RealMatrix correctedCoverage = readCounts.copy();
        GCBiasCorrector.correctGCBias(correctedCoverage, intervalGCContent);
        final double[] correctedNoiseBySample = MathUtils.rowStdDevs(correctedCoverage);
        Arrays.stream(correctedNoiseBySample).forEach(x -> Assert.assertTrue(x < NON_GC_BIAS_NOISE_LEVEL * MEAN_READ_DEPTH));

        //check that GC-bias correction is approximately idempotent -- if you correct again, very little should happen
        final RealMatrix recorrectedCoverage = correctedCoverage.copy();
        GCBiasCorrector.correctGCBias(recorrectedCoverage, intervalGCContent);
        final double correctedChange = correctedCoverage.subtract(readCounts).getFrobeniusNorm();
        final double recorrectedChange = recorrectedCoverage.subtract(correctedCoverage).getFrobeniusNorm();
        Assert.assertTrue(recorrectedChange < correctedChange / 10.);
    }
}