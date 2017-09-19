package org.broadinstitute.hellbender.tools.copynumber.utils.segmentation;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class KernelSegmenterUnitTest extends BaseTest {
    private static final int RANDOM_SEED = 1;   //reset seed before each simulated test case

    @DataProvider(name = "dataKernelSegmenter")
    public Object[][] dataKernelSegmenter() {
        final int numPoints = 1000;

        final BiFunction<Double, Double, Double> linearKernel = (x, y) -> x * y;
        final BiFunction<Double, Double, Double> gaussianKernel = (x, y) -> Math.exp(-(x - y) * (x - y));

        final Random rng = new Random(RANDOM_SEED);
        rng.setSeed(RANDOM_SEED);
        final List<Double> dataGaussian = IntStream.range(0, numPoints).boxed()
                .map(i -> Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());
        final List<Integer> changepointsExpectedGaussian = Arrays.asList(
                299, 699, 99, 899, 399, 199, 599, 499, 799);            //from python implementation

        rng.setSeed(RANDOM_SEED);
        final List<Double> dataZeroMeanMultimodal = IntStream.range(0, numPoints).boxed()
                .map(i -> 2 * (rng.nextBoolean() ? -1 : 1) * Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());
        final List<Integer> changepointsExpectedZeroMeanMultimodal = Arrays.asList(
                499, 599, 799, 699, 399, 899, 299, 199, 99);            //from python implementation

        return new Object[][]{
                {dataGaussian, linearKernel, changepointsExpectedGaussian},
                {dataZeroMeanMultimodal, gaussianKernel, changepointsExpectedZeroMeanMultimodal}
        };
    }

    @Test(dataProvider = "dataKernelSegmenter")
    public void testKernelSegmenter(final List<Double> data,
                                    final BiFunction<Double, Double, Double> kernel,
                                    final List<Integer> changepointsExpected) {
        final int maxNumChangepoints = 25;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 2.;
        final double numChangepointsPenaltyLogLinearFactor = 2.;

        final List<Integer> changepoints = new KernelSegmenter<>(data)
                .findChangepoints(maxNumChangepoints, kernel, kernelApproximationDimension, windowSizes,
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, false);
        final List<Integer> changepointsIndexSorted = new KernelSegmenter<>(data)
                .findChangepoints(maxNumChangepoints, kernel, kernelApproximationDimension, windowSizes,
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, true);

        Assert.assertEquals(changepoints, changepointsExpected);
        Assert.assertEquals(changepointsIndexSorted, changepointsExpected.stream().sorted().collect(Collectors.toList()));
    }

    @Test(dataProvider = "dataKernelSegmenter")
    public void testKernelSegmenterTruncateChangepoints(final List<Double> data,
                                                        final BiFunction<Double, Double, Double> kernel,
                                                        final List<Integer> changepointsExpected) {
        final int maxNumChangepoints = 5;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 2.;
        final double numChangepointsPenaltyLogLinearFactor = 2.;

        final List<Integer> changepoints = new KernelSegmenter<>(data)
                .findChangepoints(maxNumChangepoints, kernel, kernelApproximationDimension, windowSizes,
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, false);

        Assert.assertEquals(changepoints, changepointsExpected.subList(0, maxNumChangepoints));
    }

    @Test(dataProvider = "dataKernelSegmenter")
    public void testKernelSegmenterExtremePenalty(final List<Double> data,
                                    final BiFunction<Double, Double, Double> kernel,
                                    final List<Integer> changepointsExpected) {
        final int maxNumChangepoints = 25;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 10000.;
        final double numChangepointsPenaltyLogLinearFactor = 10000.;

        final List<Integer> changepoints = new KernelSegmenter<>(data)
                .findChangepoints(maxNumChangepoints, kernel, kernelApproximationDimension, windowSizes,
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, false);
        final List<Integer> changepointsIndexSorted = new KernelSegmenter<>(data)
                .findChangepoints(maxNumChangepoints, kernel, kernelApproximationDimension, windowSizes,
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, true);

        Assert.assertEquals(changepoints, Collections.emptyList());
        Assert.assertEquals(changepointsIndexSorted, Collections.emptyList());
    }

    @Test(dataProvider = "dataKernelSegmenter")
    public void testKernelSegmenterNoPenalty(final List<Double> data,
                                             final BiFunction<Double, Double, Double> kernel,
                                             final List<Integer> changepointsExpected) {
        final int maxNumChangepoints = 25;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 0.;
        final double numChangepointsPenaltyLogLinearFactor = 0.;

        final List<Integer> changepoints = new KernelSegmenter<>(data)
                .findChangepoints(maxNumChangepoints, kernel, kernelApproximationDimension, windowSizes,
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, false);

        Assert.assertEquals(changepoints.size(), maxNumChangepoints);
        Assert.assertEquals(changepoints.subList(0, changepointsExpected.size()), changepointsExpected);
    }
}