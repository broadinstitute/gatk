package org.broadinstitute.hellbender.tools.copynumber.utils.segmentation;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter.ChangepointSortOrder;
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
public final class KernelSegmenterUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 1;   //make sure to reset random seed to this value before each simulated test case

    /**
     * Generates data for a few test cases, including:
     * 1) Gaussian data with changepoints in the mean every 100 points,
     * 2) zero-mean multimodal data, in which the sign of each point is randomly chosen and
     * changepoints occur in the absolute value of the mean every 100 points, and
     * 3) zero-mean Gaussian data with no changepoints.
     */
    @DataProvider(name = "dataKernelSegmenter")
    public Object[][] dataKernelSegmenter() {
        final int numPoints = 1000;

        final BiFunction<Double, Double, Double> linearKernel = (x, y) -> x * y;
        final BiFunction<Double, Double, Double> gaussianKernel = (x, y) -> Math.exp(-(x - y) * (x - y));

        final Random rng = new Random(RANDOM_SEED);
        rng.setSeed(RANDOM_SEED);
        final List<Double> dataGaussian = IntStream.range(0, numPoints)
                .mapToObj(i -> Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());
        final List<Integer> changepointsExpectedGaussian = Arrays.asList(
                299, 699, 99, 899, 399, 199, 599, 499, 799);            //from python implementation

        rng.setSeed(RANDOM_SEED);
        final List<Double> dataZeroMeanMultimodal = IntStream.range(0, numPoints)
                .mapToObj(i -> 2 * (rng.nextBoolean() ? -1 : 1) * Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());
        final List<Integer> changepointsExpectedZeroMeanMultimodal = Arrays.asList(
                499, 599, 799, 699, 399, 899, 299, 199, 99);            //from python implementation

        rng.setSeed(RANDOM_SEED);
        final List<Double> dataGaussianNoSegments = IntStream.range(0, numPoints)
                .mapToObj(i -> 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());
        final List<Integer> changepointsExpectedGaussianNoSegments = Collections.emptyList();

        final List<Double> dataEmpty = Collections.emptyList();
        final List<Integer> changepointsExpectedEmpty = Collections.emptyList();

        return new Object[][]{
                {dataGaussian, linearKernel, changepointsExpectedGaussian},
                {dataZeroMeanMultimodal, gaussianKernel, changepointsExpectedZeroMeanMultimodal},
                {dataGaussianNoSegments, linearKernel, changepointsExpectedGaussianNoSegments},
                {dataEmpty, linearKernel, changepointsExpectedEmpty}
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
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, ChangepointSortOrder.BACKWARD_SELECTION);
        final List<Integer> changepointsIndexSorted = new KernelSegmenter<>(data)
                .findChangepoints(maxNumChangepoints, kernel, kernelApproximationDimension, windowSizes,
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, ChangepointSortOrder.INDEX);

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
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, ChangepointSortOrder.BACKWARD_SELECTION);

        Assert.assertEquals(changepoints, changepointsExpected.isEmpty() ? changepointsExpected : changepointsExpected.subList(0, maxNumChangepoints));
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
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, ChangepointSortOrder.BACKWARD_SELECTION);
        final List<Integer> changepointsIndexSorted = new KernelSegmenter<>(data)
                .findChangepoints(maxNumChangepoints, kernel, kernelApproximationDimension, windowSizes,
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, ChangepointSortOrder.INDEX);

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
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, ChangepointSortOrder.BACKWARD_SELECTION);

        Assert.assertEquals(changepoints.size(), data.isEmpty() ? 0 : maxNumChangepoints);
        Assert.assertEquals(changepoints.subList(0, changepointsExpected.size()), changepointsExpected);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testKernelSegmenterEmptyWindowSizes() {
        final int maxNumChangepoints = 25;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Collections.emptyList();
        final double numChangepointsPenaltyLinearFactor = 0.;
        final double numChangepointsPenaltyLogLinearFactor = 0.;
        final List<Double> data = Arrays.asList(1., 2., 3.);
        final BiFunction<Double, Double, Double> kernel = (x, y) -> x * y;

        new KernelSegmenter<>(data)
                .findChangepoints(maxNumChangepoints, kernel, kernelApproximationDimension, windowSizes,
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, ChangepointSortOrder.BACKWARD_SELECTION);
    }
}