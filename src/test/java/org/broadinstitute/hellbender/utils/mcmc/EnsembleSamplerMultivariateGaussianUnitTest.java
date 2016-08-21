package org.broadinstitute.hellbender.utils.mcmc;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedModel.EnsembleBuilder;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.WalkerPosition;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Unit test for affine-invariant ensemble sampling (Goodman & Weare 2010).  Demonstrates application of
 * {@link ModelSampler} to sample an {@link EnsembleBuilder} that is specified using {@link ParameterizedState}
 * and a class that implements {@link DataCollection}.
 * <p>
 *     Test performs sampling of a 2-dimensional Gaussian likelihood with parameters specifying the mean and covariance.
 * </p>
 * <p>
 *     Success of the test is determined by recovery of the input mean and covariance from the samples.
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class EnsembleSamplerMultivariateGaussianUnitTest extends BaseTest {
    private static final double SCALE_PARAMETER = 2.;   //value recommended by Goodman & Weare 2010

    private static final double MEAN_TRUTH_X = 1.;
    private static final double MEAN_TRUTH_Y = 1.;
    private static final double COVARIANCE_TRUTH_XX = 2.;
    private static final double COVARIANCE_TRUTH_XY = 0.5;
    private static final double COVARIANCE_TRUTH_YY = 0.5;
    private static final double[] MEAN_TRUTH = new double[]
            {MEAN_TRUTH_X, MEAN_TRUTH_Y};
    private static final double[][] COVARIANCE_TRUTH = new double[][]{
            {COVARIANCE_TRUTH_XX, COVARIANCE_TRUTH_XY},
            {COVARIANCE_TRUTH_XY, COVARIANCE_TRUTH_YY}};

    private static final double INITIAL_WALKER_BALL_SIZE = 0.1;
    private static final double[] INITIAL_WALKER_BALL_MEAN = new double[]{MEAN_TRUTH_X, MEAN_TRUTH_Y};
    private static final double[][] INITIAL_WALKER_BALL_COVARIANCE = new double[][]{
            {INITIAL_WALKER_BALL_SIZE, 0.},
            {0, INITIAL_WALKER_BALL_SIZE}};

    private static final int NUM_WALKERS = 100;
    private static final int NUM_SAMPLES = 10000;
    private static final int NUM_BURN_IN = 1000;

    private static final double RELATIVE_ERROR_THRESHOLD = 0.02;

    //Create truth distribution to sample
    private static final int RANDOM_SEED = 42;
    private static final RandomGenerator rng =
            RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));
    private static final MultivariateNormalDistribution multivariateNormalDistribution =
            new MultivariateNormalDistribution(rng, MEAN_TRUTH, COVARIANCE_TRUTH);

    //Calculates the log likelihood to be sampled
    private static double logDensity(final ParameterizedState<PointCoordinate> point) {
        return FastMath.log(multivariateNormalDistribution.density(new double[]{point.get(PointCoordinate.X, Double.class), point.get(PointCoordinate.Y, Double.class)}));
    }

    //Calculates relative error between value and truth, with respect to truth; used for checking statistics of
    //posterior samples below.
    private static double relativeError(final double value, final double truth) {
        return Math.abs((value - truth) / truth);
    }

    //Since we don't need any additional quantities to calculate the function to be sampled, we only need a trivial DataCollection
    private final class GaussianDataCollection implements DataCollection {
        public GaussianDataCollection() {}
    }

    //We enumerate the parameters of the state (point coordinates in 2-dimensional space) using an enum that implements the ParameterEnum interface.
    private enum PointCoordinate implements ParameterEnum {
        X, Y
    }

    //We create a ModelSampler helper class to initialize the model state and specify the parameters for the EnsembleBuilder.
    private final class MultivariateGaussianModeller {
        private final ParameterizedModel<PointCoordinate, ParameterizedState<PointCoordinate>, GaussianDataCollection> model;

        private MultivariateGaussianModeller(final Function<ParameterizedState<PointCoordinate>, Double> logDensity) {
            final List<WalkerPosition> initialWalkerPositions =
                    Arrays.stream(new MultivariateNormalDistribution(rng, INITIAL_WALKER_BALL_MEAN, INITIAL_WALKER_BALL_COVARIANCE).sample(NUM_WALKERS))
                    .map(doubleArray -> new WalkerPosition(Doubles.asList(doubleArray))).collect(Collectors.toList());
            final GaussianDataCollection dataCollection = new GaussianDataCollection();
            //Use identity transformation (we take walker space to be equivalent to the parameter space)
            final Function<WalkerPosition, ParameterizedState<PointCoordinate>> transformWalkerPositionToState =
                    wp -> new ParameterizedState<>(Arrays.asList(
                            new Parameter<>(PointCoordinate.X, wp.get(0)),
                            new Parameter<>(PointCoordinate.Y, wp.get(1))));

            //Build the ParameterizedModel using the EnsembleBuilder pattern.
            model = new EnsembleBuilder<>(SCALE_PARAMETER, initialWalkerPositions, dataCollection, transformWalkerPositionToState, logDensity)
                    .build();
        }
    }

    /**
     * Tests sampling of a 2-dimensional Gaussian likelihood via MCMC.  Recovery of input values for the mean and
     * covariance elements from the samples is checked to within a relative error of 1%.
     */
    @Test
    public void testRunMCMCOnMultivariateGaussianModel() {
        //Create new instance of the ModelSampler helper class, passing all quantities needed to initialize state and data.
        final MultivariateGaussianModeller modeller = new MultivariateGaussianModeller(EnsembleSamplerMultivariateGaussianUnitTest::logDensity);
        //Create a ModelSampler, passing the total number of samples (including burn-in samples)
        //and the model held by the ModelSampler.
        final ModelSampler<PointCoordinate, ParameterizedState<PointCoordinate>, GaussianDataCollection> modelSampler =
                new ModelSampler<>(NUM_SAMPLES * NUM_WALKERS, modeller.model);
        //Run the MCMC.
        modelSampler.runMCMC();

        //Get the samples of each of the parameter posteriors (discarding burn-in samples) by passing the
        //parameter name, type, and burn-in number to the getSamples method.
        final double[] xSamples = Doubles.toArray(modelSampler.getSamples(PointCoordinate.X, Double.class, NUM_BURN_IN * NUM_WALKERS));
        final double[] ySamples = Doubles.toArray(modelSampler.getSamples(PointCoordinate.Y, Double.class, NUM_BURN_IN * NUM_WALKERS));

        //Check that the mean and covariance of the samples agree with the input quantities.
        final double xMean = new Mean().evaluate(xSamples);
        Assert.assertEquals(relativeError(xMean, MEAN_TRUTH_X), 0., RELATIVE_ERROR_THRESHOLD);

        final double yMean = new Mean().evaluate(ySamples);
        Assert.assertEquals(relativeError(yMean, MEAN_TRUTH_Y), 0., RELATIVE_ERROR_THRESHOLD);

        final double xxCovariance = new Covariance().covariance(xSamples, xSamples);
        Assert.assertEquals(relativeError(xxCovariance, COVARIANCE_TRUTH_XX), 0., RELATIVE_ERROR_THRESHOLD);

        final double xyCovariance = new Covariance().covariance(xSamples, ySamples);
        Assert.assertEquals(relativeError(xyCovariance, COVARIANCE_TRUTH_XY), 0., RELATIVE_ERROR_THRESHOLD);

        final double yyCovariance = new Covariance().covariance(ySamples, ySamples);
        Assert.assertEquals(relativeError(yyCovariance, COVARIANCE_TRUTH_YY), 0., RELATIVE_ERROR_THRESHOLD);
    }
}