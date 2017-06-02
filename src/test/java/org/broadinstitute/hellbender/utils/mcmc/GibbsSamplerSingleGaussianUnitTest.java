package org.broadinstitute.hellbender.utils.mcmc;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.function.Function;

/**
 * Unit test for {@link GibbsSampler}.  Demonstrates application of {@link GibbsSampler} to a {@link ParameterizedModel}
 * that is specified using {@link ParameterizedState} and a class that implements {@link DataCollection}.
 * <p>
 *     Test performs Bayesian inference of a Gaussian model with 2 global parameters specifying the variance and the mean.
 * </p>
 * <p>
 *     Data consists of a list of 10000 datapoints drawn from a normal distribution with unity variance and mean.
 * </p>
 * <p>
 *     Success of the test is determined by recovery of the input variance and mean,
 *     as well as agreement of the standard deviations of the parameter posteriors with those given by both the
 *     python package emcee (see http://dan.iel.fm/emcee for details) and numerical evaluation in Mathematica
 *     of the analytic forms of the posteriors.
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class GibbsSamplerSingleGaussianUnitTest extends BaseTest {
    private static final int NUM_DATAPOINTS = 10000;

    private static final double VARIANCE_MIN = 0.;
    private static final double VARIANCE_MAX = Double.POSITIVE_INFINITY;
    private static final double VARIANCE_WIDTH = 0.1;
    private static final double VARIANCE_INITIAL = 5.;
    private static final double VARIANCE_TRUTH = 1.;
    private static final double VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.014;

    private static final double MEAN_WIDTH = 0.1;
    private static final double MEAN_INITIAL = 5.;
    private static final double MEAN_TRUTH = 1.;
    private static final double MEAN_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.01;

    private static final int NUM_SAMPLES = 500;
    private static final int NUM_BURN_IN = 250;

    //test specifications
    private static final double RELATIVE_ERROR_THRESHOLD_FOR_CENTERS = 0.01;
    private static final double RELATIVE_ERROR_THRESHOLD_FOR_STANDARD_DEVIATIONS = 0.1;

    //Create dataset of 10000 datapoints drawn from a normal distribution Normal(MEAN_TRUTH, VARIANCE_TRUTH)
    private static final int RANDOM_SEED = 42;
    private static final RandomGenerator rng =
            RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));
    private static final NormalDistribution normalDistribution =
            new NormalDistribution(rng, MEAN_TRUTH, Math.sqrt(VARIANCE_TRUTH));
    private static final List<Double> datapointsList =
            Doubles.asList(normalDistribution.sample(NUM_DATAPOINTS));

    //Calculates the exponent for a normal distribution; used in log-likelihood calculation below.
    private static double normalTerm(final double quantity, final double mean, final double variance) {
        return (quantity - mean) * (quantity - mean) / (2 * variance);
    }

    //Calculates relative error between x and xTrue, with respect to xTrue; used for checking statistics of
    //posterior samples below.
    private static double relativeError(final double x, final double xTrue) {
        return Math.abs((x - xTrue) / xTrue);
    }

    //The datapoints used by the samplers must be contained in a class that implements DataCollection.
    private final class GaussianDataCollection implements DataCollection {
        private final List<Double> datapoints;

        public GaussianDataCollection(final List<Double> datapoints) {
            this.datapoints = new ArrayList<>(datapoints);
        }

        public List<Double> getDatapoints() {
            return Collections.unmodifiableList(datapoints);
        }
    }

    //We enumerate the parameters of the model using an enum that implements the ParameterEnum interface.
    private enum GaussianParameter implements ParameterEnum {
        VARIANCE, MEAN
    }

    //We create a Modeller helper class to initialize the model state and specify the parameter samplers.
    private final class GaussianModeller {
        //Create fields in the Modeller for the model and samplers.
        private final ParameterizedModel<GaussianParameter, ParameterizedState<GaussianParameter>, GaussianDataCollection> model;
        private final ParameterSampler<Double, GaussianParameter, ParameterizedState<GaussianParameter>, GaussianDataCollection> varianceSampler;
        private final ParameterSampler<Double, GaussianParameter, ParameterizedState<GaussianParameter>, GaussianDataCollection> meanSampler;

        //Constructor for the Modeller takes as parameters all quantities needed to construct the ParameterizedState
        //(here, the initial variance and the initial mean) and the DataCollection (here, the list of datapoints).
        private GaussianModeller(final double varianceInitial, final double meanInitial, final List<Double> datapoints) {
            //Construct the initial ParameterizedState by passing a list of Parameters of mixed type to the constructor.
            //Initial values (and, implicitly, types) for each of the parameters are set here.
            final List<Parameter<GaussianParameter, ?>> initialParameters = Arrays.asList(
                    new Parameter<>(GaussianParameter.VARIANCE, varianceInitial),
                    new Parameter<>(GaussianParameter.MEAN, meanInitial));
            final ParameterizedState<GaussianParameter> initialState = new ParameterizedState<>(initialParameters);

            //Construct the GaussianDataCollection by passing a list of datapoints to the constructor.
            //Here, we pass 10000 datapoints, which were generated above,
            final GaussianDataCollection dataset = new GaussianDataCollection(datapoints);

            //Implement ParameterSamplers for each parameter by overriding sample().  This can be done via a lambda that takes
            //(rng, state, dataCollection) and returns a new sample of the parameter with type identical to that
            //specified during initialization above.

            //Sampler for the variance global parameter.  Assuming a uniform prior, the relevant log conditional PDF
            //is given by the log of the product of Gaussian likelihoods for each datapoint c_t:
            //      log[product_t variance^(-1/2) * exp(-(c_t - mean)^2 / (2 * variance))] + constant
            //which reduces to the form in code below.  Slice sampling is used here to generate a new sample,
            //but, in general, any method can be used; e.g., if the conditional PDF is from the exponential family,
            //one can simply sample directly from the corresponding Distribution from Apache Commons.
            varianceSampler = (rng, state, dataCollection) -> {
                final Function<Double, Double> logConditionalPDF =
                        newVariance -> -0.5 * Math.log(newVariance) * dataCollection.getDatapoints().size() +
                                dataCollection.getDatapoints().stream()
                                        .mapToDouble(c -> -normalTerm(c, state.get(GaussianParameter.MEAN, Double.class), newVariance))
                                        .sum();

                final SliceSampler sampler =
                        new SliceSampler(rng, logConditionalPDF, VARIANCE_MIN, VARIANCE_MAX, VARIANCE_WIDTH);
                return sampler.sample(state.get(GaussianParameter.VARIANCE, Double.class));
            };

            //Sampler for the mean global parameter.  Assuming a uniform prior, the relevant log conditional PDF
            //is given by the log of the product of Gaussian likelihoods for each datapoint c_t:
            //     log[product_t exp(-(c_t - mean)^2 / (2 * variance))] + constant
            //which reduces to the form in code below.
            meanSampler = (rng, state, dataCollection) -> {
                final Function<Double, Double> logConditionalPDF =
                        newMean -> dataCollection.getDatapoints().stream()
                                .mapToDouble(c -> -normalTerm(c, newMean, state.get(GaussianParameter.VARIANCE, Double.class)))
                                .sum();

                final SliceSampler sampler =
                        new SliceSampler(rng, logConditionalPDF, MEAN_WIDTH);
                return sampler.sample(state.get(GaussianParameter.MEAN, Double.class));
            };

            //Build the ParameterizedModel using the GibbsBuilder pattern.
            //Pass in the initial ParameterizedState and DataCollection, and specify the class of the ParameterizedState.
            //Add samplers for each of the parameters, with names matching those used in initialization.
            model = new ParameterizedModel.GibbsBuilder<>(initialState, dataset)
                    .addParameterSampler(GaussianParameter.VARIANCE, varianceSampler, Double.class)
                    .addParameterSampler(GaussianParameter.MEAN, meanSampler, Double.class)
                    .build();
        }
    }

    /**
     * Tests Bayesian inference of a Gaussian model via MCMC.  Recovery of input values for the variance and mean
     * global parameters is checked.  In particular, the mean and standard deviation of the posteriors for
     * both parameters must be recovered to within a relative error of 1% and 10%, respectively, in 250 samples
     * (after 250 burn-in samples have been discarded).
     */
    @Test
    public void testRunMCMCOnSingleGaussianModel() {
        //Create new instance of the Modeller helper class, passing all quantities needed to initialize state and data.
        final GaussianModeller modeller = new GaussianModeller(VARIANCE_INITIAL, MEAN_INITIAL, datapointsList);
        //Create a GibbsSampler, passing the total number of samples (including burn-in samples)
        //and the model held by the Modeller.
        final GibbsSampler<GaussianParameter, ParameterizedState<GaussianParameter>, GaussianDataCollection> gibbsSampler =
                new GibbsSampler<>(NUM_SAMPLES, modeller.model);
        //Run the MCMC.
        gibbsSampler.runMCMC();

        //Get the samples of each of the parameter posteriors (discarding burn-in samples) by passing the
        //parameter name, type, and burn-in number to the getSamples method.
        final double[] varianceSamples = Doubles.toArray(gibbsSampler.getSamples(GaussianParameter.VARIANCE, Double.class, NUM_BURN_IN));
        final double[] meanSamples = Doubles.toArray(gibbsSampler.getSamples(GaussianParameter.MEAN, Double.class, NUM_BURN_IN));

        //Check that the statistics---i.e., the means and standard deviations---of the posteriors
        //agree with those found by emcee/analytically to a relative error of 1% and 10%, respectively.
        final double variancePosteriorCenter = new Mean().evaluate(varianceSamples);
        final double variancePosteriorStandardDeviation = new StandardDeviation().evaluate(varianceSamples);
        Assert.assertEquals(relativeError(variancePosteriorCenter, VARIANCE_TRUTH),
                0., RELATIVE_ERROR_THRESHOLD_FOR_CENTERS);
        Assert.assertEquals(
                relativeError(variancePosteriorStandardDeviation, VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., RELATIVE_ERROR_THRESHOLD_FOR_STANDARD_DEVIATIONS);
        final double meanPosteriorCenter = new Mean().evaluate(meanSamples);
        final double meanPosteriorStandardDeviation = new StandardDeviation().evaluate(meanSamples);
        Assert.assertEquals(relativeError(meanPosteriorCenter, MEAN_TRUTH),
                0., RELATIVE_ERROR_THRESHOLD_FOR_CENTERS);
        Assert.assertEquals(
                relativeError(meanPosteriorStandardDeviation, MEAN_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., RELATIVE_ERROR_THRESHOLD_FOR_STANDARD_DEVIATIONS);
    }
}