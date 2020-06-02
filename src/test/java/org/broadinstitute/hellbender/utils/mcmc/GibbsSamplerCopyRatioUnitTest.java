package org.broadinstitute.hellbender.utils.mcmc;

import com.google.common.primitives.Doubles;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Unit test for {@link GibbsSampler}.  Demonstrates application of {@link GibbsSampler} to a hierarchical
 * {@link ParameterizedModel} that is specified using helper classes that extend
 * {@link ParameterizedState} or implement {@link DataCollection}.
 * <p>
 *     Test performs Bayesian inference of a simple copy-ratio model with 1 global and 100 segment-level parameters.
 *     We consider an exome with 100 segments, each of which has a number of targets that is drawn from a
 *     Poisson distribution with mean 100 in the generated data set.
 * </p>
 * <p>
 *     Data consists of a list of coverages at each target and a list of the number of targets in each segment.
 *     Coverages are assumed to be drawn from a normal distribution in each segment.
 * </p>
 * <p>
 *     The mean of each normal distribution is given by a segment-level parameter;
 *     the distribution of means across segments is taken to be uniform in (0, 10) in the generated data.
 * </p>
 * <p>
 *     The variance of each distribution is set by the global parameter, which is taken to be 1 in the generated data.
 * </p>
 * <p>
 *     Success of the test is determined by recovery of the input segment means and the global variance,
 *     as well as agreement of the standard deviations of the parameter posteriors with those given by both the
 *     python package emcee (see http://dan.iel.fm/emcee for details) and numerical evaluation in Mathematica
 *     of the analytic forms of the posteriors.
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class GibbsSamplerCopyRatioUnitTest extends GATKBaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/utils/mcmc/";

    private static final File COVERAGES_FILE =
            new File(TEST_SUB_DIR, "coverages-for-gibbs-sampler-copy-ratio-test.txt");
    private static final File NUM_TARGETS_PER_SEGMENT_FILE =
            new File(TEST_SUB_DIR, "number-of-targets-per-segment-for-gibbs-sampler-copy-ratio-test.txt");
    private static final File MEANS_TRUTH_FILE =
            new File(TEST_SUB_DIR, "means-truth-for-gibbs-sampler-copy-ratio-test.txt");

    private static final double VARIANCE_MIN = 0.;
    private static final double VARIANCE_MAX = Double.POSITIVE_INFINITY;
    private static final double VARIANCE_WIDTH = 0.1;
    private static final double VARIANCE_INITIAL = 5.;
    private static final double VARIANCE_TRUTH = 1.;
    private static final double VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.014;    //checked with emcee & analytic result

    private static final double MEAN_WIDTH = 0.1;
    private static final double MEAN_INITIAL = 5.;
    private static final double MEAN_POSTERIOR_STANDARD_DEVIATION_MEAN_TRUTH = 0.1;     //checked with emcee & analytic result

    private static final int NUM_SAMPLES = 1000;
    private static final int NUM_BURN_IN = 500;

    //test specifications
    private static final double RELATIVE_ERROR_THRESHOLD_FOR_CENTERS = 0.01;
    private static final double RELATIVE_ERROR_THRESHOLD_FOR_STANDARD_DEVIATIONS = 0.05;
    private static final int DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_1_SIGMA = 10;
    private static final int DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_2_SIGMA = 5;
    private static final int DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_3_SIGMA = 2;

    //Calculates the exponent for a normal distribution; used in log-likelihood calculation below.
    private static double normalTerm(final double quantity, final double mean, final double variance) {
        return (quantity - mean) * (quantity - mean) / (2 * variance);
    }

    //Calculates relative error between x and xTrue, with respect to xTrue; used for checking statistics of
    //posterior samples below.
    private static double relativeError(final double x, final double xTrue) {
        return Math.abs((x - xTrue) / xTrue);
    }

    //Loads test data from file
    private static <T> List<T> loadList(final File file, final Function<String, T> parse) {
        try {
            return FileUtils.readLines(file, StandardCharsets.UTF_8).stream().map(parse).collect(Collectors.toList());
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(file.getAbsolutePath(), e);
        }
    }

    //In contrast to GibbsSamplerSingleGaussianUnitTest, our model here is complicated enough to warrant
    //creating helper classes to make defining the parameter ParameterSamplers and referring to the data easier.

    //First, we enumerate the parameters of the model using an enum that implements the ParameterEnum interface.
    private enum CopyRatioParameter implements ParameterEnum {
        VARIANCE, SEGMENT_MEANS
    }

    //Next, we create a subclass SegmentMeans.  This is a simple class that represents a "parameter block" of
    //100 parameters, one for each segment-level mean.  NOTE: We don't actually enforce that the list of means
    //is of length 100 or perform other such checks, although one could if desired.
    private final class SegmentMeans extends ArrayList<Double> {
        private static final long serialVersionUID = 147369L;
        public SegmentMeans(final List<Double> segmentMeans) {
            super(new ArrayList<>(segmentMeans));
        }
    }

    //Using a helper class that extends ParameterizedState, we can make it easier to refer to the
    //100 segment-level mean parameters and the global variance parameter.  For example,
    //this will allow us to refer to the variance as simply "state.variance()", rather than
    //"state.get(CopyRatioParameter.VARIANCE, Double.class)".  Likewise, we will easily be able to refer
    //to the coverages per segment.

    //We thus create a subclass CopyRatioState.  This is a ParameterizedState<CopyRatioParameter> with
    //two parameters enumerated by <CopyRatioParameter>.  The first parameter is a double for the global variance,
    //while the second parameter is a "parameter block" object named SegmentMeans that holds all of the segment-level means.
    //We see that we can easily use such blocks to create a hierarchical model.
    private final class CopyRatioState extends ParameterizedState<CopyRatioParameter> {
        //this constructor conveniently creates a CopyRatioState from a given variance and mean (used for initialization)
        //all segment means are set to the given value
        public CopyRatioState(final double variance, final double mean, final int numSegments) {
            super(Arrays.asList(
                    new Parameter<>(CopyRatioParameter.VARIANCE, variance),
                    new Parameter<>(CopyRatioParameter.SEGMENT_MEANS, new SegmentMeans(Collections.nCopies(numSegments, mean)))));
        }

        //this getter allows us to conveniently refer to the variance when constructing the ParameterSamplers below
        public double variance() {
            return get(CopyRatioParameter.VARIANCE, Double.class);
        }

        //this getter allows us to conveniently refer to the segment means when constructing the ParameterSamplers below
        public double meanInSegment(final int segment) {
            return get(CopyRatioParameter.SEGMENT_MEANS, SegmentMeans.class).get(segment);
        }
    }

    //Likewise, we will create a CopyRatioDataCollection (which implements DataCollection) to make it easier to refer
    //to the coverages per segment when constructing the ParameterSamplers below.
    private final class CopyRatioDataCollection implements DataCollection {
        private final int numSegments;

        private final List<List<Double>> coveragesPerSegment = new ArrayList<>();

        public CopyRatioDataCollection(final List<Double> coverages, final List<Integer> numTargetsPerSegment) {
            numSegments = numTargetsPerSegment.size();
            //partition coverages by segment
            int startTarget = 0;
            for (int segment = 0; segment < numSegments; segment++) {
                final int numTargetsInSegment = numTargetsPerSegment.get(segment);
                coveragesPerSegment.add(coverages.subList(startTarget, startTarget + numTargetsInSegment));
                startTarget += numTargetsInSegment;
            }
        }

        //constructor that loads the coverages and number of targets per segment from files
        public CopyRatioDataCollection(final File coveragesFile, final File numTargetsPerSegmentFile) {
            this(loadList(coveragesFile, Double::parseDouble), loadList(numTargetsPerSegmentFile, Integer::parseInt));
        }

        //this getter allows us to conveniently refer to the coverages per segment when constructing the ParameterSamplers below
        public List<Double> getCoveragesInSegment(final int segment) {
            return coveragesPerSegment.get(segment);
        }

        //this getter allows us to conveniently refer to the number of targets per segment when constructing the ParameterSamplers below
        public int getNumTargetsInSegment(final int segment) {
            return coveragesPerSegment.get(segment).size();
        }
    }

    //We create a Modeller helper class to initialize the model state and specify the parameter samplers.
    private final class CopyRatioModeller {
        //Create fields in the Modeller for the model and samplers.
        private final ParameterizedModel<CopyRatioParameter, CopyRatioState, CopyRatioDataCollection> model;
        private final ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioDataCollection> varianceSampler;
        private final ParameterSampler<SegmentMeans, CopyRatioParameter, CopyRatioState, CopyRatioDataCollection> meansSampler;  //note that this returns a SegmentMeans sample

        //Constructor for the Modeller takes as parameters all quantities needed to construct the ParameterizedState
        //(here, the initial variance and the initial mean) and the DataCollection.
        public CopyRatioModeller(final double initialVariance, final double initalMean,
                                 final CopyRatioDataCollection dataset) {
            //Construct the initial CopyRatioState by passing all necesssary quantities.
            final CopyRatioState initialState = new CopyRatioState(initialVariance, initalMean, dataset.numSegments);

            //Implement ParameterSamplers for each parameter by overriding sample().  This can be done via a lambda that takes
            //(rng, state, dataCollection) and returns a new sample of the parameter with type identical to that
            //specified during initialization above.

            //Sampler for the variance global parameter.  Assuming a uniform prior, the relevant logConditionalPDF
            //is given by the log of the product of Gaussian likelihoods for each target t:
            //      log[product_t variance^(-1/2) * exp(-(coverage_t - mean_t)^2 / (2 * variance))] + constant
            //which reduces to the form in code below.  Note that mean_t is identical for all targets in a segment.
            varianceSampler = (rng, state, dataCollection) -> {
                final Function<Double, Double> logConditionalPDF = newVariance -> {
                    double ll = 0.;
                    for (int segment = 0; segment < dataCollection.numSegments; segment++) {
                        final double meanInSegment = state.meanInSegment(segment);
                        ll += -0.5 * Math.log(newVariance) * dataCollection.getNumTargetsInSegment(segment) +
                                dataCollection.getCoveragesInSegment(segment).stream()
                                        .mapToDouble(c -> -normalTerm(c, meanInSegment, newVariance)).sum();
                    }
                    return ll;
                };

                final SliceSampler sampler = new SliceSampler(rng, logConditionalPDF, VARIANCE_MIN, VARIANCE_MAX, VARIANCE_WIDTH);
                return sampler.sample(state.variance());
            };

            //Sampler for the segment-level mean parameters.  Assuming uniform priors, for each segment s, the relevant
            //logConditionalPDF is given by the log of the product of Gaussian likelihoods for all targets t
            //in that segment:
            //     log[product_{t in s} exp(-(coverage_t - mean_s)^2 / (2 * variance))] + constant
            //which reduces to the form in code below.  Note that a lambda is specified for each segment,
            //so that meansSampler.sample() returns a SegmentMeans object that contains the means for all segments.
            //NOTE: It is possible to shoot yourself in the foot here if you incorrectly specify the sampler!
            //For example, note that there is no check that the number of segments sampled here matches the
            //number of segments given during initialization.  Such consistency checks can be added (e.g., in the
            //constructor of SegmentMeans), if desired, but may impact performance.
            meansSampler = (rng, state, dataCollection) -> {
                final List<Double> means = new ArrayList<>();
                for (int segment = 0; segment < dataCollection.numSegments; segment++) {
                    final List<Double> coveragesInSegment = dataCollection.getCoveragesInSegment(segment);
                    final Function<Double, Double> logConditionalPDF =
                            newMean -> coveragesInSegment.stream()
                                    .mapToDouble(c -> -normalTerm(c, newMean, state.variance()))
                                    .sum();

                    final SliceSampler sampler = new SliceSampler(rng, logConditionalPDF, MEAN_WIDTH);
                    means.add(sampler.sample(state.meanInSegment(segment)));
                }
                return new SegmentMeans(means);
            };

            //Build the ParameterizedModel using the GibbsBuilder pattern.
            //Pass in the initial CopyRatioState and CopyRatioDataCollection, and specify the class of the CopyRatioState.
            //Add samplers for each of the parameters, with names matching those used in initialization.
            model = new ParameterizedModel.GibbsBuilder<>(initialState, dataset)
                    .addParameterSampler(CopyRatioParameter.VARIANCE, varianceSampler, Double.class)
                    .addParameterSampler(CopyRatioParameter.SEGMENT_MEANS, meansSampler, SegmentMeans.class)
                    .build();
        }

        //Constructor for the Modeller takes as parameters all quantities needed to construct the ParameterizedState
        //(here, the initial variance and the initial mean) and the DataCollection
        //(here, the coverage and target-number files).
        public CopyRatioModeller(final double initialVariance, final double initalMean,
                                 final File coveragesFile, final File numTargetsPerSegmentFile) {
            this(initialVariance, initalMean, new CopyRatioDataCollection(coveragesFile, numTargetsPerSegmentFile));
        }
    }

    /**
     * Tests Bayesian inference of a toy copy-ratio model via MCMC.
     * <p>
     *     Recovery of input values for the variance global parameter and the segment-level mean parameters is checked.
     *     In particular, the mean and standard deviation of the posterior for the variance must be recovered to within
     *     a relative error of 1% and 5%, respectively, in 500 samples (after 250 burn-in samples have been discarded).
     * </p>
     * <p>
     *     Furthermore, the number of truth values for the segment-level means falling outside confidence intervals of
     *     1-sigma, 2-sigma, and 3-sigma given by the posteriors in each segment should be roughly consistent with
     *     a normal distribution (i.e., ~32, ~5, and ~0, respectively; we allow for errors of 10, 5, and 2).
     *     Finally, the mean of the standard deviations of the posteriors for the segment-level means should be
     *     recovered to within a relative error of 5%.
     * </p>
     * <p>
     *     With these specifications, this unit test is not overly brittle (i.e., it should pass for a large majority
     *     of randomly generated data sets), but it is still brittle enough to check for correctness of the sampling
     *     (for example, specifying a sufficiently incorrect likelihood will cause the test to fail).
     * </p>
     */
    @Test
    public void testRunMCMCOnCopyRatioSegmentedGenome() {
        //Create new instance of the Modeller helper class, passing all quantities needed to initialize state and data.
        final CopyRatioModeller modeller =
                new CopyRatioModeller(VARIANCE_INITIAL, MEAN_INITIAL, COVERAGES_FILE, NUM_TARGETS_PER_SEGMENT_FILE);
        //Create a GibbsSampler, passing the total number of samples (including burn-in samples)
        //and the model held by the Modeller.
        final GibbsSampler<CopyRatioParameter, CopyRatioState, CopyRatioDataCollection> gibbsSampler =
                new GibbsSampler<>(NUM_SAMPLES, modeller.model);
        //Run the MCMC.
        gibbsSampler.runMCMC();

        //Check that the statistics---i.e., the mean and standard deviation---of the variance posterior
        //agree with those found by emcee/analytically to a relative error of 1% and 5%, respectively.
        final double[] varianceSamples =
                Doubles.toArray(gibbsSampler.getSamples(CopyRatioParameter.VARIANCE, Double.class, NUM_BURN_IN));
        final double variancePosteriorCenter = new Mean().evaluate(varianceSamples);
        final double variancePosteriorStandardDeviation = new StandardDeviation().evaluate(varianceSamples);
        Assert.assertEquals(relativeError(variancePosteriorCenter, VARIANCE_TRUTH),
                0., RELATIVE_ERROR_THRESHOLD_FOR_CENTERS);
        Assert.assertEquals(relativeError(variancePosteriorStandardDeviation, VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., RELATIVE_ERROR_THRESHOLD_FOR_STANDARD_DEVIATIONS);
        //Check statistics---i.e., the mean and standard deviation---of the segment-level mean posteriors.
        //In particular, check that the number of segments where the true mean falls outside confidence intervals
        //is roughly consistent with a normal distribution.
        final List<Double> meansTruth = loadList(MEANS_TRUTH_FILE, Double::parseDouble);
        final int numSegments = meansTruth.size();
        final List<SegmentMeans> meansSamples =
                gibbsSampler.getSamples(CopyRatioParameter.SEGMENT_MEANS, SegmentMeans.class, NUM_BURN_IN);
        int numMeansOutsideOneSigma = 0;
        int numMeansOutsideTwoSigma = 0;
        int numMeansOutsideThreeSigma = 0;
        final List<Double> meanPosteriorStandardDeviations = new ArrayList<>();
        for (int segment = 0; segment < numSegments; segment++) {
            final int j = segment;
            final double[] meanInSegmentSamples =
                    Doubles.toArray(meansSamples.stream().map(s -> s.get(j)).collect(Collectors.toList()));
            final double meanPosteriorCenter = new Mean().evaluate(meanInSegmentSamples);
            final double meanPosteriorStandardDeviation =
                    new StandardDeviation().evaluate(meanInSegmentSamples);
            meanPosteriorStandardDeviations.add(meanPosteriorStandardDeviation);
            final double absoluteDifferenceFromTruth = Math.abs(meanPosteriorCenter - meansTruth.get(segment));
            if (absoluteDifferenceFromTruth > meanPosteriorStandardDeviation) {
                numMeansOutsideOneSigma++;
            }
            if (absoluteDifferenceFromTruth > 2 * meanPosteriorStandardDeviation) {
                numMeansOutsideTwoSigma++;
            }
            if (absoluteDifferenceFromTruth > 3 * meanPosteriorStandardDeviation) {
                numMeansOutsideThreeSigma++;
            }
        }
        final double meanPosteriorStandardDeviationsMean =
                new Mean().evaluate(Doubles.toArray(meanPosteriorStandardDeviations));
        Assert.assertEquals(numMeansOutsideOneSigma, 100 - 68, DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_1_SIGMA);
        Assert.assertEquals(numMeansOutsideTwoSigma, 100 - 95, DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_2_SIGMA);
        Assert.assertTrue(numMeansOutsideThreeSigma <= DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_3_SIGMA);
        Assert.assertEquals(
                relativeError(meanPosteriorStandardDeviationsMean, MEAN_POSTERIOR_STANDARD_DEVIATION_MEAN_TRUTH),
                0., RELATIVE_ERROR_THRESHOLD_FOR_STANDARD_DEVIATIONS);
    }

    //constants for testing exceptions
    private final CopyRatioState SIMPLE_STATE =
            new CopyRatioState(VARIANCE_INITIAL, MEAN_INITIAL, 1);
    private final CopyRatioDataCollection SIMPLE_DATA =
            new CopyRatioDataCollection(Collections.singletonList(1.), Collections.singletonList(1));
    private final ParameterizedModel.GibbsBuilder<CopyRatioParameter, CopyRatioState, CopyRatioDataCollection> EXCEPTION_BUILDER =
            new ParameterizedModel.GibbsBuilder<>(SIMPLE_STATE, SIMPLE_DATA);
    private final CopyRatioModeller EXCEPTION_MODELLER =
            new CopyRatioModeller(VARIANCE_INITIAL, MEAN_INITIAL, SIMPLE_DATA);
    private final ParameterizedModel<CopyRatioParameter, CopyRatioState, CopyRatioDataCollection> EXCEPTION_MODEL = EXCEPTION_BUILDER
            .addParameterSampler(CopyRatioParameter.VARIANCE, EXCEPTION_MODELLER.varianceSampler, Double.class)
            .addParameterSampler(CopyRatioParameter.SEGMENT_MEANS, EXCEPTION_MODELLER.meansSampler, SegmentMeans.class)
            .build();

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativeNumSamplesException() {
        new GibbsSampler<>(-1, EXCEPTION_MODEL);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativeNumSamplesPerLogEntryException() {
        final GibbsSampler<CopyRatioParameter, CopyRatioState, CopyRatioDataCollection> gibbsSampler =
                new GibbsSampler<>(1, EXCEPTION_MODEL);
        gibbsSampler.setNumSamplesPerLogEntry(1);
        gibbsSampler.setNumSamplesPerLogEntry(-1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativeNumBurnInException() {
        final GibbsSampler<CopyRatioParameter, CopyRatioState, CopyRatioDataCollection> gibbsSampler =
                new GibbsSampler<>(1, EXCEPTION_MODEL);
        gibbsSampler.getSamples(CopyRatioParameter.VARIANCE, Double.class, -1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNumBurnInExceedsNumSamplesException() {
        final GibbsSampler<CopyRatioParameter, CopyRatioState, CopyRatioDataCollection> gibbsSampler =
                new GibbsSampler<>(5, EXCEPTION_MODEL);
        gibbsSampler.getSamples(CopyRatioParameter.VARIANCE, Double.class, 0);
        gibbsSampler.getSamples(CopyRatioParameter.VARIANCE, Double.class, 10);
    }
}