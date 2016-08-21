package org.broadinstitute.hellbender.tools.exome.copyratio;

import htsjdk.samtools.util.Log;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link CopyRatioModeller}.
 * <p>
 *     Test data consisting of coverage and segment files for 100 segments with 100 targets each
 *     was generated using a python script.
 *     The global parameters determining the variance and the outlier probability were set to 1. and 0.05,
 *     respectively.  The segment means were drawn from Uniform(-5, 5), while outlier points were drawn from
 *     Uniform(-10, 10).
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioModellerUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome";

    private static final File COVERAGES_FILE = new File(TEST_SUB_DIR, "coverages-for-copy-ratio-modeller.tsv");
    private static final File SEGMENT_FILE = new File(TEST_SUB_DIR, "segments-for-copy-ratio-modeller.seg");
    private static final File MEANS_TRUTH_FILE = new File(TEST_SUB_DIR, "segment-means-truth-for-copy-ratio-modeller.txt");
    private static final File OUTLIER_INDICATORS_TRUTH_FILE = new File(TEST_SUB_DIR, "outlier-indicators-truth-for-copy-ratio-modeller.txt");

    private static final double CREDIBLE_INTERVAL_ALPHA = 0.32;

    private static final double VARIANCE_TRUTH = 1.;
    private static final double OUTLIER_PROBABILITY_TRUTH = 0.05;

    //truths for the posterior standard deviations are based on the standard deviations of the appropriate analytic
    //posteriors, scaled appropriately for the total number of coverages or the average number of coverages per segment
    private static final double MEAN_POSTERIOR_STANDARD_DEVIATION_MEAN_TRUTH = 0.1;                 //Gaussian with 100 points for each mean and unit variance gives 1 / sqrt(100)
    private static final double VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.014;                //inverse chi-squared with 100 DOF and variance = 1
    private static final double OUTLIER_PROBABILITY_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.0022;    //Beta for alpha ~ 516 (true # of outliers + prior alpha - 1),
                                                                                                    //         beta ~ 9580 (true # of non-outliers + prior beta - 1)

    //test specifications
    private static final double MULTIPLES_OF_SD_THRESHOLD = 1.5;
    private static final double RELATIVE_ERROR_THRESHOLD = 0.15;
    private static final double FRACTION_OF_OUTLIER_INDICATORS_CORRECT_THRESHOLD = 0.98;
    private static final int DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_1_SIGMA = 10;
    private static final int DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_2_SIGMA = 5;
    private static final int DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_3_SIGMA = 2;

    private static final int NUM_SAMPLES = 500;
    private static final int NUM_BURN_IN = 250;

    //Calculates relative error between x and xTrue, with respect to xTrue; used for checking statistics of
    //posterior samples below.
    private static double relativeError(final double x, final double xTrue) {
        return Math.abs((x - xTrue) / xTrue);
    }

    //Loads test data from files
    private static <T> List<T> loadList(final File file, final Function<String, T> parse) {
        try {
            return FileUtils.readLines(file).stream().map(parse).collect(Collectors.toList());
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(file, e);
        }
    }

    /**
     * Tests Bayesian inference of the copy-ratio model via MCMC.
     * <p>
     *     Recovery of input values for the variance and outlier-probability global parameters is checked.
     *     In particular, the true input value of the variance must fall within
     *     {@link CopyRatioModellerUnitTest#MULTIPLES_OF_SD_THRESHOLD}
     *     standard deviations of the posterior mean and the standard deviation of the posterior must agree
     *     with the analytic value to within a relative error of
     *     {@link CopyRatioModellerUnitTest#RELATIVE_ERROR_THRESHOLD} for 250 samples
     *     (after 250 burn-in samples have been discarded).  Similar criteria are applied
     *     to the recovery of the true input value for the outlier probability.
     * </p>
     * <p>
     *     Furthermore, the number of truth values for the segment-level means falling outside confidence intervals of
     *     1-sigma, 2-sigma, and 3-sigma given by the posteriors in each segment should be roughly consistent with
     *     a normal distribution (i.e., ~32, ~5, and ~0, respectively; we allow for errors of
     *     {@link CopyRatioModellerUnitTest#DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_1_SIGMA},
     *     {@link CopyRatioModellerUnitTest#DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_2_SIGMA}, and
     *     {@link CopyRatioModellerUnitTest#DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_3_SIGMA}, respectively).
     *     The mean of the standard deviations of the posteriors for the segment-level means should also be
     *     recovered to within a relative error of {@link CopyRatioModellerUnitTest#RELATIVE_ERROR_THRESHOLD}.
     * </p>
     * <p>
     *     Finally, the recovered values for the latent outlier-indicator parameters should agree with those used to
     *     generate the data.  For each indicator, the recovered value (i.e., outlier or non-outlier) is taken to be
     *     that given by the majority of posterior samples.  We require that at least
     *     {@link CopyRatioModellerUnitTest#FRACTION_OF_OUTLIER_INDICATORS_CORRECT_THRESHOLD}
     *     of the 10000 indicators are recovered correctly.
     * </p>
     * <p>
     *     With these specifications, this unit test is not overly brittle (i.e., it should pass for a large majority
     *     of randomly generated data sets), but it is still brittle enough to check for correctness of the sampling
     *     (for example, specifying a sufficiently incorrect likelihood will cause the test to fail).
     * </p>
     */
    @Test
    public void testRunMCMCOnCopyRatioSegmentedGenome() throws IOException {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);

        //load data (coverages and number of targets in each segment)
        final ReadCountCollection coverage = ReadCountCollectionUtils.parse(COVERAGES_FILE);
        final Genome genome = new Genome(coverage, Collections.emptyList()); //Genome with no SNPs
        final SegmentedGenome segmentedGenome = new SegmentedGenome(SEGMENT_FILE, genome);

        //run MCMC
        final CopyRatioModeller modeller = new CopyRatioModeller(segmentedGenome);
        modeller.fitMCMC(NUM_SAMPLES, NUM_BURN_IN);

        //check statistics of global-parameter posterior samples (i.e., posterior mode and standard deviation)
        final Map<CopyRatioParameter, PosteriorSummary> globalParameterPosteriorSummaries =
                modeller.getGlobalParameterPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);

        final PosteriorSummary variancePosteriorSummary = globalParameterPosteriorSummaries.get(CopyRatioParameter.VARIANCE);
        final double variancePosteriorCenter = variancePosteriorSummary.getCenter();
        final double variancePosteriorStandardDeviation = (variancePosteriorSummary.getUpper() - variancePosteriorSummary.getLower()) / 2;
        Assert.assertEquals(Math.abs(variancePosteriorCenter - VARIANCE_TRUTH),
                0., MULTIPLES_OF_SD_THRESHOLD * VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH);
        Assert.assertEquals(relativeError(variancePosteriorStandardDeviation, VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., RELATIVE_ERROR_THRESHOLD);

        final PosteriorSummary outlierProbabilityPosteriorSummary = globalParameterPosteriorSummaries.get(CopyRatioParameter.OUTLIER_PROBABILITY);
        final double outlierProbabilityPosteriorCenter = outlierProbabilityPosteriorSummary.getCenter();
        final double outlierProbabilityPosteriorStandardDeviation = (outlierProbabilityPosteriorSummary.getUpper() - outlierProbabilityPosteriorSummary.getLower()) / 2;
        Assert.assertEquals(Math.abs(outlierProbabilityPosteriorCenter - OUTLIER_PROBABILITY_TRUTH),
                0., MULTIPLES_OF_SD_THRESHOLD * OUTLIER_PROBABILITY_POSTERIOR_STANDARD_DEVIATION_TRUTH);
        Assert.assertEquals(relativeError(outlierProbabilityPosteriorStandardDeviation,
                OUTLIER_PROBABILITY_POSTERIOR_STANDARD_DEVIATION_TRUTH), 0., RELATIVE_ERROR_THRESHOLD);

        //check statistics of segment-mean posterior samples (i.e., posterior means and standard deviations)
        final List<Double> meansTruth = loadList(MEANS_TRUTH_FILE, Double::parseDouble);
        int numMeansOutsideOneSigma = 0;
        int numMeansOutsideTwoSigma = 0;
        int numMeansOutsideThreeSigma = 0;
        final int numSegments = meansTruth.size();
        //segment-mean posteriors are expected to be Gaussian, so PosteriorSummary for
        // {@link CopyRatioModellerUnitTest#CREDIBLE_INTERVAL_ALPHA}=0.32 is
        //(posterior mean, posterior mean - posterior standard devation, posterior mean + posterior standard deviation)
        final List<PosteriorSummary> meanPosteriorSummaries =
                modeller.getSegmentMeansPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        final double[] meanPosteriorStandardDeviations = new double[numSegments];
        for (int segment = 0; segment < numSegments; segment++) {
            final double meanPosteriorCenter = meanPosteriorSummaries.get(segment).getCenter();
            final double meanPosteriorStandardDeviation =
                    (meanPosteriorSummaries.get(segment).getUpper() - meanPosteriorSummaries.get(segment).getLower()) / 2.;
            meanPosteriorStandardDeviations[segment] = meanPosteriorStandardDeviation;
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
        final double meanPosteriorStandardDeviationsMean = new Mean().evaluate(meanPosteriorStandardDeviations);
        Assert.assertEquals(numMeansOutsideOneSigma, 100 - 68, DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_1_SIGMA);
        Assert.assertEquals(numMeansOutsideTwoSigma, 100 - 95, DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_2_SIGMA);
        Assert.assertTrue(numMeansOutsideThreeSigma <= DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_3_SIGMA);
        Assert.assertEquals(relativeError(meanPosteriorStandardDeviationsMean, MEAN_POSTERIOR_STANDARD_DEVIATION_MEAN_TRUTH),
                0., RELATIVE_ERROR_THRESHOLD);

        //check accuracy of latent outlier-indicator posterior samples
        final List<CopyRatioState.OutlierIndicators> outlierIndicatorSamples =
                modeller.getOutlierIndicatorsSamples();
        int numIndicatorsCorrect = 0;
        final int numIndicatorSamples = outlierIndicatorSamples.size();
        final List<Integer> outlierIndicatorsTruthAsInt = loadList(OUTLIER_INDICATORS_TRUTH_FILE, Integer::parseInt);
        final List<Boolean> outlierIndicatorsTruth =
                outlierIndicatorsTruthAsInt.stream().map(i -> i == 1).collect(Collectors.toList());
        for (int target = 0; target < coverage.targets().size(); target++) {
            int numSamplesOutliers = 0;
            for (final CopyRatioState.OutlierIndicators sample : outlierIndicatorSamples) {
                if (sample.get(target)) {
                    numSamplesOutliers++;
                }
            }
            //take predicted state of indicator to be given by the majority of samples
            if ((numSamplesOutliers >= numIndicatorSamples / 2.) == outlierIndicatorsTruth.get(target)) {
                numIndicatorsCorrect++;
            }
        }
        final double fractionOfOutlierIndicatorsCorrect = (double) numIndicatorsCorrect / coverage.targets().size();
        Assert.assertTrue(fractionOfOutlierIndicatorsCorrect >= FRACTION_OF_OUTLIER_INDICATORS_CORRECT_THRESHOLD);
    }
}