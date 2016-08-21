package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.DecileCollection;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.PosteriorSummary;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tests for {@link TumorHeterogeneityData}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class TumorHeterogeneityDataUnitTest {
    private static final double LN2 = GATKProtectedMathUtils.LN2;
    private static final double INV_LN2 = GATKProtectedMathUtils.INV_LN2;
    private static final double EPSILON = TumorHeterogeneityUtils.EPSILON;

    private static final double REL_ERROR_THRESHOLD = 1E-4;

    private static final SimpleInterval DUMMY_INTERVAL = new SimpleInterval("1", 1, 100);
    private static final PosteriorSummary DUMMY_POSTERIOR_SUMMARY = new PosteriorSummary(Double.NaN, Double.NaN, Double.NaN);
    private static final DecileCollection DUMMY_DECILE_COLLECTION =
            new DecileCollection(Collections.singletonList(Double.NaN), DecileCollection.ConstructionMode.SAMPLES);

    private static final PloidyState NORMAL_PLOIDY_STATE = new PloidyState(1, 1);
    private static final PloidyStatePrior DUMMY_PLOIDY_STATE_PRIOR;
    private static final double DUMMY_PARAMETER = 1.;
    private static final TumorHeterogeneityPriorCollection DUMMY_PRIORS;

    static {
        DUMMY_POSTERIOR_SUMMARY.setDeciles(DUMMY_DECILE_COLLECTION);
        //need valid segment-mean posterior summary to construct TumorHeterogeneityData, but it is not used in tests
        final PosteriorSummary segmentMeanPosteriorSummary = new PosteriorSummary(0., -0.1, 0.1);
        segmentMeanPosteriorSummary.setDeciles(new DecileCollection(Arrays.asList(0., -0.1, 0.1), DecileCollection.ConstructionMode.SAMPLES));
        final ACNVModeledSegment segment1 = new ACNVModeledSegment(new SimpleInterval("1", 1, 25), segmentMeanPosteriorSummary, DUMMY_POSTERIOR_SUMMARY);
        final ACNVModeledSegment segment2 = new ACNVModeledSegment(new SimpleInterval("1", 26, 100), segmentMeanPosteriorSummary, DUMMY_POSTERIOR_SUMMARY);
        final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMap = new HashMap<>();
        unnormalizedLogProbabilityMassFunctionMap.put(NORMAL_PLOIDY_STATE, 0.);
        DUMMY_PLOIDY_STATE_PRIOR = new PloidyStatePrior(unnormalizedLogProbabilityMassFunctionMap);
        DUMMY_PRIORS = new TumorHeterogeneityPriorCollection(
                NORMAL_PLOIDY_STATE, DUMMY_PLOIDY_STATE_PRIOR,
                DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER);
    }

    @DataProvider(name = "dataNormal")
    public Object[][] dataNormal() {
        //mean and standard deviation
        return new Object[][]{
                {1., 1.},
                {5., 1.},
                {-5., 1.},
                {2., 0.1},
                {-2., 0.1},
                {10., 10.},
                {-10., 10.}
        };
    }

    @DataProvider(name = "dataBeta")
    public Object[][] dataBeta() {
        //alpha and beta
        return new Object[][]{
                {0.1, 0.1},
                {0.5, 0.1},
                {0.1, 0.5},
                {1., 1.},
                {5., 1.},
                {1., 5.},
                {10., 10.},
                {50., 10.},
                {10., 50.}
        };
    }

    @Test(dataProvider = "dataNormal")
    public void testFitNormalLogPDFToInnerDeciles(final double meanTruth, final double standardDeviationTruth) {
        //test the fitting of a normal distribution to log_2 copy-ratio deciles by setting minor-allele fraction deciles to NaN

        final double copyRatioNoiseConstant = 0.02;

        //calculate deciles for true copy-ratio posterior (normal)
        final NormalDistribution log2CopyRatioPosteriorDensityTruth =
                new NormalDistribution(Math.max(EPSILON, meanTruth + copyRatioNoiseConstant), standardDeviationTruth);
        final List<Double> log2CopyRatioDecilesTruth = IntStream.range(0, DecileCollection.NUM_DECILES).boxed()
                .map(i -> log2CopyRatioPosteriorDensityTruth.inverseCumulativeProbability(i / 10.)).collect(Collectors.toList());

        //construct TumorHeterogeneityPosteriorData from ACNVModeledSegment with true copy-ratio posterior deciles and NaN minor-allele fraction posterior
        final PosteriorSummary segmentMeanPosteriorSummary = new PosteriorSummary(
                meanTruth, meanTruth - 2 * standardDeviationTruth, meanTruth + 2 * standardDeviationTruth); //credible interval is not used in fit
        segmentMeanPosteriorSummary.setDeciles(new DecileCollection(log2CopyRatioDecilesTruth, DecileCollection.ConstructionMode.DECILES));
        final ACNVModeledSegment segment = new ACNVModeledSegment(DUMMY_INTERVAL, segmentMeanPosteriorSummary, DUMMY_POSTERIOR_SUMMARY);
        final TumorHeterogeneityData data = new TumorHeterogeneityData(Collections.singletonList(segment), DUMMY_PRIORS);

        //calculate expected log density from true distribution at point
        final double copyRatio = 1.;
        final double minorAlleleFraction = 0.25;
        final double log2CopyRatio = Math.log(Math.max(EPSILON, copyRatio + copyRatioNoiseConstant)) * INV_LN2;
        final double expectedLogDensity =
                Math.log(log2CopyRatioPosteriorDensityTruth.density(log2CopyRatio) * INV_LN2 / Math.max(EPSILON, copyRatio)) + LN2;

        //check against log density from TumorHeterogeneityData fit to true deciles
        final int segmentIndex = 0;
        final double resultLogDensity = data.logDensity(segmentIndex, copyRatio, minorAlleleFraction, copyRatioNoiseConstant);

        Assert.assertTrue(relativeError(resultLogDensity, expectedLogDensity) < REL_ERROR_THRESHOLD);
    }

    @Test(dataProvider = "dataBeta")
    public void testFitBetaLogPDFToInnerDeciles(final double alphaTruth, final double betaTruth) {
        //test the fitting of a beta distribution to minor-allele-fraction deciles,
        //assuming normal distribution is fit correctly to log_2 copy-ratio deciles

        final double copyRatio = 1.;
        final double standardDeviationTruth = 1.;
        final double copyRatioNoiseConstant = 0.02;

        //calculate deciles for true copy-ratio posterior (normal)
        final NormalDistribution log2CopyRatioPosteriorDensityTruth =
                new NormalDistribution(Math.max(EPSILON, copyRatio + copyRatioNoiseConstant), standardDeviationTruth);
        final List<Double> log2CopyRatioDecilesTruth = IntStream.range(0, DecileCollection.NUM_DECILES).boxed()
                .map(i -> log2CopyRatioPosteriorDensityTruth.inverseCumulativeProbability(i / 10.)).collect(Collectors.toList());

        //calculate deciles for true minor-allele--fraction posterior (beta)
        final BetaDistribution scaledMinorAlleleFractionPosteriorDensityTruth = new BetaDistribution(alphaTruth, betaTruth);
        final List<Double> minorAlleleFractionDecilesTruth = IntStream.range(0, DecileCollection.NUM_DECILES).boxed()
                .map(i -> scaledMinorAlleleFractionPosteriorDensityTruth.inverseCumulativeProbability(i / 10.) / 2.)
                .collect(Collectors.toList());

        //construct TumorHeterogeneityPosteriorData from ACNVModeledSegment with true posterior deciles
        final PosteriorSummary segmentMeanPosteriorSummary = new PosteriorSummary(
                copyRatio, copyRatio - 2 * standardDeviationTruth, copyRatio + 2 * standardDeviationTruth); //credible interval is not used in fit
        segmentMeanPosteriorSummary.setDeciles(new DecileCollection(log2CopyRatioDecilesTruth, DecileCollection.ConstructionMode.DECILES));
        final PosteriorSummary minorAlleleFractionPosteriorSummary = new PosteriorSummary(0.25, 0., 0.5);   //credible interval is not used in fit
        minorAlleleFractionPosteriorSummary.setDeciles(new DecileCollection(minorAlleleFractionDecilesTruth, DecileCollection.ConstructionMode.DECILES));
        final ACNVModeledSegment segment = new ACNVModeledSegment(DUMMY_INTERVAL, segmentMeanPosteriorSummary, minorAlleleFractionPosteriorSummary);
        final TumorHeterogeneityData data = new TumorHeterogeneityData(Collections.singletonList(segment), DUMMY_PRIORS);

        //calculate expected log density from true distribution at point
        final double minorAlleleFraction = 0.25;
        final double log2CopyRatio = Math.log(Math.max(EPSILON, copyRatio + copyRatioNoiseConstant)) * INV_LN2;
        final double expectedLogDensity =
                Math.log(log2CopyRatioPosteriorDensityTruth.density(log2CopyRatio) * INV_LN2 / Math.max(EPSILON, copyRatio)) +
                        Math.log(2. * scaledMinorAlleleFractionPosteriorDensityTruth.density(2. * minorAlleleFraction));

        //test log density at a point
        final int segmentIndex = 0;
        final double resultLogDensity = data.logDensity(segmentIndex, copyRatio, minorAlleleFraction, copyRatioNoiseConstant);

        Assert.assertTrue(relativeError(resultLogDensity, expectedLogDensity) < REL_ERROR_THRESHOLD);
    }

    @Test
    public void testCalculateFractionalSegmentLength() {
        //need valid segment-mean posterior summary to construct TumorHeterogeneityData, but it is not used in tests
        final PosteriorSummary segmentMeanPosteriorSummary = new PosteriorSummary(0., -0.1, 0.1);
        segmentMeanPosteriorSummary.setDeciles(new DecileCollection(Arrays.asList(0., -0.1, 0.1), DecileCollection.ConstructionMode.SAMPLES));
        final ACNVModeledSegment segment1 = new ACNVModeledSegment(new SimpleInterval("1", 1, 25), segmentMeanPosteriorSummary, DUMMY_POSTERIOR_SUMMARY);
        final ACNVModeledSegment segment2 = new ACNVModeledSegment(new SimpleInterval("1", 26, 100), segmentMeanPosteriorSummary, DUMMY_POSTERIOR_SUMMARY);
        final TumorHeterogeneityData data = new TumorHeterogeneityData(Arrays.asList(segment1, segment2), DUMMY_PRIORS);
        Assert.assertEquals(data.fractionalLength(0), 0.25);
        Assert.assertEquals(data.fractionalLength(1), 0.75);
    }

    private static double relativeError(final double x, final double xTrue) {
        return Math.abs((x - xTrue) / xTrue);
    }
}