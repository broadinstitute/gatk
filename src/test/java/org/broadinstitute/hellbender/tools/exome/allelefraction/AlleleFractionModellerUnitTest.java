package org.broadinstitute.hellbender.tools.exome.allelefraction;

import htsjdk.samtools.util.Log;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Tests the MCMC inference of the {@link AlleleFractionModeller}.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionModellerUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    // test data is a "normal" PON generated from 50 normals simulated from the allele-fraction model with alpha = 65 and beta = 60
    // and a PON with "bad SNPs" described below
    private static final File ALLELIC_PON_NORMAL_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-pon-normal.tsv");
    private static final File ALLELIC_PON_WITH_BAD_SNPS_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-pon-bad.tsv");

    private static final File SAMPLE_NORMAL_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-sample-normal.tsv");
    private static final File SAMPLE_WITH_BAD_SNPS_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-sample-bad.tsv");
    private static final File SAMPLE_WITH_EVENT_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-sample-event.tsv");

    private static final File SEGMENTS_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-segments.seg");

    private static final double CREDIBLE_INTERVAL_ALPHA = 0.05;

    /**
     * Test MCMC inference on simulated data.  Note that hyperparameter values used to generate the data should be recovered
     * along with outlier probability and minor fractions.
     */
    @Test
    public void testMCMCWithoutAllelicPON() {
        final double meanBiasSimulated = 1.2;
        final double biasVarianceSimulated = 0.04;
        testMCMC(meanBiasSimulated, biasVarianceSimulated, meanBiasSimulated, biasVarianceSimulated, AllelicPanelOfNormals.EMPTY_PON);
    }

    /**
     * Test MCMC inference on simulated data using an allelic PON.  Note that these MCMC tests were written to use
     * simulated hets before the allelic PON was introduced.  Rather than generate a simulated PON on the fly,
     * we simply use a fixed simulated PON loaded from a file and check that its MLE hyperparameters are "sampled"
     * correctly by simply taking the MLE PON values---i.e., the PON does not actually cover the simulated sites and
     * hence is not used to correct reference bias in the simulated data in any way.
     * This latter functionality is tested on fixed data loaded from files in
     * {@link AlleleFractionModellerUnitTest#testBiasCorrection} instead.
     */
    @Test
    public void testMCMCWithAllelicPON() {
        final double meanBiasSimulated = 1.2;
        final double biasVarianceSimulated = 0.04;
        final double meanBiasOfPON = 1.083;         // PON generated with alpha = 65
        final double biasVarianceOfPON = 0.0181;    // PON generated with beta = 60
        final AllelicPanelOfNormals allelicPON = new AllelicPanelOfNormals(ALLELIC_PON_NORMAL_FILE);
        testMCMC(meanBiasSimulated, biasVarianceSimulated, meanBiasOfPON, biasVarianceOfPON, allelicPON);
    }

    private void testMCMC(final double meanBiasSimulated, final double biasVarianceSimulated,
                          final double meanBiasExpected, final double biasVarianceExpected,
                          final AllelicPanelOfNormals allelicPON) {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final int numSamples = 150;
        final int numBurnIn = 50;

        final double averageHetsPerSegment = 50;
        final int numSegments = 100;
        final int averageDepth = 50;

        final double outlierProbability = 0.02;

        // note: the following tolerances could actually be made much smaller if we used more segments and/or
        // more hets -- most of the error is the sampling error of a finite simulated data set, not numerical error of MCMC
        final double minorFractionTolerance = 0.02;
        final double meanBiasTolerance = 0.02;
        final double biasVarianceTolerance = 0.01;
        final double outlierProbabilityTolerance = 0.02;
        final AlleleFractionSimulatedData simulatedData = new AlleleFractionSimulatedData(averageHetsPerSegment, numSegments,
                averageDepth, meanBiasSimulated, biasVarianceSimulated, outlierProbability);

        final AlleleFractionModeller modeller = new AlleleFractionModeller(simulatedData.getSegmentedGenome(), allelicPON);
        modeller.fitMCMC(numSamples, numBurnIn);

        final List<Double> meanBiasSamples = modeller.getmeanBiasSamples();
        Assert.assertEquals(meanBiasSamples.size(), numSamples - numBurnIn);

        final List<Double> biasVarianceSamples = modeller.getBiasVarianceSamples();
        Assert.assertEquals(biasVarianceSamples.size(), numSamples - numBurnIn);

        final List<Double> outlierProbabilitySamples = modeller.getOutlierProbabilitySamples();
        Assert.assertEquals(outlierProbabilitySamples.size(), numSamples - numBurnIn);

        final List<AlleleFractionState.MinorFractions> minorFractionsSamples = modeller.getMinorFractionsSamples();
        Assert.assertEquals(minorFractionsSamples.size(), numSamples - numBurnIn);
        for (final AlleleFractionState.MinorFractions sample : minorFractionsSamples) {
            Assert.assertEquals(sample.size(), numSegments);
        }

        final List<List<Double>> minorFractionsSamplesBySegment = modeller.getMinorFractionSamplesBySegment();

        final double mcmcMeanBias = meanBiasSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final double mcmcBiasVariance = biasVarianceSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final double mcmcOutlierProbability = outlierProbabilitySamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final List<Double> mcmcMinorFractions = minorFractionsSamplesBySegment
                .stream().map(list -> list.stream().mapToDouble(x -> x).average().getAsDouble())
                .collect(Collectors.toList());

        double totalSegmentError = 0.0;
        for (int segment = 0; segment < numSegments; segment++) {
            totalSegmentError += Math.abs(mcmcMinorFractions.get(segment) - simulatedData.getTrueState().segmentMinorFraction(segment));
        }

        Assert.assertEquals(mcmcMeanBias, meanBiasExpected, meanBiasTolerance);
        Assert.assertEquals(mcmcBiasVariance, biasVarianceExpected, biasVarianceTolerance);
        Assert.assertEquals(mcmcOutlierProbability, outlierProbability, outlierProbabilityTolerance);
        Assert.assertEquals(totalSegmentError / numSegments, 0.0, minorFractionTolerance);

        //test posterior summaries
        final Map<AlleleFractionParameter, PosteriorSummary> globalParameterPosteriorSummaries =
                modeller.getGlobalParameterPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);

        final PosteriorSummary meanBiasPosteriorSummary = globalParameterPosteriorSummaries.get(AlleleFractionParameter.MEAN_BIAS);
        final double meanBiasPosteriorCenter = meanBiasPosteriorSummary.getCenter();
        Assert.assertEquals(meanBiasPosteriorCenter, meanBiasExpected, meanBiasTolerance);

        final PosteriorSummary biasVariancePosteriorSummary = globalParameterPosteriorSummaries.get(AlleleFractionParameter.BIAS_VARIANCE);
        final double biasVariancePosteriorCenter = biasVariancePosteriorSummary.getCenter();
        Assert.assertEquals(biasVariancePosteriorCenter, biasVarianceExpected, biasVarianceTolerance);

        final PosteriorSummary outlierProbabilityPosteriorSummary = globalParameterPosteriorSummaries.get(AlleleFractionParameter.OUTLIER_PROBABILITY);
        final double outlierProbabilityPosteriorCenter = outlierProbabilityPosteriorSummary.getCenter();
        Assert.assertEquals(outlierProbabilityPosteriorCenter, outlierProbability, outlierProbabilityTolerance);

        final List<PosteriorSummary> minorAlleleFractionPosteriorSummaries =
                modeller.getMinorAlleleFractionsPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        final List<Double> minorFractionsPosteriorCenters = minorAlleleFractionPosteriorSummaries.stream().map(PosteriorSummary::getCenter).collect(Collectors.toList());
        double totalPosteriorCentersSegmentError = 0.0;
        for (int segment = 0; segment < numSegments; segment++) {
            totalPosteriorCentersSegmentError += Math.abs(minorFractionsPosteriorCenters.get(segment) - simulatedData.getTrueState().segmentMinorFraction(segment));
        }
        Assert.assertEquals(totalPosteriorCentersSegmentError / numSegments, 0.0, minorFractionTolerance);
    }

    @DataProvider(name = "biasCorrection")
    public Object[][] dataBiasCorrection() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final AllelicCountCollection sampleNormal = new AllelicCountCollection(SAMPLE_NORMAL_FILE);
        final AllelicCountCollection sampleWithBadSNPs = new AllelicCountCollection(SAMPLE_WITH_BAD_SNPS_FILE);
        final AllelicCountCollection sampleWithEvent = new AllelicCountCollection(SAMPLE_WITH_EVENT_FILE);

        final AllelicPanelOfNormals allelicPONNormal = new AllelicPanelOfNormals(ALLELIC_PON_NORMAL_FILE);
        final AllelicPanelOfNormals allelicPONWithBadSNPs = new AllelicPanelOfNormals(ALLELIC_PON_WITH_BAD_SNPS_FILE);

        final double minorFractionExpectedInMiddleSegmentNormal = 0.5;
        final double minorFractionExpectedInMiddleSegmentWithBadSNPsAndNormalPON = 0.4;
        final double minorFractionExpectedInMiddleSegmentWithEvent = 0.33;

        return new Object[][]{
                {sampleNormal, allelicPONNormal, minorFractionExpectedInMiddleSegmentNormal},
                {sampleWithBadSNPs, allelicPONNormal, minorFractionExpectedInMiddleSegmentWithBadSNPsAndNormalPON},
                {sampleWithEvent, allelicPONNormal, minorFractionExpectedInMiddleSegmentWithEvent},
                {sampleWithBadSNPs, allelicPONWithBadSNPs, minorFractionExpectedInMiddleSegmentNormal}
        };
    }

    /**
     * Tests that the allelic PoN is appropriately used to correct reference bias.  The basic set up for the test data is
     * simulated hets at 1000 sites (1:1-1000) across 3 segments.  The outer two segments are balanced with
     * minor-allele fraction = 0.5; however, in the middle segment consisting of 100 sites (1:451-550), all of the sites
     *
     * <p>
     *     1) are balanced and have biases identical to the sites in the other two segments,
     *     which are drawn from a gamma distribution with alpha = 65, beta = 60 -> mean bias = 1.083 ("SAMPLE_NORMAL")
     * </p>
     *
     * <p>
     *     2) are balanced and have relatively high biases,
     *     which are drawn from a gamma distribution with alpha = 9, beta = 6 -> mean bias = 1.5 ("SAMPLE_WITH_BAD_SNPS")
     * </p>
     *
     * <p>
     *     3) have minor-allele fraction = 0.33, copy ratio = 1.5, and biases identical to the sites in the other two segments,
     *     which are drawn from a gamma distribution with alpha = 65, beta = 60 -> mean bias = 1.083 ("SAMPLE_EVENT").
     * </p>
     *
     * In this segment, using a PON that doesn't know about the high reference bias of these sites ("ALLELIC_PON_NORMAL")
     * we should infer a minor-allele fraction of 6 / (6 + 9) = 0.40 in scenario 2; however, with a PON that does know
     * about the high bias at these sites ("ALLELIC_PON_WITH_BAD_SNPS") we correctly infer that all of the segments are balanced.
     *
     * <p>
     *     Note that alpha and beta are not actually correctly recovered in this PON via MLE because the biases are
     *     drawn from a mixture of gamma distributions (as opposed to a single gamma distribution as assumed in the model).
     *     TODO https://github.com/broadinstitute/gatk-protected/issues/421
     * </p>
     */
    @Test(dataProvider = "biasCorrection")
    public void testBiasCorrection(final AllelicCountCollection sample,
                                   final AllelicPanelOfNormals allelicPON,
                                   final double minorFractionExpectedInMiddleSegment) {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final double minorFractionTolerance = 0.025;

        final Genome genome = new Genome(AlleleFractionSimulatedData.TRIVIAL_TARGETS, sample.getCounts());
        final List<SimpleInterval> segments = SegmentUtils.readIntervalsFromSegmentFile(SEGMENTS_FILE);
        final SegmentedGenome segmentedGenome = new SegmentedGenome(segments, genome);

        final int numSamples = 100;
        final int numBurnIn = 25;
        final AlleleFractionModeller modeller = new AlleleFractionModeller(segmentedGenome, allelicPON);
        modeller.fitMCMC(numSamples, numBurnIn);

        final List<PosteriorSummary> minorAlleleFractionPosteriorSummaries =
                modeller.getMinorAlleleFractionsPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        final List<Double> minorFractionsResult = minorAlleleFractionPosteriorSummaries.stream().map(PosteriorSummary::getCenter).collect(Collectors.toList());

        final double minorFractionBalanced = 0.5;
        final List<Double> minorFractionsExpected = Arrays.asList(minorFractionBalanced, minorFractionExpectedInMiddleSegment, minorFractionBalanced);
        for (int segment = 0; segment < 3; segment++) {
            Assert.assertEquals(minorFractionsResult.get(segment), minorFractionsExpected.get(segment), minorFractionTolerance);
        }
    }
}