package org.broadinstitute.hellbender.tools.exome.detectcoveragedropout;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author lichtens &lt;lichtens@broadinstitute.org&gt;
 */
public final class CoverageDropoutDetectorTest extends BaseTest {


    @DataProvider(name = "randomUnivariateGaussianTargetsLowVariance")
    public Object[][] inputUnivariateGaussianTargetsLowVariance() {
        return getUnivariateGaussianTargets(.02);
    }

    @DataProvider(name = "randomUnivariateGaussianTargetsDropoutLowVariance")
    public Object[][] inputUnivariateGaussianTargetsDropoutLowVariance() {
        return getUnivariateGaussianTargetsWithDropout(.05, .05);
    }

    private Object[][] getUnivariateGaussianTargets(final double sigma) {
        Random rng = new Random(337);
        final RandomGenerator randomGenerator = RandomGeneratorFactory.createRandomGenerator(rng);
        NormalDistribution n = new NormalDistribution(randomGenerator, 1, sigma);
        final int numDataPoints = 10000;
        final int numEventPoints = 2000;

        final List<ReadCountRecord.SingleSampleRecord> targetList = new ArrayList<>();
        for (int i = 0; i < (numDataPoints - numEventPoints); i++){
            targetList.add(new ReadCountRecord.SingleSampleRecord(new Target("arbitrary_name", new SimpleInterval("chr1", 100 + 2*i, 101 + 2 * i)), n.sample()));
        }
        for (int i = (numDataPoints - numEventPoints); i < numDataPoints; i++){
            targetList.add(new ReadCountRecord.SingleSampleRecord(new Target("arbitrary_name", new SimpleInterval("chr1", 100 + 2 * i, 101 + 2 * i)), 0.5 + n.sample()));
        }

        HashedListTargetCollection<ReadCountRecord.SingleSampleRecord> targets = new HashedListTargetCollection<>(targetList);

        List<ModeledSegment> segments = new ArrayList<>();
        segments.add(new ModeledSegment(new SimpleInterval("chr1", 100, 16050), 8000, 1));
        segments.add(new ModeledSegment(new SimpleInterval("chr1", 16100, 20200), 2000, 1.5));

        return new Object [] []{ {targets, segments}};
    }

    private Object[][] getUnivariateGaussianTargetsWithDropout(final double sigma, final double dropoutRate) {
        Random rng = new Random(337);
        final RandomGenerator randomGenerator = RandomGeneratorFactory.createRandomGenerator(rng);
        NormalDistribution n = new NormalDistribution(randomGenerator, 1, sigma);
        final int numDataPoints = 10000;
        final int numEventPoints = 2000;

        // Randomly select dropoutRate of targets and reduce by 25%-75% (uniformly distributed)
        UniformRealDistribution uniformRealDistribution = new UniformRealDistribution(randomGenerator, 0, 1.0);
        final List<ReadCountRecord.SingleSampleRecord> targetList = new ArrayList<>();
        for (int i = 0; i < numDataPoints; i++){
            double coverage = n.sample() + (i < (numDataPoints - numEventPoints) ? 0.0 : 0.5);
            if (uniformRealDistribution.sample() < dropoutRate) {
                double multiplier = .25 + uniformRealDistribution.sample()/2;
                coverage = coverage * multiplier;
            }
            targetList.add(new ReadCountRecord.SingleSampleRecord(new Target("arbitrary_name", new SimpleInterval("chr1", 100 + 2*i, 101 + 2 * i)), coverage));
        }

        HashedListTargetCollection<ReadCountRecord.SingleSampleRecord> targets = new HashedListTargetCollection<>(targetList);

        List<ModeledSegment> segments = new ArrayList<>();
        segments.add(new ModeledSegment(new SimpleInterval("chr1", 100, 16050), 8000, 1));
        segments.add(new ModeledSegment(new SimpleInterval("chr1", 16100, 20200), 2000, 1.5));

        return new Object [] []{ {targets, segments}};
    }

    @Test(dataProvider = "randomUnivariateGaussianTargetsLowVariance")
    public void testNoError(final TargetCollection<ReadCountRecord.SingleSampleRecord> tc, final List<ModeledSegment> segments){
        CoverageDropoutDetector c = new CoverageDropoutDetector();
        Assert.assertFalse(c.determineCoverageDropoutDetected(segments, tc, 0.003, .1, 0.75, .02).isCoverageDropout());
    }

    @Test(dataProvider = "randomUnivariateGaussianTargetsDropoutLowVariance")
    public void testSmallLowVarianceDropout(final TargetCollection<ReadCountRecord.SingleSampleRecord> tc, final List<ModeledSegment> segments){
        CoverageDropoutDetector c = new CoverageDropoutDetector();
        Assert.assertTrue(c.determineCoverageDropoutDetected(segments, tc, .003, .1, 0.75, .02).isCoverageDropout(), "Coverage dropout was detected, when it should not have been (with low variance)");
    }
}
