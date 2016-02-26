package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.stream.DoubleStream;

/**
 * Unit tests for {@link TargetCoverageStats}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetCoverageStatsUnitTest extends BaseTest {

    @Test
    public void testCreation() {
        final TargetCoverageStats subject = new TargetCoverageStats(new Target("TARGET"), 1.2, 12.2, 3.0);
        Assert.assertEquals(subject.target, new Target("TARGET"));
        Assert.assertEquals(subject.mean, 1.2);
        Assert.assertEquals(subject.variance, 12.2);
        Assert.assertEquals(subject.interquartileRange, 3.0);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullTarget() {
        new TargetCoverageStats(null, 1.2, 12.2, 3.0);
    }

    @Test
    public void testNaNStats() {
        final TargetCoverageStats subject = new TargetCoverageStats(new Target("TARGET"), Double.NaN, Double.NaN, Double.NaN);
        Assert.assertEquals(subject.target, new Target("TARGET"));
        Assert.assertTrue(Double.isNaN(subject.mean));
        Assert.assertTrue(Double.isNaN(subject.variance));
        Assert.assertTrue(Double.isNaN(subject.interquartileRange));
    }

    @Test(dataProvider = "fromCoverageData")
    public void testFromCoverage(final Target target, final double[] coverage) {
        final TargetCoverageStats subject = TargetCoverageStats.fromCoverage(target, coverage);
        Assert.assertEquals(subject.target, target);
        if (coverage.length == 0) {
            Assert.assertTrue(Double.isNaN(subject.mean));
            Assert.assertTrue(Double.isNaN(subject.variance));
            Assert.assertTrue(Double.isNaN(subject.interquartileRange));
        } else if (coverage.length == 1) {
            Assert.assertEquals(subject.mean, coverage[0]);
            Assert.assertEquals(subject.variance, 0.0);
            Assert.assertEquals(subject.interquartileRange, 0.0);
        } else {
            final double mean = DoubleStream.of(coverage).sum() / ((double) coverage.length);
            Assert.assertEquals(subject.mean, mean, 0.000001);
            final double variance = DoubleStream.of(coverage).map(d -> Math.pow(d - mean, 2)).sum()
                    / ((double) coverage.length - 1);
            Assert.assertEquals(subject.variance, variance, 0.000001);
            Assert.assertEquals(subject.interquartileRange, GATKProtectedMathUtils.interquartileRange(coverage), 0.000001);
        }
    }

    @DataProvider(name = "fromCoverageData")
    public Object[][] fromCoverageData() {
        final Target testTarget = new Target("TARGET_X");
        return new Object[][] {
                { testTarget, new double[0] },
                { testTarget, new double[] { - 0.1}},
                { testTarget, new double[] { 100, 200, 300, 4.4 }}
        };
    }
}
