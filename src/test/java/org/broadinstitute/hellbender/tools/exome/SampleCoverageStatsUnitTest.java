package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link SampleCoverageStats}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class SampleCoverageStatsUnitTest extends BaseTest {

    @Test
    public void testCreation() {
        final SampleCoverageStats subject = new SampleCoverageStats("SAMPLE_X", 1.2, 12.2);
        Assert.assertEquals(subject.sample, "SAMPLE_X");
        Assert.assertEquals(subject.mean, 1.2);
        Assert.assertEquals(subject.variance, 12.2);
    }

    @Test
    public void testCreationNegativeVariance() {
        final SampleCoverageStats subject = new SampleCoverageStats("SAMPLE_X", 1.2, - 12.2);
        Assert.assertEquals(subject.sample, "SAMPLE_X");
        Assert.assertEquals(subject.mean, 1.2);
        Assert.assertEquals(subject.variance, 12.2);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullTarget() {
        new SampleCoverageStats(null, 1.2, 12.2);
    }

    @Test
    public void testNaNStats() {
        final SampleCoverageStats subject = new SampleCoverageStats("SAMPLE_X", Double.NaN, Double.NaN);
        Assert.assertEquals(subject.sample, "SAMPLE_X");
        Assert.assertTrue(Double.isNaN(subject.mean));
        Assert.assertTrue(Double.isNaN(subject.variance));
    }

    @Test(dataProvider = "fromSumsData")
    public void testFromSums(final String sample, final long size, final double sum, final double squaresSum) {
        final SampleCoverageStats subject = SampleCoverageStats.fromSums(sample, size, sum, squaresSum);
        Assert.assertEquals(subject.sample, sample);
        if (size == 0) {
            Assert.assertTrue(Double.isNaN(subject.mean));
            Assert.assertTrue(Double.isNaN(subject.variance));
        } else if (size == 1) {
            Assert.assertEquals(subject.mean, sum);
            Assert.assertEquals(subject.variance, 0.0);
        } else {
            final double mean = sum / size;
            Assert.assertEquals(subject.mean, mean, 0.000001);
            final double variance = Math.abs((squaresSum / size - Math.pow(mean, 2)) * size) / ((double) size - 1);
            Assert.assertEquals(subject.variance, variance, 0.000001);
        }
    }


    @Test(dataProvider = "fromSumsWrongData", expectedExceptions = IllegalArgumentException.class)
    public void testFromSumsWrong(final String sample, final long size, final double sum, final double squaresSum) {
        SampleCoverageStats.fromSums(sample, size, sum, squaresSum);
    }

    @DataProvider(name = "fromSumsData")
    public Object[][] fromSumsData() {
        final String sample = "SAMPLE_X";
        return new Object[][]{
                {sample, 0L, 0, 0},
                {sample, 1L, 100, 10000},
                { sample, 4L, 604.4D, 140019.4D }
        };
    }

    @DataProvider(name = "fromSumsWrongData")
    public Object[][] fromSumsWrongData() {
        final String sample = "SAMPLE_X";
        return new Object[][]{
                {null, 0L, 0.0, 0.0},
                {sample, -100L, 0, 0},
                {sample, 1L, 10, 10000},
                { sample, 1L, 604.4D, 1400190.4D },
                { sample, 4L, 10, -10000 },
        };
    }
}
