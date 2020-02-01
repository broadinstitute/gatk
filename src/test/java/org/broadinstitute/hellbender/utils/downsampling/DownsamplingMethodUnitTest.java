package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class DownsamplingMethodUnitTest {
    @Test
    public void testNone() throws Exception {
        final DownsamplingMethod none = DownsamplingMethod.NONE;
        Assert.assertEquals(none.toCoverage, null);
        Assert.assertEquals(none.toFraction, null);
        Assert.assertEquals(none.type, DownsampleType.NONE);
        Assert.assertNotNull(none.toString());
    }

    @Test
    public void testToCov() throws Exception {
        final DownsamplingMethod none = new DownsamplingMethod(DownsamplingMethod.DEFAULT_DOWNSAMPLING_TYPE, 100, null);
        Assert.assertEquals((int)none.toCoverage, 100);
        Assert.assertEquals(none.toFraction, null);
        Assert.assertEquals(none.type, DownsampleType.BY_SAMPLE);
        Assert.assertNotNull(none.toString());
    }

    @Test
    public void testToFraction() throws Exception {
        final DownsamplingMethod bySample = new DownsamplingMethod(DownsamplingMethod.DEFAULT_DOWNSAMPLING_TYPE, null, 0.1);
        Assert.assertEquals(bySample.toCoverage, null);
        Assert.assertEquals((double) bySample.toFraction, 0.1);
        Assert.assertEquals(bySample.type, DownsampleType.BY_SAMPLE);
        Assert.assertNotNull(bySample.toString());
    }

    @Test(expectedExceptions = UserException.class)
    public void testInvalid_nullCountAndFraction() throws Exception {
        new DownsamplingMethod(DownsampleType.BY_SAMPLE, null, null);
    }

    @Test(expectedExceptions = UserException.class)
    public void testInvalid_NotNullCountAndFraction() throws Exception {
        new DownsamplingMethod(DownsampleType.BY_SAMPLE, 100, 0.0);
    }

    @Test(expectedExceptions = UserException.class)
    public void testInvalid_NotPositiveCount() throws Exception {
        new DownsamplingMethod(DownsampleType.BY_SAMPLE, 0, null);
    }

    @Test(expectedExceptions = UserException.class)
    public void testInvalid_NegativeFraction() throws Exception {
        new DownsamplingMethod(DownsampleType.BY_SAMPLE, null, -0.1);
    }

    @Test(expectedExceptions = UserException.class)
    public void testInvalid_TooHighFraction() throws Exception {
        new DownsamplingMethod(DownsampleType.BY_SAMPLE, null, 1.1);
    }

}
