package org.broadinstitute.hellbender.utils.locusiterator;

import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class LIBSDownsamplingInfoUnitTest extends GATKBaseTest {

    @Test
    public void testToLibsInfoNone() throws Exception {
        final DownsamplingMethod none = DownsamplingMethod.NONE;
        final LIBSDownsamplingInfo libsDownsamplingInfo = LIBSDownsamplingInfo.toDownsamplingInfo(none);
        Assert.assertFalse(libsDownsamplingInfo.isPerformDownsampling());
        Assert.assertEquals(libsDownsamplingInfo.getToCoverage(), 0);
    }

    @Test
    public void testToLibsInfoBySampleToCoverage() throws Exception {
        final DownsamplingMethod bySample = new DownsamplingMethod(DownsamplingMethod.DEFAULT_DOWNSAMPLING_TYPE, 100, null);
        final LIBSDownsamplingInfo libsDownsamplingInfo = LIBSDownsamplingInfo.toDownsamplingInfo(bySample);
        Assert.assertTrue(libsDownsamplingInfo.isPerformDownsampling());
        Assert.assertEquals(libsDownsamplingInfo.getToCoverage(), 100);
    }

    @Test
    public void testToLibsInfoBySampleToFraction() throws Exception {
        final DownsamplingMethod bySample = new DownsamplingMethod(DownsamplingMethod.DEFAULT_DOWNSAMPLING_TYPE, null, 0.1);
        final LIBSDownsamplingInfo libsDownsamplingInfo = LIBSDownsamplingInfo.toDownsamplingInfo(bySample);
        Assert.assertFalse(libsDownsamplingInfo.isPerformDownsampling());
        Assert.assertEquals(libsDownsamplingInfo.getToCoverage(), 0);
    }
}
