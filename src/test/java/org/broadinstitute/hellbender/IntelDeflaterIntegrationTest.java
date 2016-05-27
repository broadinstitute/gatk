package org.broadinstitute.hellbender;

import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.hellbender.utils.NativeUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

/**
 * Test that it's possible to load libIntelDeflater
 */
public class IntelDeflaterIntegrationTest extends BaseTest {

    // TODO: re-enable once intel-gkl is on maven central and the IntelDeflater
    // TODO: is re-enabled at the htsjdk level
    @Test(enabled = false)
    public void testIntelDeflaterIsAvailable(){
        if ( ! NativeUtils.runningOnLinux() ) {
            throw new SkipException("IntelDeflater not available on this platform");
        }

        if ( NativeUtils.runningOnPPCArchitecture() ) {
            throw new SkipException("IntelDeflater not available for this architecture");
        }

        Assert.assertTrue(DeflaterFactory.usingIntelDeflater(), "libIntelDeflater.so was not loaded. " +
                "This could be due to a configuration error, or your system might not support it.");
    }

}
