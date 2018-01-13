package org.broadinstitute.hellbender.tools.copynumber.formats;

import org.testng.Assert;
import org.testng.annotations.Test;

public final class CopyNumberFormatsUtilsUnitTest {
    @Test
    public void testFormatDouble() {
        Assert.assertEquals(CopyNumberFormatsUtils.formatDouble(Double.NaN), "NaN");
        Assert.assertEquals(CopyNumberFormatsUtils.formatDouble(Double.POSITIVE_INFINITY), "Infinity");
        Assert.assertEquals(CopyNumberFormatsUtils.formatDouble(1.0000001), "1.000000");
    }
}