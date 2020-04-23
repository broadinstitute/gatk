package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class NuMTFilterToolTest {

    @DataProvider(name = "cutoffData")
    public Object[][] getNuMTFilterTestData() {
        return new Object[][] {
                {10.0, 3, 25},
                {10.0, NuMTFilterTool.DEFAULT_MAX_NUMT_AUTOSOMAL_COPIES, 31},
                {10.0, 5, 37},
                {22.0, 3, 47},
                {22.0, NuMTFilterTool.DEFAULT_MAX_NUMT_AUTOSOMAL_COPIES, 60},
                {22.0, 5, 73},
                {25.0, 3, 52},
                {25.0, NuMTFilterTool.DEFAULT_MAX_NUMT_AUTOSOMAL_COPIES, 67},
                {25.0, 5, 82},
                {27.0, 3, 56},
                {27.0, NuMTFilterTool.DEFAULT_MAX_NUMT_AUTOSOMAL_COPIES, 72},
                {27.0, 5, 87}
        };
    }

    @Test(dataProvider = "cutoffData")
    public void testGetMaxAltDepthCutoff(final double medianAutosomalCoverage, final double maxCopies, final int expected ) {
        Assert.assertEquals(NuMTFilterTool.getMaxAltDepthCutoff(maxCopies, medianAutosomalCoverage), expected);
    }
}