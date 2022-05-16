package org.broadinstitute.hellbender.tools.sv;

import org.testng.Assert;
import org.testng.annotations.Test;

public class SiteDepthtoBAFUnitTest {
    @Test
    public void testBAFCalc() {
        final SiteDepthtoBAF instance = new SiteDepthtoBAF();
        final String tig = "1";
        final int pos = 1;
        final int refIndex = 0;
        final int altIndex = 1;
        final String smpl = "s";
        // too shallow
        Assert.assertNull(instance.calcBAF(new SiteDepth(tig, pos, smpl, 4, 4, 0, 0), refIndex, altIndex));
        // just deep enough
        Assert.assertEquals(instance.calcBAF(new SiteDepth(tig, pos, smpl, 5, 5, 0, 0), refIndex, altIndex).getValue(), .5);
        // too unequal
        Assert.assertNull(instance.calcBAF(new SiteDepth(tig, pos, smpl, 7, 3, 0, 0), refIndex, altIndex));
        // sufficiently equal
        Assert.assertEquals(instance.calcBAF(new SiteDepth(tig, pos, smpl, 6, 4, 0, 0), refIndex, altIndex).getValue(), .4);
    }
}
