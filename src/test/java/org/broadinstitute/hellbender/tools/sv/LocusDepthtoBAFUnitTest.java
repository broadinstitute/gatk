package org.broadinstitute.hellbender.tools.sv;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.testng.Assert;
import org.testng.annotations.Test;

public class LocusDepthtoBAFUnitTest {
    @Test
    public void testBAFCalc() {
        final LocusDepthtoBAF instance = new LocusDepthtoBAF();
        final String tig = "1";
        final int pos = 1;
        final int refIndex = 0;
        final int altIndex = 1;
        final String smpl = "s";
        // too shallow
        Assert.assertEquals(instance.calcBAF(new LocusDepth(tig, pos, smpl, refIndex, altIndex, 4, 4, 0, 0)), 0.);
        // just deep enough
        Assert.assertEquals(instance.calcBAF(new LocusDepth(tig, pos, smpl, refIndex, altIndex, 5, 5, 0, 0)), .5);
        // too unequal
        Assert.assertEquals(instance.calcBAF(new LocusDepth(tig, pos, smpl, refIndex, altIndex, 7, 3, 0, 0)), 0.);
        // sufficiently equal
        Assert.assertEquals(instance.calcBAF(new LocusDepth(tig, pos, smpl, refIndex, altIndex, 6, 4, 0, 0)), .4);
    }
}
