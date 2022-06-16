package org.broadinstitute.hellbender.tools.sv;

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
        Assert.assertEquals(instance.calcBAF(new LocusDepth(tig, pos, smpl, 4, 4, 0, 0), refIndex, altIndex), BafEvidence.MISSING_VALUE);
        // just deep enough
        Assert.assertEquals(instance.calcBAF(new LocusDepth(tig, pos, smpl, 5, 5, 0, 0), refIndex, altIndex), .5);
        // too unequal
        Assert.assertEquals(instance.calcBAF(new LocusDepth(tig, pos, smpl, 7, 3, 0, 0), refIndex, altIndex), BafEvidence.MISSING_VALUE);
        // sufficiently equal
        Assert.assertEquals(instance.calcBAF(new LocusDepth(tig, pos, smpl, 6, 4, 0, 0), refIndex, altIndex), .4);
    }
}
