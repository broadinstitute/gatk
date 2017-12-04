package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.CombineTenxSVVCF.BreakendAdjacency;
import org.testng.Assert;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public final class CombineTenxSVVCFTest extends GATKBaseTest {

    @Test(groups = "sv")
    public void testParseBreakendAllele() {
        BreakendAdjacency breakendAdjacency = CombineTenxSVVCF.parseBreakendAllele("G]17:198982]");
        Assert.assertEquals(breakendAdjacency.before, false);
        Assert.assertEquals(breakendAdjacency.revComp, true);
        Assert.assertEquals(breakendAdjacency.contig, "17");
        Assert.assertEquals(breakendAdjacency.position, 198982);
        Assert.assertEquals(breakendAdjacency.localBases, "G");
    }
}