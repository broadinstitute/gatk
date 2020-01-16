package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.testng.Assert;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class InferOriginalReadEngineUnitTest {

    @Test
    public void testVarianceToIndelQual(){
        Assert.assertEquals(InferOriginalReadEngine.varianceToIndelQuality(30), (byte)15);
    }

    @Test
    public void testComputeVarianceAroundMostCommon(){
    }

}