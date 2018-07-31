package org.broadinstitute.hellbender.tools.walkers.annotator;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class StrandBiasBySampleUnitTest extends GATKBaseTest {

    @Test
    public void testGetContingencyArray() throws Exception {
        final int[][] t = new int[2][2];
        t[0][0] = 1; t[0][1] = 2; t[1][0] = 3; t[1][1] = 4;
        final List<Integer> tList = StrandBiasBySample.getContingencyArray(t);
        final List<Integer> truthList = new ArrayList<>();
        truthList.add(1); truthList.add(2); truthList.add(3); truthList.add(4);
        Assert.assertEquals(tList, truthList);
    }
}
