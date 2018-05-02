package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth.LargeSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.utils.IntrachromosomalBreakpointPair;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;

public class LargeSimpleSVTest {

    @Test(groups = "sv")
    public void testGetters() {
        final List<String> firstContigs = Collections.singletonList("assembly_contig_1");
        final List<String> secondContigs = Collections.singletonList("assembly_contig_2");
        final IntrachromosomalBreakpointPair breakpoints = new IntrachromosomalBreakpointPair(0, 90, 205, firstContigs, secondContigs);
        final LargeSimpleSV sv = new LargeSimpleSV(SimpleSVType.TYPES.DEL, 100, 200, 0, 1, 2, 3, 4, breakpoints);

        Assert.assertEquals(sv.getEventType(), SimpleSVType.TYPES.DEL);
        Assert.assertEquals(sv.getStart(), 100);
        Assert.assertEquals(sv.getEnd(), 200);
        Assert.assertEquals(sv.getContigId(), 0);
        Assert.assertEquals(sv.getReadPairEvidence(), 1);
        Assert.assertEquals(sv.getSplitReadEvidence(), 2);
        Assert.assertEquals(sv.getBreakpoints(), breakpoints);
    }

    @Test(groups = "sv")
    public void testComputeScore() {
        final LargeSimpleSV sv = new LargeSimpleSV(SimpleSVType.TYPES.DEL, 100, 200, 0, 1, 2, 3, 4, null);
        final double expectedScore = (1 + 2) / (double) (3 + 4);
        Assert.assertEquals(sv.getScore(1.0), expectedScore);
        Assert.assertEquals(LargeSimpleSV.computeScore(1, 2, 3, 4, 1.0), expectedScore);

        final LargeSimpleSV sv2 = new LargeSimpleSV(SimpleSVType.TYPES.DEL, 100, 200, 0, 1, 2, 0, 0, null);
        final double expectedScore2 = 3.0;
        Assert.assertEquals(sv2.getScore(1.0), expectedScore2);
        Assert.assertEquals(LargeSimpleSV.computeScore(1, 2, 0, 0, 1.0), expectedScore2);
    }
}