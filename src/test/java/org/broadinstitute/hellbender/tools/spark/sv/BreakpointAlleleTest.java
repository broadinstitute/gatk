package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;

public class BreakpointAlleleTest {

    @Test
    public void testIsInversion() throws Exception {

        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("100M"), true, new SimpleInterval("1", 10000, 10100), 60, 1, 100, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("100M"), false, new SimpleInterval("1", 20100, 20200), 60, 101, 200, 0);
        final BreakpointAlignment breakpoint1 = new BreakpointAlignment("contig-1", region1, region2, "", "", new ArrayList<>());

        Assert.assertTrue(breakpoint1.getBreakpointAllele().isInversion());

        final AlignmentRegion region3 = new AlignmentRegion("4", "contig-7", TextCigarCodec.decode("137M141S"), true, new SimpleInterval("19", 38343346, 38343483), 60, 1, 137, 0);
        final AlignmentRegion region4 = new AlignmentRegion("4", "contig-7", TextCigarCodec.decode("137S141M"), false, new SimpleInterval("10", 38342908, 38343049), 60, 138, 278, 0);
        final BreakpointAlignment breakpoint2 = new BreakpointAlignment("contig-7", region3, region4, "", "", new ArrayList<>());

        Assert.assertFalse(breakpoint2.getBreakpointAllele().isInversion());

        final AlignmentRegion region5 = new AlignmentRegion("3", "contig-7", TextCigarCodec.decode("137M141S"), true, new SimpleInterval("19", 38343346, 38343483), 60, 1, 137, 0);
        final AlignmentRegion region6 = new AlignmentRegion("3", "contig-7", TextCigarCodec.decode("137S141M"), false, new SimpleInterval("19", 38342908, 38343049), 60, 138, 278, 0);
        final BreakpointAlignment breakpoint3 = new BreakpointAlignment("contig-7", region5, region6, "", "", new ArrayList<>());

        Assert.assertTrue(breakpoint3.getBreakpointAllele().isInversion());

    }

}