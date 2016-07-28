package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AlignmentRegion;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AssembledBreakpoint;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.ArrayList;


public class CallVariantsFromAlignedContigsSparkTest extends BaseTest {

    @Test
    public void testInversionFilter() throws Exception {

        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("100M"), true, new SimpleInterval("1", 10000, 10100), 60, 1, 100, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("100M"), false, new SimpleInterval("1", 20100, 20200), 60, 101, 200, 0);
        final AssembledBreakpoint breakpoint1 = new AssembledBreakpoint("contig-1", region1, region2, "", "", new ArrayList<>());

        final Tuple2 breakpointTuple1 = new Tuple2<>(breakpoint1.getBreakpointAllele(), new Tuple2<>(new Tuple2<>("1", "contig-1"), breakpoint1));

        Assert.assertTrue(CallVariantsFromAlignedContigsSpark.inversionBreakpointFilter(breakpointTuple1));

        final AlignmentRegion region3 = new AlignmentRegion("4", "contig-7", TextCigarCodec.decode("137M141S"), true, new SimpleInterval("19", 38343346, 38343483), 60, 1, 137, 0);
        final AlignmentRegion region4 = new AlignmentRegion("4", "contig-7", TextCigarCodec.decode("137S141M"), false, new SimpleInterval("10", 38342908, 38343049), 60, 138, 278, 0);
        final AssembledBreakpoint breakpoint2 = new AssembledBreakpoint("contig-7", region3, region4, "", "", new ArrayList<>());

        final Tuple2 breakpointTuple2 = new Tuple2<>(breakpoint2.getBreakpointAllele(), new Tuple2<>(new Tuple2<>("14399","contig-7"), breakpoint2));

        Assert.assertFalse(CallVariantsFromAlignedContigsSpark.inversionBreakpointFilter(breakpointTuple2));

        final AlignmentRegion region5 = new AlignmentRegion("3", "contig-7", TextCigarCodec.decode("137M141S"), true, new SimpleInterval("19", 38343346, 38343483), 60, 1, 137, 0);
        final AlignmentRegion region6 = new AlignmentRegion("3", "contig-7", TextCigarCodec.decode("137S141M"), false, new SimpleInterval("19", 38342908, 38343049), 60, 138, 278, 0);
        final AssembledBreakpoint breakpoint3 = new AssembledBreakpoint("contig-7", region5, region6, "", "", new ArrayList<>());

        final Tuple2 breakpointTuple3 = new Tuple2<>(breakpoint3.getBreakpointAllele(), new Tuple2<>(new Tuple2<>("14399","contig-7"), breakpoint3));

        Assert.assertTrue(CallVariantsFromAlignedContigsSpark.inversionBreakpointFilter(breakpointTuple3));

    }
}