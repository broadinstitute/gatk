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
    public void testCreateVariant() throws Exception {
        //new ReferenceMultiSource(null, b37_reference_20_21)
    }

    @Test
    public void testInversionFilter() throws Exception {

        final AlignmentRegion region1 = new AlignmentRegion(TextCigarCodec.decode("100M"), true, new SimpleInterval("1", 10000, 10100), 60, 1, 100, 0);
        final AlignmentRegion region2 = new AlignmentRegion(TextCigarCodec.decode("100M"), false, new SimpleInterval("1", 20100, 20200), 60, 101, 200, 0);
        final AssembledBreakpoint breakpoint1 = new AssembledBreakpoint("contig-1", region1, region2, "", "", new ArrayList<>());

        final Tuple2<String, AssembledBreakpoint> breakpointTuple1 = new Tuple2<>("1", breakpoint1);

        Assert.assertTrue(CallVariantsFromAlignedContigsSpark.inversionBreakpointFilter(breakpointTuple1));

        final AlignmentRegion region3 = new AlignmentRegion(TextCigarCodec.decode("105M2D17M2D49M131S"), true, new SimpleInterval("20", 48450747, 48450922), 60, 1, 171, 0);
        final AlignmentRegion region4 = new AlignmentRegion(TextCigarCodec.decode("137S31M2I121M8D11M"), false, new SimpleInterval("20", 48450712, 48450883), 60, 138, 302, 0);
        final AssembledBreakpoint breakpoint2 = new AssembledBreakpoint("contig-13", region3, region4, "", "TATATATATATACACAGTATATATATATATATAC", new ArrayList<>());

        final Tuple2<String, AssembledBreakpoint> breakpointTuple2 = new Tuple2<>("14880", breakpoint2);

        Assert.assertFalse(CallVariantsFromAlignedContigsSpark.inversionBreakpointFilter(breakpointTuple2));

        final AlignmentRegion region5 = new AlignmentRegion(TextCigarCodec.decode("137M141S"), true, new SimpleInterval("19", 38343346, 38343483), 60, 1, 137, 0);
        final AlignmentRegion region6 = new AlignmentRegion(TextCigarCodec.decode("137S141M"), false, new SimpleInterval("19", 38342908, 38343049), 60, 138, 278, 0);
        final AssembledBreakpoint breakpoint3 = new AssembledBreakpoint("contig-7", region5, region6, "", "", new ArrayList<>());

        final Tuple2<String, AssembledBreakpoint> breakpointTuple3 = new Tuple2<>("14399", breakpoint3);

        Assert.assertTrue(CallVariantsFromAlignedContigsSpark.inversionBreakpointFilter(breakpointTuple3));

    }
}