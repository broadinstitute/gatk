package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;

public class ChimericAlignmentUnitTest extends BaseTest {
    private static final PipelineOptions dummyOptions = null;
    private static final SAMSequenceDictionary seqDict = new ReferenceMultiSource(dummyOptions, b37_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION).getReferenceSequenceDictionary(null);
    @Test
    public void testGetLeftAlignedLeftBreakpointOnAssembledContig() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", TextCigarCodec.decode("100M100S"), true, new SimpleInterval("20", 100, 200), 60, 1, 100, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", TextCigarCodec.decode("100M100S"), false, new SimpleInterval("20", 500, 600), 60, 101, 200, 0);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, "", "", new ArrayList<>());
        Assert.assertEquals(chimericAlignment.getLeftJustifiedBreakpoints(seqDict)._1(), new SimpleInterval("20", 200, 200));
    }

    @Test
    public void testGetLeftAlignedLeftBreakpointOnAssembledContigWithHomology() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", TextCigarCodec.decode("105M100S"), true, new SimpleInterval("20", 100, 205), 60, 1, 105, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", TextCigarCodec.decode("105M100S"), false, new SimpleInterval("20", 500, 605), 60, 95, 200, 0);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, "", "ACACA", new ArrayList<>());
        Assert.assertEquals(chimericAlignment.getLeftJustifiedBreakpoints(seqDict)._1(), new SimpleInterval("20", 200, 200));
    }

    @Test
    public void testAlignedBreakpointBreakpointAllele() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("146M51S"), true, new SimpleInterval("21", 108569148, 108569294), 60, 1, 146, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("147S50M"), false, new SimpleInterval("21", 108569314, 108569364), 60, 148, 197, 0);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, "TC", "", new ArrayList<>());
        final SimpleInterval leftAlignedLeftBreakpointOnAssembledContig = chimericAlignment.getLeftJustifiedBreakpoints(seqDict)._1();
        Assert.assertEquals(leftAlignedLeftBreakpointOnAssembledContig, new SimpleInterval("21", 108569294, 108569294));
        final SimpleInterval leftAlignedRightBreakpointOnAssembledContig = chimericAlignment.getLeftJustifiedBreakpoints(seqDict)._2();
        Assert.assertEquals(leftAlignedRightBreakpointOnAssembledContig, new SimpleInterval("21", 108569364, 108569364));
    }
}