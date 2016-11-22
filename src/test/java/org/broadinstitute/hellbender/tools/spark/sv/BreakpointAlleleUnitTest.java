package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.Strandedness.THREE_TO_FIVE;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.Strandedness.FIVE_TO_THREE;

public class BreakpointAlleleUnitTest extends BaseTest{

    private static final PipelineOptions dummyOptions = null;
    private static final SAMSequenceDictionary seqDict = new ReferenceMultiSource(dummyOptions, b37_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION).getReferenceSequenceDictionary(null);

    // -----------------------------------------------------------------------------------------------
    // Tests for generic functions on the base class
    // -----------------------------------------------------------------------------------------------

    @Test
    public void testStrandedness() {

        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 10000, 10100), TextCigarCodec.decode("100M"), true, 60, 0, 1, 100);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 20100, 20200), TextCigarCodec.decode("100M"), false, 60, 0, 101, 200);
        final ChimericAlignment breakpoint1 = new ChimericAlignment(region1, region2, "", "", new ArrayList<>());

        Assert.assertTrue(breakpoint1.involvesStrandSwitch());
        Assert.assertNotEquals(new BreakpointAllele(breakpoint1, seqDict).determineStrandedness(), BreakpointAllele.Strandedness.SAME_STRAND);

        final AlignmentRegion region3 = new AlignmentRegion("4", "contig-7", new SimpleInterval("21", 38343346, 38343483), TextCigarCodec.decode("137M141S"), true, 60, 0, 1, 137);
        final AlignmentRegion region4 = new AlignmentRegion("4", "contig-7", new SimpleInterval("20", 38342908, 38343049), TextCigarCodec.decode("137S141M"), false, 60, 0, 138, 278);
        final ChimericAlignment breakpoint2 = new ChimericAlignment(region3, region4, "", "", new ArrayList<>());

        Assert.assertTrue(breakpoint2.involvesStrandSwitch());
        Assert.assertEquals(new BreakpointAllele(breakpoint2, seqDict).determineStrandedness(), BreakpointAllele.Strandedness.FIVE_TO_THREE);

        final AlignmentRegion region5 = new AlignmentRegion("3", "contig-7", new SimpleInterval("21", 38343346, 38343483), TextCigarCodec.decode("137M141S"), true, 60, 0, 1, 137);
        final AlignmentRegion region6 = new AlignmentRegion("3", "contig-7", new SimpleInterval("21", 38342908, 38343049), TextCigarCodec.decode("137S141M"), false, 60, 0, 138, 278);
        final ChimericAlignment breakpoint3 = new ChimericAlignment(region5, region6, "", "", new ArrayList<>());

        Assert.assertTrue(breakpoint3.involvesStrandSwitch());
        Assert.assertEquals(new BreakpointAllele(breakpoint3, seqDict).determineStrandedness(), BreakpointAllele.Strandedness.FIVE_TO_THREE);
    }

    @Test
    public void testEqualsAndHashCode() throws Exception {

        final BreakpointAllele breakpointAllele1 = getTestBreakpointAllele("1", "contig-1", "foo");

        final BreakpointAllele breakpointAllele2 = getTestBreakpointAllele("2", "contig-2", "bar");

        Assert.assertEquals(breakpointAllele1, breakpointAllele2);
        Assert.assertEquals(breakpointAllele1.hashCode(), breakpointAllele2.hashCode());
    }

    @Test(groups = "sv")
    void testKryoSerializer() {
        // uses inversion subclass for testing
        final BreakpointAllele breakpointAllele1 = getTestBreakpointAllele("1", "contig-1", "foo");
        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, breakpointAllele1);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final BreakpointAllele roundTrip = (BreakpointAllele)kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, breakpointAllele1);
    }

    private static BreakpointAllele getTestBreakpointAllele(final String assemblyId, final String contigId, final String insertionMapping) {
        final AlignmentRegion region1 = new AlignmentRegion(assemblyId, contigId, new SimpleInterval("20", 10000, 10100), TextCigarCodec.decode("100M"), true, 60, 0, 1, 100);
        final AlignmentRegion region2 = new AlignmentRegion(assemblyId, contigId, new SimpleInterval("20", 20100, 20200), TextCigarCodec.decode("100M"), false, 60, 0, 101, 200);
        final ArrayList<String> insertionMappings = new ArrayList<>();
        insertionMappings.add(insertionMapping);
        final ChimericAlignment breakpoint = new ChimericAlignment(region1, region2, "ACAC", "TGTGTGT", insertionMappings);
        return new BreakpointAllele(breakpoint, seqDict);
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for inversion
    // -----------------------------------------------------------------------------------------------

    @Test
    public void test5to3InversionCtor_1() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", new SimpleInterval("20", 100, 205), TextCigarCodec.decode("105M100S"), true, 60, 0, 1, 105);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", new SimpleInterval("20", 500, 605), TextCigarCodec.decode("100S105M"), false, 60, 0, 95, 200);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, "ACACA", "", new ArrayList<>());
        final BreakpointAllele breakpointAllele = new BreakpointAllele(chimericAlignment, seqDict);
        Assert.assertEquals(breakpointAllele.leftJustifiedLeftBreakpoint, new SimpleInterval("20", 200, 200));
        Assert.assertEquals(breakpointAllele.leftJustifiedRightBreakpoint, new SimpleInterval("20", 605, 605));
        Assert.assertEquals(breakpointAllele.homology, "ACACA");
        Assert.assertEquals(breakpointAllele.determineStrandedness(), FIVE_TO_THREE);
    }

    @Test
    public void test5to3InversionCtor_2() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", new SimpleInterval("20", 500, 605), TextCigarCodec.decode("105M100S"), true, 60, 0, 1, 105);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", new SimpleInterval("20", 100, 205), TextCigarCodec.decode("100S105M"), false, 60, 0, 95, 200);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, "ACACA", "", new ArrayList<>());
        final BreakpointAllele breakpointAllele = new BreakpointAllele(chimericAlignment, seqDict);
        Assert.assertEquals(breakpointAllele.leftJustifiedLeftBreakpoint, new SimpleInterval("20", 200, 200));
        Assert.assertEquals(breakpointAllele.leftJustifiedRightBreakpoint, new SimpleInterval("20", 605, 605));
        Assert.assertEquals(breakpointAllele.homology, "ACACA");
        Assert.assertEquals(breakpointAllele.determineStrandedness(), FIVE_TO_THREE);
    }

    @Test
    public void test3to5InversionCtor_1() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", new SimpleInterval("20", 200, 305), TextCigarCodec.decode("100S105M"), false, 60, 0, 95, 200);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", new SimpleInterval("20", 600, 705), TextCigarCodec.decode("105M100S"), true, 60, 0, 1, 105);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, "ACACA", "", new ArrayList<>());
        final BreakpointAllele breakpointAllele = new BreakpointAllele(chimericAlignment, seqDict);
        Assert.assertEquals(breakpointAllele.leftJustifiedLeftBreakpoint, new SimpleInterval("20", 200, 200));
        Assert.assertEquals(breakpointAllele.leftJustifiedRightBreakpoint, new SimpleInterval("20", 605, 605));
        Assert.assertEquals(breakpointAllele.homology, "ACACA");
        Assert.assertEquals(breakpointAllele.determineStrandedness(), THREE_TO_FIVE);
    }

    @Test
    public void test3to5InversionCtor_2() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", new SimpleInterval("20", 600, 705), TextCigarCodec.decode("105M100S"), false, 60, 0, 1, 105);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", new SimpleInterval("20", 200, 305), TextCigarCodec.decode("100S105M"), true, 60, 0, 95, 200);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, "ACACA", "", new ArrayList<>());
        final BreakpointAllele breakpointAllele = new BreakpointAllele(chimericAlignment, seqDict);
        Assert.assertEquals(breakpointAllele.leftJustifiedLeftBreakpoint, new SimpleInterval("20", 200, 200));
        Assert.assertEquals(breakpointAllele.leftJustifiedRightBreakpoint, new SimpleInterval("20", 605, 605));
        Assert.assertEquals(breakpointAllele.homology, "ACACA");
        Assert.assertEquals(breakpointAllele.determineStrandedness(), THREE_TO_FIVE);
    }
}