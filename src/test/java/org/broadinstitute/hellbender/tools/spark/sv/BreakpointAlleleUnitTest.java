package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion.InversionType.INV_3_TO_5;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion.InversionType.INV_5_TO_3;

public class BreakpointAlleleUnitTest {

    // -----------------------------------------------------------------------------------------------
    // Tests for generic functions on the base class
    // -----------------------------------------------------------------------------------------------

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
        final BreakpointAlleleInversion roundTrip = (BreakpointAlleleInversion)kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, breakpointAllele1);
    }

    private static BreakpointAllele getTestBreakpointAllele(final String assemblyId, final String contigId, final String insertionMapping) {
        final AlignmentRegion region1 = new AlignmentRegion(assemblyId, contigId, TextCigarCodec.decode("100M"), true, new SimpleInterval("1", 10000, 10100), 60, 1, 100, 0);
        final AlignmentRegion region2 = new AlignmentRegion(assemblyId, contigId, TextCigarCodec.decode("100M"), false, new SimpleInterval("1", 20100, 20200), 60, 101, 200, 0);
        final ArrayList<String> insertionMappings = new ArrayList<>();
        insertionMappings.add(insertionMapping);
        final ChimericAlignment breakpoint = new ChimericAlignment(region1, region2, "TGTGTGT", "ACAC", insertionMappings);
        return new BreakpointAlleleInversion(breakpoint);
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for inversion
    // -----------------------------------------------------------------------------------------------

    @Test
    public void test5to3InversionCtor_1() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", TextCigarCodec.decode("105M100S"), true, new SimpleInterval("1", 100, 205), 60, 1, 105, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", TextCigarCodec.decode("100S105M"), false, new SimpleInterval("1", 500, 605), 60, 95, 200, 0);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, "", "ACACA", new ArrayList<>());
        final BreakpointAllele breakpointAllele = new BreakpointAlleleInversion(chimericAlignment);
        Assert.assertEquals(breakpointAllele.leftAlignedLeftBreakpoint, new SimpleInterval("1", 200, 200));
        Assert.assertEquals(breakpointAllele.leftAlignedRightBreakpoint, new SimpleInterval("1", 605, 605));
        Assert.assertEquals(breakpointAllele.homology, "ACACA");
        Assert.assertEquals(((BreakpointAlleleInversion)breakpointAllele).getInversionType(), INV_5_TO_3);
    }

    @Test
    public void test5to3InversionCtor_2() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", TextCigarCodec.decode("105M100S"), true, new SimpleInterval("1", 500, 605), 60, 1, 105, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", TextCigarCodec.decode("100S105M"), false, new SimpleInterval("1", 100, 205), 60, 95, 200, 0);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, "", "ACACA", new ArrayList<>());
        final BreakpointAllele breakpointAllele = new BreakpointAlleleInversion(chimericAlignment);
        Assert.assertEquals(breakpointAllele.leftAlignedLeftBreakpoint, new SimpleInterval("1", 200, 200));
        Assert.assertEquals(breakpointAllele.leftAlignedRightBreakpoint, new SimpleInterval("1", 605, 605));
        Assert.assertEquals(breakpointAllele.homology, "ACACA");
        Assert.assertEquals(((BreakpointAlleleInversion)breakpointAllele).getInversionType(), INV_5_TO_3);
    }

    @Test
    public void test3to5InversionCtor_1() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", TextCigarCodec.decode("100S105M"), false, new SimpleInterval("1", 200, 305), 60, 95, 200, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", TextCigarCodec.decode("105M100S"), true, new SimpleInterval("1", 600, 705), 60, 1, 105, 0);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, "", "ACACA", new ArrayList<>());
        final BreakpointAllele breakpointAllele = new BreakpointAlleleInversion(chimericAlignment);
        Assert.assertEquals(breakpointAllele.leftAlignedLeftBreakpoint, new SimpleInterval("1", 200, 200));
        Assert.assertEquals(breakpointAllele.leftAlignedRightBreakpoint, new SimpleInterval("1", 605, 605));
        Assert.assertEquals(breakpointAllele.homology, "ACACA");
        Assert.assertEquals(((BreakpointAlleleInversion)breakpointAllele).getInversionType(), INV_3_TO_5);
    }

    @Test
    public void test3to5InversionCtor_2() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", TextCigarCodec.decode("105M100S"), false, new SimpleInterval("1", 600, 705), 60, 1, 105, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", TextCigarCodec.decode("100S105M"), true, new SimpleInterval("1", 200, 305), 60, 95, 200, 0);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, "", "ACACA", new ArrayList<>());
        final BreakpointAllele breakpointAllele = new BreakpointAlleleInversion(chimericAlignment);
        Assert.assertEquals(breakpointAllele.leftAlignedLeftBreakpoint, new SimpleInterval("1", 200, 200));
        Assert.assertEquals(breakpointAllele.leftAlignedRightBreakpoint, new SimpleInterval("1", 605, 605));
        Assert.assertEquals(breakpointAllele.homology, "ACACA");
        Assert.assertEquals(((BreakpointAlleleInversion)breakpointAllele).getInversionType(), INV_3_TO_5);
    }
}