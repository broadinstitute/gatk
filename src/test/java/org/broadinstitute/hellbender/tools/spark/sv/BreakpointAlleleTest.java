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
import java.util.List;

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


    @Test
    public void testEqualsAndHashCode() throws Exception {

        final BreakpointAllele breakpointAllele1 = getTestBreakpointAllele("1", "contig-1", "foo");

        final BreakpointAllele breakpointAllele2 = getTestBreakpointAllele("2", "contig-2", "bar");


        Assert.assertEquals(breakpointAllele1, breakpointAllele2);
        Assert.assertEquals(breakpointAllele1.hashCode(), breakpointAllele2.hashCode());
    }

    private BreakpointAllele getTestBreakpointAllele(final String assemblyId, final String contigId, final String insertionMapping) {
        final AlignmentRegion region1 = new AlignmentRegion(assemblyId, contigId, TextCigarCodec.decode("100M"), true, new SimpleInterval("1", 10000, 10100), 60, 1, 100, 0);
        final AlignmentRegion region2 = new AlignmentRegion(assemblyId, contigId, TextCigarCodec.decode("100M"), false, new SimpleInterval("1", 20100, 20200), 60, 101, 200, 0);
        final ArrayList<String> insertionMappings = new ArrayList<>();
        insertionMappings.add(insertionMapping);
        final BreakpointAlignment breakpoint = new BreakpointAlignment(contigId, region1, region2, "TGTGTGT", "ACAC", insertionMappings);
        return breakpoint.getBreakpointAllele();
    }

    @Test(groups = "spark")
    void serializationTest() {
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
}