package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class AmbiguousBaseReadFilterTest {

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][]{
                {"CGTTTTTCAGTGATTTCTTCATTTTTCAATTCGTCAAGTGGATGTTTCTCATTTTCCATGATTTTCAGTTTTCTTGCCATATTCCACGTCCTACAGTGGA", Boolean.TRUE},
                {"NNNNNNTCAGTGATTTCTTCATTTTTCAATTCGTCAAGTGGATGTTTCTCATTTTCCATGATTTTCAGTTTTCTTGCCATATTCCACGTCCTACAGTGGA", Boolean.FALSE},
                {"CGTTTTTCAGTGATTTCTTCATTTTTCAATTCGTCAAGTGGATGTTTCTCATTTTCCATGATTTTNNNNNNTCTTGCCATATTCCACGTCCTACAGTGGA", Boolean.FALSE},
                {"CGTTTTTCAGTGATTTCTTCATTNNNNNATTCGTCAAGTGGATGTTTCTCATTTTCCATGATTTTCAGTTTTCTTGCCATATTCCACGTCCTACAGTGGA", Boolean.TRUE}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testMaxBases(final String seq_in, final Boolean test_out) throws Exception {
        AmbiguousBaseReadFilter filter = new AmbiguousBaseReadFilter(5);
        final byte[] qual = new byte[seq_in.length()];
        Arrays.fill(qual, (byte) 30);
        GATKRead read_in = ArtificialReadUtils.createArtificialRead(seq_in.getBytes(), qual, "*");
        final boolean test_i = filter.test(read_in);
        Assert.assertEquals(test_out.booleanValue(), test_i);
    }

    @Test(dataProvider = "sequenceStrings")
    public void testMaxFraction(final String seq_in, final Boolean test_out) throws Exception {
        AmbiguousBaseReadFilter filter = new AmbiguousBaseReadFilter();
        filter.maxAmbiguousBaseFraction = 0.05;
        final byte[] qual = new byte[seq_in.length()];
        Arrays.fill(qual, (byte) 30);
        GATKRead read_in = ArtificialReadUtils.createArtificialRead(seq_in.getBytes(), qual, "*");
        final boolean test_i = filter.test(read_in);
        Assert.assertEquals(test_out.booleanValue(), test_i);
    }

}