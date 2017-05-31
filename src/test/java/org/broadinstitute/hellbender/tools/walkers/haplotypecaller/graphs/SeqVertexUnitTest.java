package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public final class SeqVertexUnitTest extends BaseTest {
    @Test
    public void testEqualsAndHashCode() {
        final byte[] bases = "ACT".getBytes();
        final SeqVertex v1 = new SeqVertex(bases);
        final SeqVertex v1_neq = new SeqVertex(bases);

        Assert.assertEquals(v1, v1);
        Assert.assertEquals(v1.hashCode(), v1.hashCode());
        Assert.assertFalse(v1.equals(v1_neq));
        Assert.assertFalse(v1_neq.equals(v1));
        Assert.assertFalse(v1_neq.hashCode() == v1.hashCode());
    }

    @DataProvider(name = "WithoutSuffixData")
    public Object[][] makeWithoutSuffixData() {
        List<Object[]> tests = new ArrayList<>();

        final String bases = "ACGTACGTACGT";
        final int l = bases.length();
        for ( int suffixLength = 0; suffixLength <= l; suffixLength++ ) {
            final int suffixStart = l - suffixLength;
            final String prefix = suffixLength == l ? null : bases.substring(0, suffixStart);
            final String suffix = suffixStart == l ? "" : bases.substring(suffixStart, l);
            tests.add(new Object[]{bases, suffix, prefix});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "WithoutSuffixData")
    public void testWithoutSuffix(final String bases, final String suffix, final String expected) {
        final SeqVertex basesSV = new SeqVertex(bases);
        if ( expected == null )
            Assert.assertNull(basesSV.withoutSuffix(suffix.getBytes()), "Failed for bases " + bases + " with suffix " + suffix + " != " + expected);
        else
            Assert.assertEquals(basesSV.withoutSuffix(suffix.getBytes()).getSequenceString(), expected, "Failed for bases " + bases + " with suffix " + suffix + " != " + expected);
    }
}
