package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class BaseVertexUnitTest extends BaseTest {
    @Test
    public void testBasic() {
        final byte[] bases = "ACT".getBytes();
        final BaseVertex v = new BaseVertex(bases);
        Assert.assertEquals(v.getSequence(), bases);
        Assert.assertEquals(v.getAdditionalSequence(false), bases);
        Assert.assertEquals(v.getAdditionalSequence(true), bases);
        Assert.assertEquals(v.getSequenceString(), new String(bases));
        Assert.assertEquals(v.toString(), v.getSequenceString());
        Assert.assertEquals(v.length(), bases.length);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCreationNull() {
        new BaseVertex((byte[])null);
    }

    @Test()
    public void testCreationEmptySeq() {
        final BaseVertex v = new BaseVertex(new byte[0]);
        Assert.assertTrue(v.isEmpty(), "Version with length == 0 should be empty");
    }

    @Test
    public void testEqualsAndHashCode() {
        final BaseVertex v1 = new BaseVertex("ACT".getBytes());
        final BaseVertex v1_eq = new BaseVertex("ACT".getBytes());
        final BaseVertex v2 = new BaseVertex("ACG".getBytes());

        Assert.assertEquals(v1, v1);
        Assert.assertEquals(v1.hashCode(), v1.hashCode());
        Assert.assertEquals(v1, v1_eq);
        Assert.assertEquals(v1.hashCode(), v1_eq.hashCode());
        Assert.assertFalse(v1.equals(v2));
        Assert.assertFalse(v2.equals(v1));
        Assert.assertFalse(v2.hashCode() == v1.hashCode());
        Assert.assertFalse(v2.equals(v1));
    }
}
