package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class MultiDeBruijnVertexUnitTest {

    @Test
    public void testMergingIdenticals() {
        final MultiDeBruijnVertex v1 = new MultiDeBruijnVertex("fred".getBytes(), true);
        final MultiDeBruijnVertex v2 = new MultiDeBruijnVertex("fred".getBytes(), true);
        Assert.assertEquals(v1, v2);
    }

    @Test
    public void test(){
       final MultiDeBruijnVertex v1 = new MultiDeBruijnVertex("fred".getBytes());
       final MultiDeBruijnVertex v2 = new MultiDeBruijnVertex("fred".getBytes());
       Assert.assertNotEquals(v1, v2);

       Assert.assertEquals(v1.getKmerSize(), 4);
       Assert.assertEquals(v1.getSuffix(), (byte)'d');
       Assert.assertTrue(Arrays.equals(v1.getSequence(), "fred".getBytes()));
       Assert.assertEquals(v1.getSuffixString(), "d");
       Assert.assertTrue(Arrays.equals(v1.getAdditionalSequence(true), "fred".getBytes()));
       Assert.assertTrue(Arrays.equals(v1.getAdditionalSequence(false), "d".getBytes()));

       Assert.assertNotEquals(v1.hashCode(), v2.hashCode());
       Assert.assertEquals(v1.getKmerSize(), v2.getKmerSize());
       Assert.assertEquals(v1.getAdditionalInfo(), v2.getAdditionalInfo());
       Assert.assertTrue(Arrays.equals(v1.getSequence(), v2.getSequence()));
       Assert.assertEquals(v1.hasAmbiguousSequence(), v2.hasAmbiguousSequence());
       Assert.assertNotNull(v1.toString());//not blow up - we dont check the string

        v1.addRead("fred");
        Assert.assertEquals(v1.getAdditionalInfo(), "");
        v1.setAdditionalInfo("some info");
        Assert.assertEquals(v1.getAdditionalInfo(), "some info");
   }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullRead(){
        final MultiDeBruijnVertex v1 = new MultiDeBruijnVertex("fred".getBytes());
        v1.addRead(null);
    }

    @Test
    public void testBasic() {
        final byte[] bases = "ACT".getBytes();
        final MultiDeBruijnVertex v = new MultiDeBruijnVertex(bases);
        Assert.assertEquals(v.getSequence(), bases);
        Assert.assertEquals(v.getSequenceString(), new String(bases));
        Assert.assertEquals(v.length(), bases.length);
        Assert.assertEquals(v.getSuffix(), (byte) 'T');
        Assert.assertEquals(v.getSuffixString(), "T");

        Assert.assertEquals(v.getAdditionalSequence(true), bases);
        Assert.assertEquals(v.getAdditionalSequence(false).length, 1);
        Assert.assertEquals(v.getAdditionalSequence(false)[0], (byte) 'T');
    }
}
