package org.broadinstitute.hellbender.utils.genotyper;

import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Testing the basic functionality of the diploid genotype class
 */
public final class DiploidGenotypeUnitTest extends BaseTest {

    @Test
    public void testCreateDiploidFromString() {
        String genotype = "AA";
        DiploidGenotype g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));
        Assert.assertFalse(g.isHet());
        Assert.assertFalse(g.isHetRef(BaseUtils.Base.A.base));
        Assert.assertFalse(g.isHetRef(BaseUtils.Base.C.base));
        Assert.assertFalse(g.isHetRef(BaseUtils.Base.G.base));
        Assert.assertFalse(g.isHetRef(BaseUtils.Base.T.base));

        Assert.assertTrue(g.isHomRef(BaseUtils.Base.A.base));
        Assert.assertFalse(g.isHomRef(BaseUtils.Base.C.base));
        Assert.assertFalse(g.isHomRef(BaseUtils.Base.G.base));
        Assert.assertFalse(g.isHomRef(BaseUtils.Base.T.base));

        Assert.assertFalse(g.isHomVar(BaseUtils.Base.A.base));
        Assert.assertTrue(g.isHomVar(BaseUtils.Base.C.base));
        Assert.assertTrue(g.isHomVar(BaseUtils.Base.G.base));
        Assert.assertTrue(g.isHomVar(BaseUtils.Base.T.base));

        genotype = "AC";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));
        Assert.assertTrue(g.isHet());

        Assert.assertTrue(g.isHetRef(BaseUtils.Base.A.base));
        Assert.assertTrue(g.isHetRef(BaseUtils.Base.C.base));
        Assert.assertFalse(g.isHetRef(BaseUtils.Base.G.base));
        Assert.assertFalse(g.isHetRef(BaseUtils.Base.T.base));

        Assert.assertFalse(g.isHomRef(BaseUtils.Base.A.base));
        Assert.assertFalse(g.isHomRef(BaseUtils.Base.C.base));
        Assert.assertFalse(g.isHomRef(BaseUtils.Base.G.base));
        Assert.assertFalse(g.isHomRef(BaseUtils.Base.T.base));

        Assert.assertFalse(g.isHomVar(BaseUtils.Base.A.base));
        Assert.assertFalse(g.isHomVar(BaseUtils.Base.C.base));
        Assert.assertFalse(g.isHomVar(BaseUtils.Base.G.base));
        Assert.assertFalse(g.isHomVar(BaseUtils.Base.T.base));

        genotype = "AG";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        g = DiploidGenotype.createDiploidGenotype(BaseUtils.Base.A.base, BaseUtils.Base.G.base);
        Assert.assertTrue(genotype.equals(g.toString()));

        g = DiploidGenotype.createDiploidGenotype(0, 2);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "AT";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "CC";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "CG";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "CT";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "GG";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "GT";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "TT";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testErrorIndex1() throws Exception {
        DiploidGenotype.createDiploidGenotype(-1, BaseUtils.Base.G.base);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testErrorIndex2() throws Exception {
        DiploidGenotype.createDiploidGenotype(BaseUtils.Base.A.base, -1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testErrorBase1() throws Exception {
        DiploidGenotype.createDiploidGenotype((byte)'f', BaseUtils.Base.G.base);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testErrorBase2() throws Exception {
        DiploidGenotype.createDiploidGenotype(BaseUtils.Base.A.base, (byte)'f');
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testErrorHomBase() throws Exception {
        DiploidGenotype.createHomGenotype((byte)'f');
    }

    @Test
    public void testIsHet() {
        String genotype = "AG";
        DiploidGenotype g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(g.isHet());

        genotype = "AA";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(!g.isHet());
    }

     @Test
    public void testIsHom() {
        String genotype = "AA";
        DiploidGenotype g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(g.isHom());

        genotype = "AG";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(!g.isHom());
    }

    @Test
      public void testCreateGenotype() {
        byte ref = 'A';
        DiploidGenotype g = DiploidGenotype.createHomGenotype(ref);
        Assert.assertTrue("AA".equals(g.toString()));

        ref = 'a';
        g = DiploidGenotype.createHomGenotype(ref);
        Assert.assertTrue("AA".equals(g.toString()));

        ref = 't';
        g = DiploidGenotype.createHomGenotype(ref);
        Assert.assertTrue("TT".equals(g.toString()));

        ref = 'T';
        g = DiploidGenotype.createHomGenotype(ref);
        Assert.assertTrue("TT".equals(g.toString()));

    }


}
