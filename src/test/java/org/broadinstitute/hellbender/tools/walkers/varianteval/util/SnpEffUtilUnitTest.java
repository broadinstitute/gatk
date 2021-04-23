package org.broadinstitute.hellbender.tools.walkers.varianteval.util;

/**
 * Created with IntelliJ IDEA.
 * User: farjoun
 * Date: 6/5/13
 * Time: 2:31 PM
 * To change this template use File | Settings | File Templates.
 */


import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class SnpEffUtilUnitTest {


    @DataProvider(name="effects")
    public Object[][] childParentpairs() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{SnpEffUtil.EffectType.GENE,SnpEffUtil.EffectType.CHROMOSOME});
        tests.add(new Object[]{SnpEffUtil.EffectType.UTR_3_PRIME,SnpEffUtil.EffectType.TRANSCRIPT});
        tests.add(new Object[]{SnpEffUtil.EffectType.CODON_CHANGE,SnpEffUtil.EffectType.CDS});
        tests.add(new Object[]{SnpEffUtil.EffectType.STOP_GAINED,SnpEffUtil.EffectType.EXON});
        tests.add(new Object[]{SnpEffUtil.EffectType.SYNONYMOUS_START,SnpEffUtil.EffectType.TRANSCRIPT});
        tests.add(new Object[]{SnpEffUtil.EffectType.FRAME_SHIFT,SnpEffUtil.EffectType.CDS});
        tests.add(new Object[]{SnpEffUtil.EffectType.UPSTREAM,SnpEffUtil.EffectType.INTERGENIC});
        tests.add(new Object[]{SnpEffUtil.EffectType.SPLICE_SITE_DONOR,SnpEffUtil.EffectType.INTRON});
        tests.add(new Object[]{SnpEffUtil.EffectType.SPLICE_SITE_ACCEPTOR,SnpEffUtil.EffectType.INTRON});
        tests.add(new Object[]{SnpEffUtil.EffectType.STOP_LOST,SnpEffUtil.EffectType.NON_SYNONYMOUS_CODING});
        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name="self")
    public Object[][] childEqualsParentpairs() {
        List<Object[]> tests = new ArrayList<>();

        for(SnpEffUtil.EffectType type:SnpEffUtil.EffectType.values()){
            tests.add(new Object[]{type,type});
        }
        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name="noneffects")
    public Object[][] nonchildParentpairs() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{SnpEffUtil.EffectType.START_GAINED,SnpEffUtil.EffectType.NON_SYNONYMOUS_CODING});
        tests.add(new Object[]{SnpEffUtil.EffectType.GENE,SnpEffUtil.EffectType.NONE});
        tests.add(new Object[]{SnpEffUtil.EffectType.UTR_3_PRIME,SnpEffUtil.EffectType.CDS});
        tests.add(new Object[]{SnpEffUtil.EffectType.CODON_CHANGE,SnpEffUtil.EffectType.REGULATION});
        tests.add(new Object[]{SnpEffUtil.EffectType.DOWNSTREAM,SnpEffUtil.EffectType.REGULATION});
        tests.add(new Object[]{SnpEffUtil.EffectType.SPLICE_SITE_ACCEPTOR,SnpEffUtil.EffectType.EXON});
        tests.add(new Object[]{SnpEffUtil.EffectType.START_GAINED,SnpEffUtil.EffectType.SYNONYMOUS_START});
        tests.add(new Object[]{SnpEffUtil.EffectType.NON_SYNONYMOUS_CODING,SnpEffUtil.EffectType.DOWNSTREAM});
        tests.add(new Object[]{SnpEffUtil.EffectType.CODON_DELETION,SnpEffUtil.EffectType.INTRON});
        tests.add(new Object[]{SnpEffUtil.EffectType.UTR_5_PRIME,SnpEffUtil.EffectType.EXON_DELETED});
        tests.add(new Object[]{SnpEffUtil.EffectType.INTRON,SnpEffUtil.EffectType.NONE});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "effects")
    public void testSubType(SnpEffUtil.EffectType subType,SnpEffUtil.EffectType parentType) {
        Assert.assertTrue(SnpEffUtil.isSubTypeOf(subType,parentType),String.format("testing that %s is subtype of %s.",subType,parentType));
    }
    @Test(dataProvider = "self")
    public void testSubTypeSelf(SnpEffUtil.EffectType subType,SnpEffUtil.EffectType parentType) {
        Assert.assertTrue(SnpEffUtil.isSubTypeOf(subType,parentType),String.format("testing that %s is subtype of %s.",subType,parentType));
    }
    @Test(dataProvider = "effects")
    public void testNonSubTypeSelf(SnpEffUtil.EffectType parentType,SnpEffUtil.EffectType subType) {
        Assert.assertTrue(!SnpEffUtil.isSubTypeOf(subType,parentType),String.format("testing that %s is subtype of %s.",subType,parentType));
    }
    @Test(dataProvider = "noneffects")
    public void testNonSubType(SnpEffUtil.EffectType subType,SnpEffUtil.EffectType parentType) {
        Assert.assertTrue(!SnpEffUtil.isSubTypeOf(subType, parentType), String.format("testing that %s is NOT subtype of %s.", subType, parentType));
        Assert.assertTrue(!SnpEffUtil.isSubTypeOf(parentType,subType), String.format("testing that %s is NOT subtype of %s.", parentType,subType));
    }
}
