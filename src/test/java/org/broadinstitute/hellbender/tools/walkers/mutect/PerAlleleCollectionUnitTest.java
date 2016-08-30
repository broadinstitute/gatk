package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

/**
 * Created by tsato on 6/21/16.
 */
public class PerAlleleCollectionUnitTest {
    private static final Allele refA = Allele.create("A", true);
    private static final Allele altA = Allele.create("A", false);
    private static final Allele altC = Allele.create("C", false);
    private static final Allele altG = Allele.create("G", false);
    private static final Allele altT = Allele.create("T", false);

    @Test
    public void testSet() throws Exception {
        PerAlleleCollection<Integer> alleleCounts = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        alleleCounts.set(refA, 40);
        alleleCounts.set(altT, 10);
        Assert.assertEquals(alleleCounts.getRef().intValue(), 40);
        Assert.assertEquals(alleleCounts.getAlt(altT).intValue(), 10);
    }

    @Test
    public void testGet() throws Exception {
        PerAlleleCollection<Integer> alleleCounts = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        alleleCounts.set(refA, 40);
        alleleCounts.set(altT, 10);
        Assert.assertEquals(alleleCounts.get(refA).intValue(), 40);
        Assert.assertEquals(alleleCounts.get(altT).intValue(), 10);
    }

    @Test
    public void testGetAltAlleles() throws Exception {
        PerAlleleCollection<Integer> alleleCounts = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        final List<Allele> altAlleles = Arrays.asList(altA, altC, altG, altT);
        altAlleles.forEach(a -> alleleCounts.set(a,3));
        Assert.assertTrue(alleleCounts.getAltAlleles().containsAll(altAlleles));
        Assert.assertFalse(alleleCounts.getAltAlleles().contains(refA));
    }
}