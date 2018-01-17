package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;


/**
 * Created by gauthier on 6/1/18.
 */
public class GnarlyGenotyperUnitTest {
    @Test
    public void summarizePLs() throws Exception {
        final Allele Aref = Allele.create("A", true);
        final Allele C = Allele.create("C");
        final Allele T = Allele.create("T");

        //AA AB BB AC BC CC
        final int[] homRefPLs1 = {0, 30, 60, 90, 200, 210};  //made up values, but all unique
        final int[] homRefPLs2 = {0, 90, 210, 30, 200, 220};
        final int[] homVarPLs1 = {60, 30, 0, 200, 90, 210};
        final int[] homVarPLs2 = {210, 200, 90, 30, 150, 0};
        final int[] ABHetPLs = {46, 0, 45, 90, 91, 200};
        final int[] ACHetPLs = {46, 90, 200, 0, 30, 45};
        final int[] BCHetPLs = {200, 90, 60, 85, 0, 61};

        final Genotype g1 = VariantContextTestUtils.makeG("g1", Aref, Aref, homRefPLs1);
        final Genotype g2 = VariantContextTestUtils.makeG("g2", Aref, Aref, homRefPLs2);
        final Genotype g3 = VariantContextTestUtils.makeG("g3", C, C, homVarPLs1);
        final Genotype g4 = VariantContextTestUtils.makeG("g4", T, T, homVarPLs2);
        final Genotype g5 = VariantContextTestUtils.makeG("g5", Aref, C, ABHetPLs);
        final Genotype g6 = VariantContextTestUtils.makeG("g6", Aref, T, ACHetPLs);
        final Genotype g7 = VariantContextTestUtils.makeG("g7", C, T, BCHetPLs);

        final VariantContext vc = VariantContextTestUtils.makeVC("test", Arrays.asList(Aref, C, T), g1, g2, g3, g4, g5, g6, g7);

        //ABGQ compares het to hom (g2 is homozygous and allele is in g1) and hom to het (g2 heterozygous and contains g1 allele)
        //ALTGQ takes the best over GTs that have GT's alt filtered out -- het-non-ref is the best over all other GTs

        //For AA ABGQ should be best of A/!A = AB, AC (one allele in genotype, one allele not in genotype)
        //For AA ALTGQ = ABGQ (best across genotypes with GT's alt removed)
        //For AB ABGQ should be the best of AA, BB (is homozygous and allele is in GT)
        //For AB ALTGQ should be the best of B filtered out (AA, AC, CC)
        //For BB ABGQ should be the best of AB, BC
        //For BB ALTGQ should be the best of AA, AC, CC
        //For BC ABGQ should be the best of BB, CC
        //For BC ALTGQ should be the best of B filtered out (AA, AC, CC) or C filtered out (AA, AB, BB) -- anything except BC

        checkIntAttribute(g1, vc, "ABGQ", 30);
        checkIntAttribute(g1, vc, "ALTGQ", 30);
        checkIntAttribute(g2, vc, "ABGQ", 30);
        checkIntAttribute(g2, vc, "ALTGQ", 30);
        checkIntAttribute(g3, vc, "ABGQ", 30);
        checkIntAttribute(g3, vc, "ALTGQ", 60);
        checkIntAttribute(g4, vc, "ABGQ", 30);
        checkIntAttribute(g4, vc, "ALTGQ", 90);
        checkIntAttribute(g5, vc, "ABGQ", 45);
        checkIntAttribute(g5, vc, "ALTGQ", 46);
        checkIntAttribute(g6, vc, "ABGQ", 45);
        checkIntAttribute(g6, vc, "ALTGQ", 46);
        checkIntAttribute(g7, vc, "ABGQ", 60);
        checkIntAttribute(g7, vc, "ALTGQ", 60);
    }

    /**
     *
     * @param g Genotype to calculate attributes for
     * @param vc VariantContext to get allele order in PLs
     * @param attribute must be int attribute (unchecked cast)
     * @param expected
     */
    @SuppressWarnings("unchecked")
    private static void checkIntAttribute(final Genotype g, final VariantContext vc, final String attribute, final int expected) {
        final GenotypeBuilder gb1 = new GenotypeBuilder(g);
        GnarlyGenotyper.summarizePLs(gb1, g, vc);
        Assert.assertEquals((int)gb1.make().getAnyAttribute(attribute), expected);
    }

}