package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AlleleSpecificAnnotationDataUnitTest extends GATKBaseTest {

    @Test
    public void testCreate() throws Exception {
        final Allele Aref= Allele.create("A", true);
        final Allele T= Allele.create("T", false);
        final List<Allele> alleles= Arrays.asList(Aref, T);
        String rawData= "1|2";
        final AlleleSpecificAnnotationData<Integer> asad = new AlleleSpecificAnnotationData<>(alleles, rawData);
        Assert.assertEquals(asad.getAlleles(), alleles);
        Assert.assertEquals(asad.getRefAllele(), Aref);
        Assert.assertNull(asad.getAttribute(Aref));
        Assert.assertNull(asad.getAttribute(T));
        Assert.assertEquals(asad.getRawData(), rawData);

        final Map<Allele, Integer> map= new HashMap<>();
        map.put(Aref, 10);
        map.put(T, 11);
        asad.setAttributeMap(map);
        Assert.assertEquals(asad.getAttribute(Aref), (Integer)10);
        Assert.assertEquals(asad.getAttribute(T), (Integer)11);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNoRef() throws Exception {
        final Allele A= Allele.create("A", false);
        final Allele T= Allele.create("T", false);
        final List<Allele> alleles= Arrays.asList(A, T);
        String rawData= "1|2";
        new AlleleSpecificAnnotationData<>(alleles, rawData);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testTwoRef() throws Exception {
        final Allele Aref= Allele.create("A", true);
        final Allele Tref= Allele.create("T", true);
        final List<Allele> alleles= Arrays.asList(Aref, Tref);
        String rawData= "1|2";
        new AlleleSpecificAnnotationData<>(alleles, rawData);
    }

}
