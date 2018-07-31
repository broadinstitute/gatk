package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public final class ReducibleAnnotationDataUnitTest extends GATKBaseTest {
    @Test(expectedExceptions = GATKException.class)
    public void testValidateAllelesListMultipleReferenceAlleles() throws Exception {
        final Allele Aref= Allele.create("A", true);
        final Allele T= Allele.create("T", false);
        final Allele Gref= Allele.create("G", true);
        String rawData= "1|2";
        final ReducibleAnnotationData<Integer> asad = new ReducibleAnnotationData<>(rawData);
        final Map<Allele, Integer> map= new HashMap<>();
        map.put(Aref, 10);
        map.put(T, 11);
        map.put(Gref, 10);
        asad.setAttributeMap(map);
        asad.validateAllelesList();
    }

    @Test
    public void testCreate() throws Exception {
        final Allele Aref= Allele.create("A", true);
        final Allele T= Allele.create("T", false);
        String rawData= "1|2";
        final ReducibleAnnotationData<Integer> asad = new ReducibleAnnotationData<>(rawData);
        Assert.assertEquals(asad.getAlleles(), Arrays.asList(Allele.NO_CALL));
        Assert.assertNull(asad.getAttribute(Aref));
        Assert.assertNull(asad.getAttribute(T));
        Assert.assertFalse(asad.hasAttribute(Aref));
        Assert.assertFalse(asad.hasAttribute(T));
        Assert.assertEquals(asad.getRawData(), rawData);

        final Map<Allele, Integer> map= new HashMap<>();
        map.put(Aref, 10);
        map.put(T, 11);
        asad.setAttributeMap(map);
        Assert.assertEquals(asad.getAttribute(Aref), (Integer)10);
        Assert.assertEquals(asad.getAttribute(T), (Integer)11);
        Assert.assertEquals(asad.getAttributeMap(), map);

        asad.putAttribute(T, 19);
        Assert.assertEquals(asad.getAttribute(T), (Integer)19);

    }



}
