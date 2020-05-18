package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public abstract class AS_StrandBiasUnitTest extends ReducibleAnnotationBaseTest {

    @DataProvider
    public Object[][] getAlleleSpecificAnnotationMaps() {
        List<Object[]> tests = new ArrayList<>();

        //well-behaved data
        final Map<Allele, Double> alleleSpecificAnnotations = new LinkedHashMap<>();
        alleleSpecificAnnotations.put(Allele.REF_A, Double.NaN);
        alleleSpecificAnnotations.put(Allele.ALT_C, 1.1);
        alleleSpecificAnnotations.put(Allele.NON_REF_ALLELE, 2.2);
        final VariantContext goodVC = new VariantContextBuilder().chr("chr1").start(10001).stop(10001).alleles(Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.NON_REF_ALLELE)).make();
        tests.add(new Object[]{goodVC, alleleSpecificAnnotations});

        final VariantContext bigVC = new VariantContextBuilder().chr("chr1").start(10001).stop(10001).alleles(Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_T, Allele.NON_REF_ALLELE)).make();
        tests.add(new Object[]{bigVC, alleleSpecificAnnotations});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getAlleleSpecificAnnotationMaps")
    public void testReducedAnnotationString(final VariantContext vc, final Map<Allele,Double> perAltsStrandCounts) {
        final String annotationString = AS_StrandBiasTest.makeReducedAnnotationString(vc, perAltsStrandCounts);
        Assert.assertTrue(StringUtils.splitPreserveAllTokens(annotationString, ",").length == vc.getAlternateAlleles().size());
        Assert.assertTrue(StringUtils.countMatches(annotationString, AnnotationUtils.ALLELE_SPECIFIC_REDUCED_DELIM) == vc.getAlternateAlleles().size()-1);  //make sure there are the right number of commas
        Assert.assertFalse(StringUtils.contains(annotationString, AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM));
        final List<Allele> remainingAlleles = new ArrayList<>(vc.getAlternateAlleles());
        remainingAlleles.removeAll(perAltsStrandCounts.keySet());
        if (!remainingAlleles.isEmpty()) {
            Assert.assertTrue(StringUtils.contains(annotationString, "."));
        }
    }
}
