package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class GenotypingGivenAllelesUtilsUnitTest {

    @Test
    public void testComposeGivenAllelesVariantContextFromRod() {

        final SimpleInterval loc = new SimpleInterval("20", 10093568, 10093568);

        final List<Allele> sameLocDelAlleles1 = Arrays.asList(Allele.create("GT", true), Allele.create("G"));
        final VariantContext sameLocDelVc1 = new VariantContextBuilder("a", "20", 10093568, 10093569, sameLocDelAlleles1).make();

        final List<Allele> sameLocSnpAlleles1 = Arrays.asList(Allele.create("G", true), Allele.create("C"));
        final VariantContext sameLocSnpVc1 = new VariantContextBuilder("a", "20", 10093568, 10093568, sameLocSnpAlleles1).make();

        final List<Allele> spanningDelAlleles1 = Arrays.asList(Allele.create("ATGTA", true), Allele.create("A"));
        final VariantContext spanningDelVc1 = new VariantContextBuilder("a", "20", 10093566, 10093570, spanningDelAlleles1).make();

        final List<Allele> expectedAlleles = Arrays.asList(Allele.create("GT", true), Allele.create("G"), Allele.create("CT"));
        final VariantContext expectedVC = new VariantContextBuilder("a", "20", 10093568, 10093569, expectedAlleles).make();

        final VariantContext variantContext =
                GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromVariantList(Arrays.asList(spanningDelVc1, sameLocDelVc1, sameLocSnpVc1),
                        loc,
                        false, true);

        VariantContextTestUtils.assertVariantContextsAreEqual(variantContext, expectedVC, Collections.emptyList());
    }
}